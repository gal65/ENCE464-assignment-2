#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// Global flag
// set to true when operating in debug mode to enable verbose logging
static bool debug = false;

#define IDX(n, i, j, k) (((k) * (n) + (j)) * (n) + (i))

#define CACHE_LINE_SIZE 64
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

// Define implementation options:
#define PTR_OPTIMIZATION
#define SINGLE_BOUND_LOOP

#define UNSAFE_ASSERT(x)                                                                           \
    if (!(x))                                                                                      \
    __builtin_unreachable()

typedef struct {
    float* source;
    float* curr;
    float* next;
    int n;
    int k_start;
    int k_end;
    int iterations;
    float delta_squared;
} thread_args_t;

typedef struct {
    float* source;
    float* curr;
    float* next;
    int n;
    int iterations;
    float delta_squared;
} boundary_thread_args_t;

// barrier is used for all threads to wait for eachother between iterations
pthread_barrier_t barrier;

void* worker(void* vargs)
{
    thread_args_t* args = (thread_args_t*)vargs;
    float* curr = args->curr;
    float* next = args->next;
    int n = args->n;
#ifdef PTR_OPTIMIZATION
    #define BUF_POINTERS                                                                           \
        X(top)                                                                                     \
        X(bottom)                                                                                  \
        X(front)                                                                                   \
        X(back)                                                                                    \
        X(middle)                                                                                  \
        X(source_ptr)                                                                              \
        X(next_ptr)

    #define X(name) float* name;
    BUF_POINTERS
    #undef X
#endif

    for (int iter = 0; iter < args->iterations; iter++) {
#ifdef PTR_OPTIMIZATION
        middle = &curr[IDX(n, 1, 1, args->k_start)];
        source_ptr = &args->source[IDX(n, 1, 1, args->k_start)];
        next_ptr = &next[IDX(n, 1, 1, args->k_start)];

        // top and bottom start above and below (j=0, 2) the middle pointer
        top = &curr[IDX(n, 1, 0, args->k_start)];
        bottom = &curr[IDX(n, 1, 2, args->k_start)];

        // front and back start into and out of the page (k = start-1, start+1) from the middle
        // pointer
        front = &curr[IDX(n, 1, 1, args->k_start - 1)];
        back = &curr[IDX(n, 1, 1, args->k_start + 1)];
#endif
        for (int k = args->k_start; k < args->k_end; k++) {
            // Tell the compiler that k is not in boundary
            UNSAFE_ASSERT(k > 0);
            UNSAFE_ASSERT(k < n - 1);
            for (int j = 1; j < n - 1; j++) {
                for (int i = 1; i < n - 1; i++) {
#ifdef PTR_OPTIMIZATION
                    float source_term = args->delta_squared * (*source_ptr++);
                    middle++;
                    *next_ptr++ = 1.0 / 6.0 *
                                  ((*top++) + (*bottom++) + (*front++) + (*back++) +
                                   (*(middle - 2)) + (*middle) - source_term);

#endif
#ifdef INDEXING
                    float source_term = args->delta_squared * args->source[IDX(n, i, j, k)];
                    next[IDX(n, i, j, k)] =
                        1.0 / 6.0 *
                        (curr[IDX(n, i + 1, j, k)] + curr[IDX(n, i - 1, j, k)] +
                         curr[IDX(n, i, j + 1, k)] + curr[IDX(n, i, j - 1, k)] +
                         curr[IDX(n, i, j, k + 1)] + curr[IDX(n, i, j, k - 1)] - source_term);
#endif
                    // do_cell(args->source, curr, next, args->delta_squared, n, i, j, k);
                }

#ifdef PTR_OPTIMIZATION
    #define X(name) name += 2;
                BUF_POINTERS
    #undef X
#endif
            }

#ifdef PTR_OPTIMIZATION
    #define X(name) name += 2 * n;
            BUF_POINTERS
    #undef X
#endif
        }

        float* temp = curr;
        curr = next;
        next = temp;

#ifdef PRINT_TIMESTAMP
        struct timeval t;
        gettimeofday(&t, NULL);
        printf("worker wait %ld\n", t.tv_usec);
#endif

        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

/**
 * Do cell update while checking all boundary conditions
 */
// TODO do the pointer optimization like in the normal worker task here as well
inline void
do_cell(float* source, float* curr, float* next, float delta_squared, int n, int i, int j, int k)
{
    int ip = (i == n - 1) ? -1 : 1;
    int in = (i == 0) ? -1 : 1;
    int jp = (j == n - 1) ? -1 : 1;
    int jn = (j == 0) ? -1 : 1;
    int kp = (k == n - 1) ? -1 : 1;
    int kn = (k == 0) ? -1 : 1;

    float source_term = delta_squared * source[IDX(n, i, j, k)];
    next[IDX(n, i, j, k)] = 1.0 / 6.0 *
                            (curr[IDX(n, i + ip, j, k)] + curr[IDX(n, i - in, j, k)] +
                             curr[IDX(n, i, j + jp, k)] + curr[IDX(n, i, j - jn, k)] +
                             curr[IDX(n, i, j, k + kp)] + curr[IDX(n, i, j, k - kn)] - source_term);
}

// TODO at large number of threads or large sizes, boundary thread will surpass worker threads in
// work required. Need to allocate a proportional amount of threads to the boundary so other threads
// are not stuck waiting. OR: each thread does it's own boundary scan
void* boundary_worker(void* vargs)
{
    boundary_thread_args_t* args = (boundary_thread_args_t*)vargs;
    float* curr = args->curr;
    float* next = args->next;
    int n = args->n;

    for (int iter = 0; iter < args->iterations; iter++) {
#ifdef SINGLE_BOUND_LOOP
        for (int b = 0; b < n; b++) {
            for (int a = 0; a < n; a++) {
                do_cell(args->source, curr, next, args->delta_squared, n, a, b, 0);
                do_cell(args->source, curr, next, args->delta_squared, n, a, b, n - 1);
                do_cell(args->source, curr, next, args->delta_squared, n, a, 0, b);
                do_cell(args->source, curr, next, args->delta_squared, n, a, n - 1, b);
                do_cell(args->source, curr, next, args->delta_squared, n, 0, a, b);
                do_cell(args->source, curr, next, args->delta_squared, n, n - 1, a, b);
            }
        }
#endif
#ifdef MANY_BOUND_LOOP
        // top k-slice
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                do_cell(args->source, curr, next, args->delta_squared, n, i, j, 0);
            }
        }
        // bottom k-slice
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                do_cell(args->source, curr, next, args->delta_squared, n, i, j, n - 1);
            }
        }

        // left j-slice
        for (int k = 1; k < n - 1; k++) {
            for (int i = 0; i < n; i++) {
                do_cell(args->source, curr, next, args->delta_squared, n, i, 0, k);
            }
        }
        // right j-slice
        for (int k = 1; k < n - 1; k++) {
            for (int i = 0; i < n; i++) {
                do_cell(args->source, curr, next, args->delta_squared, n, i, n - 1, k);
            }
        }

        // front i-slice
        for (int k = 1; k < n - 1; k++) {
            for (int j = 1; j < n - 1; j++) {
                do_cell(args->source, curr, next, args->delta_squared, n, 0, j, k);
            }
        }
        // back i-slice
        for (int k = 1; k < n - 1; k++) {
            for (int j = 1; j < n - 1; j++) {
                do_cell(args->source, curr, next, args->delta_squared, n, n - 1, j, k);
            }
        }
#endif

        float* temp = curr;
        curr = next;
        next = temp;

#ifdef PRINT_TIMESTAMP
        struct timeval t;
        gettimeofday(&t, NULL);
        printf("boundary wait %ld\n", t.tv_usec);
#endif
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

float* poisson_neumann(int n, float* source, int iterations, int num_threads, float delta)
{
    if (debug) {
        printf(
            "Starting solver with:\n"
            "n = %i\n"
            "iterations = %i\n"
            "num_threads = %i\n"
            "delta = %f\n",
            n,
            iterations,
            num_threads,
            delta);
    }

    float *curr, *next;
    // Allocate curr and next to be cache aligned
    posix_memalign((void**)&curr, CACHE_LINE_SIZE, n * n * n * sizeof(float));
    posix_memalign((void**)&next, CACHE_LINE_SIZE, n * n * n * sizeof(float));

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL) {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    thread_args_t* thread_args = (thread_args_t*)malloc(num_threads * sizeof(thread_args_t));
    float delta_squared = delta * delta;

    // Init the boundary thread
    pthread_barrier_init(&barrier, NULL, num_threads + 1);
    boundary_thread_args_t args = {
        .source = source,
        .curr = curr,
        .next = next,
        .n = n,
        .iterations = iterations,
        .delta_squared = delta_squared};
    pthread_t boundary_thread;
    pthread_create(&boundary_thread, NULL, &boundary_worker, &args);

    // init the worker threads
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
        int block_size = n / num_threads;
        // if the first thread, it's block should start at 1 (not 0) because of boundary
        int k_start = thread_idx == 0 ? 1 : thread_idx * block_size;

        // if the last thread, it's block should go all the way to the end of the array (minus
        // boundary)
        int k_end = thread_idx == num_threads - 1 ? n - 1 : (thread_idx + 1) * block_size;
        // If using boundary thread, the main threads should not touch the boundaries of k
        if (thread_idx == 0) {
            k_start = 1;
        }

        thread_args[thread_idx] = (thread_args_t){
            .source = source,
            .curr = curr,
            .next = next,
            .n = n,
            .k_start = k_start,
            .k_end = k_end,
            .iterations = iterations,
            .delta_squared = delta_squared};

        pthread_create(&threads[thread_idx], NULL, &worker, &thread_args[thread_idx]);
    }

    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
        pthread_join(threads[thread_idx], NULL);
    }
    pthread_join(boundary_thread, NULL);

    if (iterations % 2 != 0) {
        float* temp = curr;
        curr = next;
        next = temp;
    }

    free(next);
    free(threads);
    free(thread_args);

    if (debug) {
        printf("Finished solving.\n");
    }

    return curr;
}

int main(int argc, char** argv)
{

    // default settings for solver
    int iterations = 10;
    int n = 5;
    int threads = 1;
    float delta = 1;

    // parse the command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("usage: poisson [-n size] [-i iterations] [-t threads] "
                   "[--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp(argv[i], "-n") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-i") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-t") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "--debug") == 0) {
            debug = true;
        }
    }

    // ensure we have an odd sized cube
    if (n % 2 == 0) {
        fprintf(stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    float* source;
    posix_memalign((void**)&source, CACHE_LINE_SIZE, n * n * n * sizeof(float));
    memset(source, 0, n * n * n * sizeof(float));

    source[(n * n * n) / 2] = 1;
// #define TIME_RUN
#ifdef TIME_RUN
    struct timeval start, end;
    gettimeofday(&start, NULL);
#endif
    // Calculate the resulting field with Neumann conditions
    float* result = poisson_neumann(n, source, iterations, threads, delta);
#ifdef TIME_RUN
    gettimeofday(&end, NULL);
    int elapsed = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    printf("%d elapsed\n", elapsed);
#else

    // Print out the middle slice of the cube for validation
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            printf("%0.5f ", result[((n / 2) * n + y) * n + x]);
        }
        printf("\n");
    }
#endif

// #define SECOND_SLICE
#ifdef SECOND_SLICE
    printf("\n");
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            printf("%0.5f ", result[((n / 2 + 1) * n + y) * n + x]);
        }
        printf("\n");
    }
#endif

    free(source);
    free(result);

    return EXIT_SUCCESS;
}
