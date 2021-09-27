#include <pthread.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Global flag
// set to true when operating in debug mode to enable verbose logging
static bool debug = false;

#define IDX(n, i, j, k) (((k) * (n) + (j)) * (n) + (i))

#define BOUNDARY_THREAD

#ifdef BOUNDARY_THREAD
// start of iteration when only iterating over non-boundary cells
#define INNER_START (1)
// end of iteration when only iterating over non-boundary cells
#define INNER_END (n - 1)
#define LOOKAHEAD(x) 1
#define LOOKBEHIND(x) -1
#else
// Not using boundary thread
#define INNER_START (0)
#define INNER_END (n)
#define LOOKAHEAD(x) ((x == n - 1) ? -1 : 1)
#define LOOKBEHIND(x) ((x == 0) ? 1 : -1)
#endif

#define BLOCK_SIZE 8
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

typedef struct {
    double* source;
    double* curr;
    double* next;
    int n;
    int k_start;
    int k_end;
    int iterations;
    double delta_squared;
} thread_args_t;

typedef struct {
    double* source;
    double* curr;
    double* next;
    int n;
    int iterations;
    double delta_squared;
} boundary_thread_args_t;

pthread_barrier_t barrier;

void* worker(void* vargs)
{
    thread_args_t* args = (thread_args_t*)vargs;
    double* curr = args->curr;
    double* next = args->next;
    int n = args->n;

    // top/bottom: +-j, front/back: +-k
    double *top, *bottom, *front, *back, *middle, *source_ptr, *next_ptr;

    for (int iter = 0; iter < args->iterations; iter++) {

        // reset pointers to their staring position
        top = &curr[IDX(n, 1, 0, args->k_start)];
        bottom = &curr[IDX(n, 1, 2, args->k_start)];

        front = &curr[IDX(n, 1, 1, args->k_start - 1)];
        back = &curr[IDX(n, 1, 1, args->k_start + 1)];

        middle = &curr[IDX(n, 1, 1, args->k_start)];
        source_ptr = &args->source[IDX(n, 1, 1, args->k_start)];
        next_ptr = &next[IDX(n, 1, 1, args->k_start)];

        for (int k = args->k_start; k < args->k_end; k++) {
            // loop tiling: iterate over j in blocks
            for (int j_block = 1; j_block < INNER_END; j_block += BLOCK_SIZE) {
                int j_next = MIN(n - 1, j_block + BLOCK_SIZE);
                for (int j = j_block; j < j_next; j++) {
                    // note that i,j,k in this worker thread are always from 1 to n-1 (inclusive).
                    // ie. they never touch the boundary.
                    for (int i = 1; i < n - 1; i++) {
                        double source_term = args->delta_squared * (*source_ptr++);
                        middle++;
                        // REMEMBER: "middle" is already pointing to the right neighbor when the
                        // relaxation is calulated.
                        *next_ptr++ = 1.0 / 6.0 *
                                      ((*top++) + (*bottom++) + (*front++) + (*back++) +
                                       (*(middle - 2)) + (*middle) - source_term);
                    }
                    // Since the boundary thread is handling the boundaries (one column on each
                    // side), skip these by incrementing by 2
                    top += 2;
                    bottom += 2;
                    front += 2;
                    back += 2;
                    middle += 2;
                    source_ptr += 2;
                    next_ptr += 2;
                }
            }
            // boundary thread is handling one slice on top and bottom (k == 0 and k == n-1), so
            // skip over these
            top += 2 * n;
            bottom += 2 * n;
            front += 2 * n;
            back += 2 * n;
            middle += 2 * n;
            source_ptr += 2 * n;
            next_ptr += 2 * n;
        }

        double* temp = curr;
        curr = next;
        next = temp;
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

/**
 * Do cell update while checking all boundary conditions
 */
// TODO do the pointer optimization like in the normal worker task here as well:
inline void do_cell(
    double* source, double* curr, double* next, double delta_squared, int n, int i, int j, int k)
{
    int ip = (i == n - 1) ? -1 : 1;
    int in = (i == 0) ? -1 : 1;
    int jp = (j == n - 1) ? -1 : 1;
    int jn = (j == 0) ? -1 : 1;
    int kp = (k == n - 1) ? -1 : 1;
    int kn = (k == 0) ? -1 : 1;

    double source_term = delta_squared * source[IDX(n, i, j, k)];
    next[IDX(n, i, j, k)] = 1.0 / 6.0 *
                            (curr[IDX(n, i + ip, j, k)] + curr[IDX(n, i - in, j, k)] +
                             curr[IDX(n, i, j + jp, k)] + curr[IDX(n, i, j - jn, k)] +
                             curr[IDX(n, i, j, k + kp)] + curr[IDX(n, i, j, k - kn)] - source_term);
}

void* boundary_worker(void* vargs)
{
    boundary_thread_args_t* args = (boundary_thread_args_t*)vargs;
    double* curr = args->curr;
    double* next = args->next;
    int n = args->n;

    for (int iter = 0; iter < args->iterations; iter++) {
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
        double* temp = curr;
        curr = next;
        next = temp;
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

double* poisson_neumann(int n, double* source, int iterations, int num_threads, float delta)
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

    // Allocate some buffers to calculate the solution in
    double* curr = (double*)calloc(n * n * n, sizeof(double));
    double* next = (double*)calloc(n * n * n, sizeof(double));

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL) {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    thread_args_t* thread_args = (thread_args_t*)malloc(num_threads * sizeof(thread_args_t));
    double delta_squared = delta * delta;

#ifdef BOUNDARY_THREAD
    // Init the boundary thread if enabled
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
#else
    pthread_barrier_init(&barrier, NULL, num_threads);
#endif

    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
        int block_size = n / num_threads;
        int k_start = thread_idx * block_size;

        // if the last thread, it's block should go all the way to the end of the array
        int k_end = thread_idx == num_threads - 1 ? n : (thread_idx + 1) * block_size;
#ifdef BOUNDARY_THREAD
        // If using boundary thread, the main threads should not touch the boundaries of k
        if (thread_idx == 0) {
            k_start = 1;
        }
        if (thread_idx == num_threads - 1) {
            k_end = n - 1;
        }
#endif

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

#ifdef BOUNDARY_THREAD
    pthread_join(boundary_thread, NULL);
#endif

    if (iterations % 2 != 0) {
        double* temp = curr;
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
    double* source = (double*)calloc(n * n * n, sizeof(double));
    if (source == NULL) {
        fprintf(stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    // Calculate the resulting field with Neumann conditions
    double* result = poisson_neumann(n, source, iterations, threads, delta);

    // Print out the middle slice of the cube for validation
    for (int x = 0; x < n; ++x) {
        for (int y = 0; y < n; ++y) {
            printf("%0.5f ", result[((n / 2) * n + y) * n + x]);
        }
        printf("\n");
    }

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
