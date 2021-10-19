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
#define MIN(a, b) (((a) < (b)) ? (a) : (b))

#define CACHE_LINE_SIZE 64

#define COLUMN_BLOCK_SIZE 8

typedef float cell_t;

typedef struct {
    cell_t* source;
    cell_t* curr;
    cell_t* next;

    int n;
    // inclusive
    int j_block_start;
    // exclusiive
    int j_block_end;
    int iterations;
    cell_t delta_squared;
} thread_args_t;

// barrier is used for all threads to wait for eachother between iterations
pthread_barrier_t barrier;

void* worker(void* vargs)
{
    thread_args_t* args = (thread_args_t*)vargs;
    cell_t* curr = args->curr;
    cell_t* next = args->next;
    int n = args->n;
    const cell_t delta_squared = args->delta_squared;

    for (int iter = 0; iter < args->iterations; iter++) {
        for (int k = 0; k < n; k++) {
            for (int j = args->j_block_start; j < args->j_block_end; j++) {
                for (int i = 0; i < n; i++) {

                    int ip = (i == n - 1) ? -1 : 1;
                    int in = (i == 0) ? -1 : 1;
                    int jp = (j == n - 1) ? -1 : 1;
                    int jn = (j == 0) ? -1 : 1;
                    int kp = (k == n - 1) ? -1 : 1;
                    int kn = (k == 0) ? -1 : 1;

                    cell_t source_term = delta_squared * args->source[IDX(n, i, j, k)];
                    next[IDX(n, i, j, k)] = 1.0 / 6.0 * (curr[IDX(n, i + ip, j, k)] + curr[IDX(n, i - in, j, k)] + curr[IDX(n, i, j + jp, k)] + curr[IDX(n, i, j - jn, k)] + curr[IDX(n, i, j, k + kp)] + curr[IDX(n, i, j, k - kn)] - source_term);
                }
            }
        }

        cell_t* temp = curr;
        curr = next;
        next = temp;
        pthread_barrier_wait(&barrier);
    }
    return NULL;
}

/**
 * Do cell update while checking all boundary conditions
 */
// TODO do the pointer optimization like in the normal worker task here as well
inline void do_cell(
    cell_t* source, cell_t* curr, cell_t* next, cell_t delta_squared, int n,
    int i, int j, int k)
{
    int ip = (i == n - 1) ? -1 : 1;
    int in = (i == 0) ? -1 : 1;
    int jp = (j == n - 1) ? -1 : 1;
    int jn = (j == 0) ? -1 : 1;
    int kp = (k == n - 1) ? -1 : 1;
    int kn = (k == 0) ? -1 : 1;

    cell_t source_term = delta_squared * source[IDX(n, i, j, k)];
    next[IDX(n, i, j, k)] = 1.0 / 6.0 * (curr[IDX(n, i + ip, j, k)] + curr[IDX(n, i - in, j, k)] + curr[IDX(n, i, j + jp, k)] + curr[IDX(n, i, j - jn, k)] + curr[IDX(n, i, j, k + kp)] + curr[IDX(n, i, j, k - kn)] - source_term);
}

cell_t* poisson_neumann(
    int n, cell_t* source, int iterations, int num_threads, cell_t delta)
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

    cell_t *curr, *next;
    // Allocate curr and next to be cache aligned
    posix_memalign((void**)&curr, CACHE_LINE_SIZE, n * n * n * sizeof(cell_t));
    posix_memalign((void**)&next, CACHE_LINE_SIZE, n * n * n * sizeof(cell_t));

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL) {
        fprintf(
            stderr,
            "Error: ran out of memory when trying to allocate %i sized cube\n",
            n);
        exit(EXIT_FAILURE);
    }

    pthread_t* threads = (pthread_t*)malloc(num_threads * sizeof(pthread_t));
    thread_args_t* thread_args = (thread_args_t*)malloc(num_threads * sizeof(thread_args_t));
    cell_t delta_squared = delta * delta;

    // Init the boundary thread
    pthread_barrier_init(&barrier, NULL, num_threads);

    // init the worker threads
    // printf("num_blocks: %d, per_thread: %d\n", num_blocks,
    // blocks_per_thread);
    int per_thread = n / num_threads;
    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
        int start = thread_idx * per_thread;
        int end = thread_idx == num_threads - 1 ? n : (thread_idx + 1) * per_thread;

        thread_args[thread_idx] = (thread_args_t) { .source = source,
            .curr = curr,
            .next = next,
            .n = n,
            .j_block_start = start,
            .j_block_end = end,
            .iterations = iterations,
            .delta_squared = delta_squared };

        pthread_create(
            &threads[thread_idx],
            NULL,
            &worker,
            &thread_args[thread_idx]);
    }

    for (int thread_idx = 0; thread_idx < num_threads; thread_idx++) {
        pthread_join(threads[thread_idx], NULL);
    }
    // pthread_join(boundary_thread, NULL);

    if (iterations % 2 != 0) {
        cell_t* temp = curr;
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
    cell_t delta = 1;

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

#ifdef CACHE_ALIGN_BUFFERS
    cell_t* source;
    posix_memalign(
        (void**)&source,
        CACHE_LINE_SIZE,
        n * n * n * sizeof(cell_t));
    memset(source, 0, n * n * n * sizeof(cell_t));
#else
    cell_t* source = (cell_t*)calloc(n * n * n, sizeof(cell_t));
#endif

    source[(n * n * n) / 2] = 1;
// #define TIME_RUN
#ifdef TIME_RUN
    struct timeval start, end;
    gettimeofday(&start, NULL);
    cell_t* result = poisson_neumann(n, source, iterations, threads, delta);
    gettimeofday(&end, NULL);
    int elapsed = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_usec - start.tv_usec);
    printf("%d elapsed\n", elapsed);
#else
    cell_t* result = poisson_neumann(n, source, iterations, threads, delta);
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
