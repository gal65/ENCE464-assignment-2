#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * poisson.c
 * Implementation of a Poisson solver with Neumann boundary conditions.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * BUILDING:
 * g++ -o poisson poisson.c -lpthread
 *
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 *
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */

// Global flag
// set to true when operating in debug mode to enable verbose logging
static bool debug = false;

#define INDEX(n, i, j, k) ((k) * (n) * (n) + (j) * (n) + i)

/**
 * @brief Solve Poissons equation for a given cube with Neumann boundary
 * conditions on all sides.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing
 * function.
 * @param iterations    Number of iterations to solve with.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @return double*      Solution to Poissons equation. Caller must free().
 */
double *poisson_neumann(int n, double *source, int iterations, int threads,
                        float delta) {
  if (debug) {
    printf("Starting solver with:\n"
           "n = %i\n"
           "iterations = %i\n"
           "threads = %i\n"
           "delta = %f\n",
           n, iterations, threads, delta);
  }

  // Allocate some buffers to calculate the solution in
  double *curr = (double *)calloc(n * n * n, sizeof(double));
  double *next = (double *)calloc(n * n * n, sizeof(double));

  // Ensure we haven't run out of memory
  if (curr == NULL || next == NULL) {
    fprintf(stderr,
            "Error: ran out of memory when trying to allocate %i sized cube\n",
            n);
    exit(EXIT_FAILURE);
  }

  memcpy(curr, source, n * n * n * sizeof(double));

  for (int iter = 0; iter < iterations; iter++) {
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
          int ip = i == n - 1 ? i - 1 : i + 1;
          int in = i == 0 ? 1 : i - 1;
          int jp = j == n - 1 ? j - 1 : j + 1;
          int jn = j == 0 ? 1 : j - 1;
          int kp = k == n - 1 ? k - 1 : k + 1;
          int kn = k == 0 ? 1 : k - 1;

          double source_term = delta * delta * source[INDEX(n, i, j, k)];
          next[INDEX(n, i, j, k)] =
              1.0 / 6.0 *
              (curr[INDEX(n, ip, j, k)] + curr[INDEX(n, in, j, k)] +
               curr[INDEX(n, i, jp, k)] + curr[INDEX(n, i, jn, k)] +
               curr[INDEX(n, i, j, kp)] + curr[INDEX(n, i, j, kn)] -
               source_term);
        }
      }
    }
    double *temp = curr;
    curr = next;
    next = temp;
  }

  // Free one of the buffers and return the correct answer in the other.
  // The caller is now responsible for free'ing the returned pointer.
  free(next);

  if (debug) {
    printf("Finished solving.\n");
  }

  return curr;
}

int main(int argc, char **argv) {

  // default settings for solver
  int iterations = 10;
  int n = 5;
  int threads = 1;
  float delta = 1;

  // parse the command line arguments
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
      printf(
          "usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
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
  double *source = (double *)calloc(n * n * n, sizeof(double));
  if (source == NULL) {
    fprintf(stderr, "Error: failed to allocated source term (n=%i)\n", n);
    return EXIT_FAILURE;
  }

  source[(n * n * n) / 2] = 1;

  // Calculate the resulting field with Neumann conditions
  double *result = poisson_neumann(n, source, iterations, threads, delta);

  // Print out the middle slice of the cube for validation
  for (int x = 0; x < n; ++x) {
    for (int y = 0; y < n; ++y) {
      printf("%0.5f ", result[((n / 2) * n + y) * n + x]);
    }
    printf("\n");
  }

  free(source);
  free(result);

  return EXIT_SUCCESS;
}
