## Inspected cache misses with no optimization:
- Inner loop register load (into xmm) misses cache approximately 6% of the time
- This is expected since cache lines are 64 bytes, and floats are 4 bytes
- 16 floats can fit in cache line, so 1/16th of the time (6.25%) a new line
  will have to be loaded.
- This assumes only one cache line is used per float being read in the inner loop. This is the case without optimizations turned on, probably because function arguments (worker_args_t*) which are loaded from stack into cache.

## Pointer arith optimization:
- Tried with pointer arithmetic, where each pointer increments (++) after the
  inner loop body, and is incremented by 2 to skip over the boundary layer
```c
for (int i = 1; i < n - 1; i++) {
    double source_term = args->delta_squared * (*source_ptr++);
    middle++;
    // REMEMBER: "middle" is already pointing to the right neighbor when the
    // relaxation is calulated.
    *next_ptr++ = 1.0 / 6.0 *
        ((*top++) + (*bottom++) + (*front++) + (*back++) +
        (*(middle - 2)) + (*middle) - source_term);
}
```

- Also tried with "naive" indexing:
```c
for (int i = 1; i < n - 1; i++) {
    float source_term = args->delta_squared * args->source[IDX(n, i, j, k)];
    next[IDX(n, i, j, k)] =
        1.0 / 6.0 *
        (curr[IDX(n, i + 1, j, k)] + curr[IDX(n, i - 1, j, k)] +
         curr[IDX(n, i, j + 1, k)] + curr[IDX(n, i, j - 1, k)] +
         curr[IDX(n, i, j, k + 1)] + curr[IDX(n, i, j, k - 1)] - source_term);
}
```

On average:
### Ptr arith
- cycles: 2.19T
## indexing:
- cycles: 2.25T

## Aligning to cache lines:
- malloc/calloc is not algined to cache line (64 bytes) by default. This could
  cause up to an extra 2 cache lines needed per row of the array to be loaded,
  if there was unalignment at both ends. To remedy: use:
```c
float *curr, *next;
// Allocate curr and next to be cache aligned
posix_memalign((void**)&curr, CACHE_LINE_SIZE, n * n * n * sizeof(float));
posix_memalign((void**)&next, CACHE_LINE_SIZE, n * n * n * sizeof(float));
```

## Calculating boundary conditions
Inside the main loop, we can check for boundaries to re-use code from the boundary worker (DRY). Or we can special case it. special casing is only a tiny bit better.

### Check boundaries even though unnecessary
- cycle est = 2.176T
### Dont check boundaries
- cycle est = 2.177T (better ish?)

## Non-overlapping boundary condition
### Impl 1
- do multiple for loops over surface of cube, such that no loops overlap and
  calculate the same value twice
- instruction fetches: 412 675 326
- L1 cache misses: 8 081 252
- cycle est = 524.9M

### Impl 2
- dont care about calculating the same boundary condition twice, just loop over
  full faces of the cube, even if they overlap at the corners
- instruction fetches 694 405 338
- L1 cache misses: 6 324 050
- cycle est = 434.8M

Verdict: "naive" approach is once again better

## Pre-fetching cache lines
use `__builtin_prefetch` in inner loop of worker thread to load the next 64
bytes into cache before they are needed.

Decreased performance: it made the vectorization code (SIMD) more complicated,
and therefore slower. The CPU is better at prefetching than we are

## Marking do_cell as inline
Only takes effect in boundary worker, as it is always inlines in `worker`.

Effects in boundary worker:
### with inline
- instr fetch = 163M
- L1 miss = 6 305 000
- L3 miss = 2 005 873
- cycle estimation = 426M

### without inline
- instr fetch = 122M
- L1 miss = 6 323 908
- L3 miss = 2 023 103
- cycle estimation = 659M

## also cover in report
- compiler flags
- AB test makefile
- preprocessor args
- loop tiling and why it doesn't work