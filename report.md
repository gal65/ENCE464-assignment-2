## Inspected cache misses with no optimization:
- Inner loop register load (into xmm) misses cache approximately 6% of the time
- This is expected since cache lines are 64 bytes, and floats are 4 bytes
- 16 floats can fit in cache line, so 1/16th of the time (6.25%) a new line
  will have to be loaded.

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

Both produce identical ASM output, so the more readable form was kept. Compiler
knew that the indices were incrementing each loop and changed the indexing into
load/increment.

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
Originally, there were two versions of the calculations: one with boundary
check and one without. It was found that the boundary check branches were
getting optimized out anyway, so a single function was used. This gets inlined
anyway so its nbd.