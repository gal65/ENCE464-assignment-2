#include <stdlib.h>
#include <stdio.h>


void main()
{
    int n = 301;
    double* in = (double*)calloc(n * n * n, sizeof(double));
    double* out = (double*)calloc(n * n * n, sizeof(double));
    for (int iter = 0; iter < 10; iter++) {
                for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
        for (int k = 1; k < n - 1; k++) {
                    out[k*n*n + j*n + i] = 1.0 / 6.0 * (
                        in[k*n*n + j*n + (i+1)] +
                        in[k*n*n + j*n + (i-1)] +
                        in[k*n*n + (j+1)*n + i] +
                        in[k*n*n + (j-1)*n + i] +
                        in[(k+1)*n*n + j*n + i] +
                        in[(k-1)*n*n + j*n + i] );
                }
            }
        }
    }
    printf("%f\n", out[n*n*n / 2]);
}