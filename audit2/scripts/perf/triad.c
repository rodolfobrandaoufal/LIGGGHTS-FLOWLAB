/* STREAM-like triad a = b + s*c, N doubles per array (3 arrays), best of 10; OpenMP threads via OMP_NUM_THREADS */
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
int main(int argc, char **argv) {
  long n = argc > 1 ? atol(argv[1]) : 40000000L;
  double *a = malloc(n*8), *b = malloc(n*8), *c = malloc(n*8);
  #pragma omp parallel for
  for (long i = 0; i < n; i++) { a[i] = 0; b[i] = 1; c[i] = 2; }
  double best = 1e30;
  for (int r = 0; r < 10; r++) {
    double t = omp_get_wtime();
    #pragma omp parallel for
    for (long i = 0; i < n; i++) a[i] = b[i] + 3.0*c[i];
    t = omp_get_wtime() - t; if (t < best) best = t;
  }
  printf("threads %d triad %.1f GB/s (a[1]=%g)\n", omp_get_max_threads(), 24.0*n/best/1e9, a[1]);
  return 0;
}
