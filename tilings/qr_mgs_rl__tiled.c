#include <math.h>
#include <stdio.h>

// qr_mgs_rl

void qr_mgs_rl__tiled (int M, int N, int B, double A[M][N], double R[N][N] )
{
int i, j, k;
int i0;

// the coefficients of R are computed, used only once and right after
// their computation, so we write them as soon as we are done with them
// we write the elements of R on the fly, and so 
// (1) the total number of writes due to R is N * (N+1) / 2
// (2) the requirement of cache for R is 1 for the whole computation, so
// we do not take R at all for the cache useage.

#pragma scop
for (i0 = 0; i0 < N; i0+=B) {
////////////////////////////////////////////////////////////
// read A(1:M,i0:i0+B)
   for (i = i0; ((i < i0+B)&&(i < N)); i++) {
      R[i][i] = 0.0e+00; 
      for (k = 0; k < M; k++)
         R[i][i] += A[k][i] * A[k][i];
      R[i][i] = sqrt(R[i][i]);
      for (k = 0; k < M; k++)
         A[k][i] /= R[i][i];
      for (j = i + 1; ((j < i0+B)&&(j < N)); j++) {
         R[i][j] = 0.0e+00;
         for (k = 0; k < M; k++)
            R[i][j] += A[k][i] * A[k][j];
         for (k = 0; k < M; k++)
            A[k][j] -= A[k][i] * R[i][j];
      }
   }
////////////////////////////////////////////////////////////
   for (j = i0+B; j < N; j++) {
//    load A(1:M,j)
      for (i = i0; ((i < i0+B)&&(i < N)); i++) {
         R[i][j] = 0.0e+00;
         for (k = 0; k < M; k++)
            R[i][j] += A[k][i] * A[k][j];
         for (k = 0; k < M; k++)
            A[k][j] -= A[k][i] * R[i][j];
      }
//    write A(1:M,j)
   }
// write A(1:M,i0:i0+B)
}
#pragma endscop
}
