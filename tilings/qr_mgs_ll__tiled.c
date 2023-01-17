#include <math.h>

// qr_mgs_ll__tiled
// this is the left-looking variant of MGS

void qr_mgs_ll__tiled (int M, int N, int B, double A[M][N], double R[N][N] )
{
int i, j, k;
int j0;
#pragma scop
for (j0 = 0; j0 < N; j0+=B) {
// read A(1:M,i0:i0+B)
   for (i = 0; i < j0; i++) {
//    read A(1:M,i)
      for (j = j0; ((j < j0+B)&&(j < N)); j++) {
         R[i][j] = 0.0e+00;
         for (k = 0; k < M; k++)
            R[i][j] += A[k][i] * A[k][j];
         for (k = 0; k < M; k++)
            A[k][j] -= A[k][i] * R[i][j];
      }
//    discard A(1:M,i)
   }
   for (j = j0; ((j < j0+B)&&(j < N)); j++) {
      for (i = j0; i < j; i++) {
         R[i][j] = 0.0e+00;
         for (k = 0; k < M; k++)
            R[i][j] += A[k][i] * A[k][j];
         for (k = 0; k < M; k++)
            A[k][j] -= A[k][i] * R[i][j];
      }
      R[j][j] = 0.0e+00;
      for (k = 0; k < M; k++)
         R[j][j] += A[k][j] * A[k][j];
      R[j][j] = sqrt(R[j][j]);
      for (k = 0; k < M; k++)
         A[k][j] /= R[j][j];
   }
// write A(1:M,i0:i0+B)
}
#pragma endscop
}
