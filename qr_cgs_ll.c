#include <math.h>

// qr_cgs_ll 

void qr_cgs_ll (int M, int N, double A[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for (j = 0; j < N; j++) {
   for (i = 0; i < j; i++) {
      R[i][j] = 0.0e+00;
      for (k = 0; k < M; k++)
         R[i][j] += A[k][i] * A[k][j];
   }
   for (i = 0; i < j; i++)
      for (k = 0; k < M; k++)
         A[k][j] -= A[k][i] * R[i][j];
   R[j][j] = 0.0e+00;
   for (k = 0; k < M; k++)
      R[j][j] += A[k][j] * A[k][j];
   R[j][j] = sqrt(R[j][j]);
   for (k = 0; k < M; k++)
      A[k][j] /= R[j][j];
}
#pragma endscop
}
