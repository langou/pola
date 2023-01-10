#include <math.h>

// qr_mgs_v1

void qr_mgs_v1 (int M, int N, double A[M][N], double R[N][N] )
{
  int i, j, k;
  double nrm;

#pragma scop
for (k = 0; k < N; k++) {
   nrm = 0.0e+00;
   for (i = 0; i < M; i++)
      nrm += A[i][k] * A[i][k];
   R[k][k] = sqrt(nrm);
   for (i = 0; i < M; i++)
      A[i][k] /= R[k][k];
   for (j = k + 1; j < N; j++) {
      R[k][j] = 0.0e+00;
      for (i = 0; i < M; i++)
         R[k][j] += A[i][k] * A[i][j];
      for (i = 0; i < M; i++)
         A[i][j] -= A[i][k] * R[k][j];
   }
}
#pragma endscop

}

