#include <math.h>

// qr_cgs2_v1

void qr_cgs2_v1 (int M, int N, double Q[M][N], double R[N][N], double tmp[N] )
{
  int i, j, k;
  double nrm;

#pragma scop
for (k = 0; k < N; k++) {
   nrm = 0.0e+00;
   for (i = 0; i < M; i++)
      nrm += Q[i][k] * Q[i][k];
   R[k][k] = sqrt(nrm);
   for (i = 0; i < M; i++)
      Q[i][k] /= R[k][k];
   for (j = k + 1; j < N; j++) {
      R[k][j] = 0.0e+00;
      for (i = 0; i < M; i++)
         R[k][j] += Q[i][k] * Q[i][j];
   }
   for (j = k + 1; j < N; j++) {
      for (i = 0; i < M; i++)
         Q[i][j] -= Q[i][k] * R[k][j];
   }
   for (j = k + 1; j < N; j++) {
      tmp[j] = 0.0e+00;
      for (i = 0; i < M; i++)
         tmp[j] += Q[i][k] * Q[i][j];
      R[k][j] += tmp[j];
   }
   for (j = k + 1; j < N; j++) {
      for (i = 0; i < M; i++)
         Q[i][j] -= Q[i][k] * tmp[j];
   }
}
#pragma endscop

}

