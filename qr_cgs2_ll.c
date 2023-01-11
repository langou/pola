#include <math.h>

// qr_cgs2_ll

void qr_cgs2_ll (int M, int N, double Q[M][N], double R[N][N], double tmp[N] )
{
  int i, j, k;

#pragma scop
for (j = 0; j < N; j++) {
   for (i = 0; i < j; i++) {
      R[i][j] = 0.0e+00;
      for (k = 0; k < M; k++)
         R[i][j] += Q[k][i] * Q[k][j];
   }
   for (i = 0; i < j; i++)
      for (k = 0; k < M; k++)
         Q[k][j] -= Q[k][i] * R[i][j];
   for (i = 0; i < j; i++) {
      tmp[i] = 0.0e+00;
      for (k = 0; k < M; k++)
         tmp[i] += Q[k][i] * Q[k][j];
      R[i][j] += tmp[i];
   }
   for (i = 0; i < j; i++)
      for (k = 0; k < M; k++)
         Q[k][j] -= Q[k][i] * tmp[i];
   R[j][j] = 0.0e+00;
   for (k = 0; k < M; k++)
      R[j][j] += Q[k][j] * Q[k][j];
   R[j][j] = sqrt(R[j][j]);
   for (k = 0; k < M; k++)
      Q[k][j] /= R[j][j];
}
#pragma endscop

}

