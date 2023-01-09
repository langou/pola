#include <math.h>

// qr_mgs_v1

//void qr_mgs_v1 (int M, int N, double **Q, double **R )
void qr_mgs_v1 (int M, int N, double Q[M][N], double R[N][N] )
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
      Q[i][k] = Q[i][k] / R[k][k];
   for (j = k + 1; j < N; j++) {
      R[k][j] = 0.0e+00;
      for (i = 0; i < M; i++)
         R[k][j] += Q[i][k] * Q[i][j];
      for (i = 0; i < M; i++)
         Q[i][j] = Q[i][j] - Q[i][k] * R[k][j];
   }
}
#pragma endscop

}

