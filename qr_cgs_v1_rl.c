#include <math.h>

// qr_cgs_v1_rl

void qr_cgs_v1_rl (int M, int N, double A[M][N], double Q[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for (i = 0; i < N; i++) {
   R[i][i] = 0.0e+00;
   for (k = 0; k < M; k++)
      R[i][i] += Q[k][i] * Q[k][i];
   R[i][i] = sqrt(R[i][i]);
   for (k = 0; k < M; k++)
      Q[k][i] /= R[i][i];
   for (j = i + 1; j < N; j++) {
      R[i][j] = 0.0e+00;
      for (k = 0; k < M; k++)
//       R[i][j] += Q[k][i] * Q[k][j]; // MGS
         R[i][j] += Q[k][i] * A[k][j]; // CGS
   }
   for (j = i + 1; j < N; j++) {
      for (k = 0; k < M; k++)
         Q[k][j] -= Q[k][i] * R[i][j];
   }
}
#pragma endscop
}
