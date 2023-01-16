#include <math.h>

// qr_cgs_rl
//
// this is the right-looking (RL) version of Classical Gram-Schmidt (CGS)
//
// array A is input only
// arrays Q and R do not need to be initialized in input, they are ouput only
// the strictly lower part of R is not referenced

void qr_cgs_rl__tiled (int M, int N, int B, double A[M][N], double Q[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for(i = 0; i < M; i++) for(j = 0; j < N; j++) Q[i][j] = A[i][j];
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
         R[i][j] += Q[k][i] * A[k][j]; // CGS
   }
   for (j = i + 1; j < N; j++) {
      for (k = 0; k < M; k++)
         Q[k][j] -= Q[k][i] * R[i][j];
   }
}
#pragma endscop
}
