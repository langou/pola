#include <math.h>

// qr_mgs_v1_rl
// this is the right-looking variant of MGS
// this code is ``same`` as polybench

void qr_mgs_v1_rl (int M, int N, double A[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for (i = 0; i < N; i++) {
   R[i][i] = 0.0e+00;                     // RL2  LL5
   for (k = 0; k < M; k++)
      R[i][i] += A[k][i] * A[k][i];       // RL3  LL6
   R[i][i] = sqrt(R[i][i]);               // RL4  LL7
   for (k = 0; k < M; k++)
      A[k][i] /= R[i][i];                 // RL5  LL8
   for (j = i + 1; j < N; j++) {
      R[i][j] = 0.0e+00;                  // RL6  LL2
      for (k = 0; k < M; k++)
         R[i][j] += A[k][i] * A[k][j];    // RL7  LL3
      for (k = 0; k < M; k++)
         A[k][j] -= A[k][i] * R[i][j];    // RL8  LL4
   }
}
#pragma endscop
}
