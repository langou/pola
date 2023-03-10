#include <math.h>

// qr_mgs_rl
// this is the right-looking (RL) variant of Modified Gram-Schmidt (MGS)
// this code is ``same`` as polybench, just that we require only one array for A

void qr_mgs_rl (int M, int N, double A[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for (i = 0; i < N; i++) {
// read A(1:M,i)
   R[i][i] = 0.0e+00;                     // RL2  LL5
   for (k = 0; k < M; k++)
      R[i][i] += A[k][i] * A[k][i];       // RL3  LL6
   R[i][i] = sqrt(R[i][i]);               // RL4  LL7
   for (k = 0; k < M; k++)
      A[k][i] /= R[i][i];                 // RL5  LL8
// write R[i][i]
   for (j = i + 1; j < N; j++) {
//    read A(1:M,j)
      R[i][j] = 0.0e+00;                  // RL6  LL2
      for (k = 0; k < M; k++)
         R[i][j] += A[k][i] * A[k][j];    // RL7  LL3
      for (k = 0; k < M; k++)
         A[k][j] -= A[k][i] * R[i][j];    // RL8  LL4
//    write R[i][j]
//    write A(1:M,j)
   }
// write A(1:M,i)
}
#pragma endscop
}
