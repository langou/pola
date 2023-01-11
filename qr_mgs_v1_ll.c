#include <math.h>

// qr_mgs_v1_ll
// this is the left-looking variant of MGS

void qr_mgs_v1_ll (int M, int N, double A[M][N], double R[N][N] )
{
int i, j, k;
#pragma scop
for (j = 0; j < N; j++) {
   for (i = 0; i < j; i++) {
      R[i][j] = 0.0e+00;                // LL2  RL6
      for (k = 0; k < M; k++)
         R[i][j] += A[k][i] * A[k][j];  // LL3  RL7
      for (k = 0; k < M; k++)           
         A[k][j] -= A[k][i] * R[i][j];  // LL4  RL8
   }                                    
   R[j][j] = 0.0e+00;                   // LL5  RL2   
   for (k = 0; k < M; k++)              
      R[j][j] += A[k][j] * A[k][j];     // LL6  RL3    
   R[j][j] = sqrt(R[j][j]);             // LL7  RL4
   for (k = 0; k < M; k++)
      A[k][j] /= R[j][j];               // LL8  RL5
}
#pragma endscop
}
