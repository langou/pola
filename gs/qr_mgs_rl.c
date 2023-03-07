#include <math.h>

// qr_mgs_rl
// this is the right-looking (RL) variant of Modified Gram-Schmidt (MGS)
// this code is ``same`` as polybench, just that we require only one array for A

void qr_mgs_rl (int m, int n, double *A, int lda, double *R, int ldr )
{
int i, j, k;
#pragma scop
for (i = 0; i < n; i++) {
// read A(1:M,i)
   R[i+i*ldr] = 0.0e+00;                     // RL2  LL5
   for (k = 0; k < m; k++)
      R[i+i*ldr] += A[k+i*lda] * A[k+i*lda];       // RL3  LL6
   R[i+i*ldr] = sqrt(R[i+i*ldr]);               // RL4  LL7
   for (k = 0; k < m; k++)
      A[k+i*lda] /= R[i+i*ldr];                 // RL5  LL8
// write R[i][i]
   for (j = i + 1; j < n; j++) {
//    read A(1:M,j)
      R[i+j*ldr] = 0.0e+00;                  // RL6  LL2
      for (k = 0; k < m; k++)
         R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];    // RL7  LL3
      for (k = 0; k < m; k++)
         A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];    // RL8  LL4
//    write R[i][j]
//    write A(1:M,j)
   }
// write A(1:M,i)
}
#pragma endscop
}
