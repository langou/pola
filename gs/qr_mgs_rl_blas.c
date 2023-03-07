#include <math.h>
#include <cblas.h>

// qr_mgs_rl
// this is the right-looking (RL) variant of Modified Gram-Schmidt (MGS)
// this code is ``same`` as polybench, just that we require only one array for A

void qr_mgs_rl_blas (int m, int n, double *A, int lda, double *R, int ldr )
{
int i, j;
#pragma scop
for (i = 0; i < n; i++) {
   R[i+i*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[i*lda]), 1 );
   R[i+i*ldr] = sqrt(R[i+i*ldr]);
   cblas_dscal( m, 1.0/R[i+i*ldr], &(A[i*lda]), 1 );
   for (j = i + 1; j < n; j++) {
      R[i+j*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[j*lda]), 1 );
      cblas_daxpy( m, -R[i+j*ldr], &(A[i*lda]), 1, &(A[j*lda]), 1 );
   }
}
#pragma endscop
}
