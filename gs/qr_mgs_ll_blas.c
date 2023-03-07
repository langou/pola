#include <math.h>
#include <cblas.h>

// qr_mgs_ll
// this is the left-looking variant of MGS

void qr_mgs_ll_blas (int m, int n, double *A, int lda, double *R, int ldr )
{
int i, j;
for (j = 0; j < n; j++) {
   for (i = 0; i < j; i++) {
      R[i+j*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[j*lda]), 1 );
      cblas_daxpy( m, -R[i+j*ldr], &(A[i*lda]), 1, &(A[j*lda]), 1 );
   }                                    
   R[j+j*ldr] = cblas_ddot( m, &(A[j*lda]), 1, &(A[j*lda]), 1 );
   R[j+j*ldr] = sqrt(R[j+j*ldr]);
   cblas_dscal( m, 1.0/R[j+j*ldr], &(A[j*lda]), 1 );
}
}
