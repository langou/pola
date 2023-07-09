#include <math.h>
#include <cblas.h>

void qr_householder_v2q_ll________blas ( int M, int N, double *A, int lda, double *tau )
{
   int j, k;
   double tmp;

#pragma scop

for(k = N-1; k > -1; k--){

   A[k+lda*k] = 1.0e+00 - tau[k];

   cblas_dscal( (M-k-1), -tau[k], &(A[(k+1)+k*lda]), 1 );

   for(j = k-1; j > -1; j--){

      tmp = cblas_ddot( (M-j-1), &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );

      tmp *= tau[j];

      A[j+lda*k] = -tmp;

      cblas_daxpy( (M-j-1), -tmp, &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );

   }

}

#pragma endscop
}
