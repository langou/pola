#include <math.h>
#include <cblas.h>

void qr_householder_a2vll_blas ( int M, int N, double *A, int lda, double *tau )
{
  int j, k;
  double norma2, norma, tmp;

#pragma scop
for(k = 0; k < N; k++){

   for(j = 0; j < k; j++){
      tmp = A[j+lda*k];
      tmp += cblas_ddot( (M-j-1), &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );
      tmp = tau[j] * tmp;
      A[j+lda*k] = A[j+lda*k] - tmp;
      cblas_daxpy( (M-j-1), -tmp, &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );
   }

   norma2 = cblas_ddot( (M-j-1), &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );
   norma = sqrt( A[k+lda*k] * A[k+lda*k] + norma2 );
   
   A[k+lda*k] = ( A[k+lda*k] > 0 ) ? ( A[k+lda*k] + norma ) : ( A[k+lda*k] - norma ) ;
   
   tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k+lda*k] * A[k+lda*k] ) ) ;
   
   cblas_dscal( (M-k-1), 1.0/A[k+lda*k], &(A[(k+1)+k*lda]), 1 );
   A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;
   

}
#pragma endscop
}
