#include <math.h>
#include <cblas.h>

void qr_householder_a2v_ll__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau )
{
  int j, k;
  double norma2, norma, tmp;
  int k0;

#pragma scop
for (k0 = 0; k0 < N; k0+=B) {
// read A(1:M,k0:k0+B)
   int l; if ( k0+B < N ) l = B; else l = N-k0;
   for(j = 0; j < k0; j++){
//    read A(j:M,j)
      for (k = k0; k < k0+l; k++) tau[k] = A[j+lda*k];

      cblas_dgemv( CblasColMajor, CblasTrans, (M-j-1), l,
                  1.0, &(A[(j+1)+k0*lda]), lda, &(A[(j+1)+j*lda]), 1,
                  1.0, &(tau[k0]), 1 );

      for (k = k0; k < k0+l; k++) tau[k] = tau[j] * tau[k];
      for (k = k0; k < k0+l; k++) A[j+lda*k] = A[j+lda*k] - tau[k];

      cblas_dger( CblasColMajor, (M-j-1), l,
                -1.0, &(A[(j+1)+j*lda]), 1, &(tau[k0]), 1, &(A[(j+1)+k0*lda]), lda );

//    discard A(j:M,j)
   }

   for (k = k0; ((k < k0+B)&&(k < N)); k++) {
      for(j = k0; j < k; j++){
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
   
      cblas_dscal( (M-j-1), 1.0/A[k+lda*k], &(A[(k+1)+k*lda]), 1 );
      A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;
   }
// write A(1:M,k0:k0+B)
}
#pragma endscop
}
