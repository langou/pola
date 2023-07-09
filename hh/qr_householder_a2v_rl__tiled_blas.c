#include <math.h>
#include <cblas.h>


void qr_householder_a2v_rl__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau )
{
  int j, k, k0;
  double norma2, norma, tmp;

#pragma scop
for(k0 = 0; k0 < N; k0+=B){

   int l; if ( k0+B < N ) l = B; else l = N-k0;

   for (k = k0; ((k < k0+B)&&(k < N)); k++) {

      norma2 = cblas_ddot( (M-k-1), &(A[(k+1)+k*lda]), 1, &(A[(k+1)+k*lda]), 1 );
      norma = sqrt( A[k+lda*k] * A[k+lda*k] + norma2 );
   
      A[k+lda*k] = ( A[k+lda*k] > 0 ) ? ( A[k+lda*k] + norma ) : ( A[k+lda*k] - norma ) ;
   
      tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k+lda*k] * A[k+lda*k] ) ) ;
   
      cblas_dscal( (M-k-1), 1.0/A[k+lda*k], &(A[(k+1)+k*lda]), 1 );
      A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;

      for(j = k+1; j < k0+l; j++) tau[j] = A[k+lda*j];
   
      cblas_dgemv( CblasColMajor, CblasTrans, (M-k-1), (k0+l-k-1),
                   1.0, &(A[(k+1)+(k+1)*lda]), lda, &(A[(k+1)+k*lda]), 1,
                   1.0, &(tau[k+1]), 1 );

      for(j = k+1; j < k0+l; j++) tau[j] = tau[k] * tau[j];
      for(j = k+1; j < k0+l; j++) A[k+lda*j] = A[k+lda*j] - tau[j];

      cblas_dger( CblasColMajor, (M-k-1), (k0+l-k-1),
                  -1.0, &(A[(k+1)+k*lda]), 1, &(tau[k+1]), 1, &(A[(k+1)+(k+1)*lda]), lda );
   
   }

   for(j = k0+B; j < N; j++){
      for (k = k0; ((k < k0+B)&&(k < N)); k++) {
         tmp = A[k+lda*j];
         tmp += cblas_ddot( (M-k-1), &(A[(k+1)+k*lda]), 1, &(A[(k+1)+j*lda]), 1 );
         tmp = tau[k] * tmp;
         A[k+lda*j] = A[k+lda*j] - tmp;
         cblas_daxpy( (M-k-1), -tmp, &(A[(k+1)+k*lda]), 1, &(A[(k+1)+j*lda]), 1 );
      }
   }

}
#pragma endscop
}
