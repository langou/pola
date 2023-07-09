#include <math.h>
#include <cblas.h>

void qr_householder_a2v_rl________blas ( int M, int N, double *A, int lda, double *tau )
{
  int j, k;
  double norma2, norma;

#pragma scop
for(k = 0; k < N; k++){

   norma2 = cblas_ddot( (M-k-1), &(A[(k+1)+k*lda]), 1, &(A[(k+1)+k*lda]), 1 );
   norma = sqrt( A[k+lda*k] * A[k+lda*k] + norma2 );
   
   A[k+lda*k] = ( A[k+lda*k] > 0 ) ? ( A[k+lda*k] + norma ) : ( A[k+lda*k] - norma ) ;
   
   tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k+lda*k] * A[k+lda*k] ) ) ;
   
   cblas_dscal( (M-k-1), 1.0/A[k+lda*k], &(A[(k+1)+k*lda]), 1 );
   A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;

   for(j = k+1; j < N; j++) tau[j] = A[k+lda*j];
   
   cblas_dgemv( CblasColMajor, CblasTrans, (M-k-1), (N-k-1),
                1.0, &(A[(k+1)+(k+1)*lda]), lda, &(A[(k+1)+k*lda]), 1,
                1.0, &(tau[k+1]), 1 );

   for(j = k+1; j < N; j++) tau[j] = tau[k] * tau[j];
   for(j = k+1; j < N; j++) A[k+lda*j] = A[k+lda*j] - tau[j];

  cblas_dger( CblasColMajor, (M-k-1), (N-k-1),
     -1.0, &(A[(k+1)+k*lda]), 1, &(tau[k+1]), 1, &(A[(k+1)+(k+1)*lda]), lda );

   
}
#pragma endscop
}
