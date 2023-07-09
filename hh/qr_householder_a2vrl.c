#include <math.h>

void qr_householder_a2vrl ( int M, int N, double *A, int lda, double *tau )
{
  int i, j, k;
  double norma2, norma, tmp;

#pragma scop
for(k = 0; k < N; k++){

   norma2 = 0.e+00;
   for(i = k+1; i < M; i++){
      norma2 += A[i+lda*k] * A[i+lda*k];
   }
   norma = sqrt( A[k+lda*k] * A[k+lda*k] + norma2 );
   
   A[k+lda*k] = ( A[k+lda*k] > 0 ) ? ( A[k+lda*k] + norma ) : ( A[k+lda*k] - norma ) ;
   
   tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k+lda*k] * A[k+lda*k] ) ) ;
   
   for(i = k+1; i < M; i++){
      A[i+lda*k] /= A[k+lda*k];
   }
   A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;

   for(j = k+1; j < N; j++){
      tmp = A[k+lda*j];
      for(i = k+1; i < M; i++){
            tmp += A[i+lda*k] * A[i+lda*j];
      }
      tmp = tau[k] * tmp;
      A[k+lda*j] = A[k+lda*j] - tmp;
      for(i = k+1; i < M; i++){
         A[i+lda*j] = A[i+lda*j] - A[i+lda*k] * tmp;
      }
   }
   
}
#pragma endscop
}
