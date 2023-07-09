#include <math.h>

void qr_householder_a2v_ll__tiled_____ ( int M, int N, int B, double *A, int lda, double *tau )
{
  int i, j, k;
  double norma2, norma, tmp;
  int k0;

#pragma scop
for (k0 = 0; k0 < N; k0+=B) {
// read A(1:M,k0:k0+B)
   for(j = 0; j < k0; j++){
//    read A(j:M,j)
      for (k = k0; ((k < k0+B)&&(k < N)); k++) {
         tmp = A[j+lda*k];
         for(i = j+1; i < M; i++)
            tmp += A[i+lda*j] * A[i+lda*k];
         tmp = tau[j] * tmp;
         A[j+lda*k] = A[j+lda*k] - tmp;
         for(i = j+1; i < M; i++)
            A[i+lda*k] = A[i+lda*k] - A[i+lda*j] * tmp;
      }
//    discard A(j:M,j)
   }

   for (k = k0; ((k < k0+B)&&(k < N)); k++) {
      for(j = k0; j < k; j++){
         tmp = A[j+lda*k];
         for(i = j+1; i < M; i++)
            tmp += A[i+lda*j] * A[i+lda*k];
         tmp = tau[j] * tmp;
         A[j+lda*k] = A[j+lda*k] - tmp;
         for(i = j+1; i < M; i++)
            A[i+lda*k] = A[i+lda*k] - A[i+lda*j] * tmp;
      }

      norma2 = 0.e+00;
      for(i = k+1; i < M; i++)
         norma2 += A[i+lda*k] * A[i+lda*k];
      norma = sqrt( A[k+lda*k] * A[k+lda*k] + norma2 );
   
      A[k+lda*k] = ( A[k+lda*k] > 0 ) ? ( A[k+lda*k] + norma ) : ( A[k+lda*k] - norma ) ;
   
      tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k+lda*k] * A[k+lda*k] ) ) ;
   
      for(i = k+1; i < M; i++)
         A[i+lda*k] /= A[k+lda*k];
      A[k+lda*k]= ( A[k+lda*k] > 0 ) ? ( - norma ) : ( norma ) ;
   }
// write A(1:M,k0:k0+B)
}
#pragma endscop
}
