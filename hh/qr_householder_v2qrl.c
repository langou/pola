#include <math.h>

void qr_householder_v2qrl ( int M, int N, double *A, int lda, double *tau )
{
  int i, j, k;

#pragma scop

for(k = N-1; k > -1; k--){

   for(j = k+1; j < N; j++){

      tau[j] = 0.e+00;

         for(i = k+1; i < M; i++){

            tau[j] += A[i+lda*k] * A[i+lda*j];

         }

   }

   for(j = k+1; j < N; j++){ 

      tau[j] *= tau[k];

   }

   for(j = k+1; j < N; j++){

      A[k+lda*j] = -tau[j];

   }
   for(j = k+1; j < N; j++){

      for(i = k+1; i < M; i++){

         A[i+lda*j] -= A[i+lda*k] * tau[j];

      }
   }
   
   A[k+lda*k] = 1.0e+00 - tau[k];

   for(i = k+1; i < M; i++){
      A[i+lda*k] = - A[i+lda*k] * tau[k];
   }

}
#pragma endscop
}
