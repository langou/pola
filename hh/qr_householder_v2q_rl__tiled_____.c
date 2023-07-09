#include <math.h>

void qr_householder_v2q_rl__tiled_____ ( int M, int N, int B, double *A, int lda, double *tau )
{
  int i, j, k, k0;
  double tmp;

#pragma scop

for(k0 = N-1; k0 > -1; k0-=B){

   for(j = k0+1; j < N; j++){

      for (k = k0; ((k > k0-B)&&(k > -1)); k--) {

         tmp = 0.e+00;

         for(i = k+1; i < M; i++){

            tmp += A[i+lda*k] * A[i+lda*j];

         }

         tmp *= tau[k];


         A[k+lda*j] = -tmp;


         for(i = k+1; i < M; i++){

            A[i+lda*j] -= A[i+lda*k] * tmp;

         }
      }
   
   }

   for (k = k0; ((k > k0-B)&&(k > -1)); k--) {

      for(j = k+1; j < k0+1; j++){

         tau[j] = 0.e+00;

         for(i = k+1; i < M; i++){

            tau[j] += A[i+lda*k] * A[i+lda*j];

         }

      }

      for(j = k+1; j < k0+1; j++){ 

         tau[j] *= tau[k];

      }

      for(j = k+1; j < k0+1; j++){

         A[k+lda*j] = -tau[j];

      }

      for(j = k+1; j < k0+1; j++){

         for(i = k+1; i < M; i++){

            A[i+lda*j] -= A[i+lda*k] * tau[j];

         }
      }
   
      A[k+lda*k] = 1.0e+00 - tau[k];

      for(i = k+1; i < M; i++){

         A[i+lda*k] = - A[i+lda*k] * tau[k];

      }

   }

}
#pragma endscop
}
