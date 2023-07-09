#include <math.h>

void qr_householder_v2q_ll____________ ( int M, int N, double *A, int lda, double *tau )
{
  int i, j, k;
  double tmp;

#pragma scop

for(k = N-1; k > -1; k--){

   A[k+lda*k] = 1.0e+00 - tau[k];

   for(i = k+1; i < M; i++){

      A[i+lda*k] = - A[i+lda*k] * tau[k];

   }

   for(j = k-1; j > -1; j--){

      tmp = 0.e+00;

      for(i = j+1; i < M; i++){

         tmp += A[i+lda*j] * A[i+lda*k];

      }

      tmp *= tau[j];

      A[j+lda*k] = -tmp;

      for(i = j+1; i < M; i++){

         A[i+lda*k] -= A[i+lda*j] * tmp;

      }

   }

}

#pragma endscop
}
