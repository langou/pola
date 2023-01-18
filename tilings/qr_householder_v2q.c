#include <math.h>

void qr_householder_v2q ( int M, int N, double A[M][N], double tau[N] )
{
  int i, j, k;

#pragma scop

for(k = N-1; k > -1; k--){

   for(j = k+1; j < N; j++){

      tau[j] = 0.e+00;

         for(i = k+1; i < M; i++){

            tau[j] += A[i][k] * A[i][j];

         }

   }

   for(j = k+1; j < N; j++){ 

      tau[j] *= tau[k];

   }

   A[k][k] = 1.0e+00 - tau[k];

   for(j = k+1; j < N; j++){

      A[k][j] = -tau[j];

   }
   for(j = k+1; j < N; j++){

      for(i = k+1; i < M; i++){

         A[i][j] -= A[i][k] * tau[j];

      }
   }
   
   for(i = k+1; i < M; i++){
      A[i][k] = - A[i][k] * tau[k];
   }


}
#pragma endscop
}
