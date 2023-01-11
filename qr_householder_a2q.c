#include <math.h>

void qr_householder_a2q ( int M, int N, double A[M][N], double R[N][N], double work[N] )
{
  int i, j, k;
  double norma2, norma;

#pragma scop

for(k = 0; k < N; k++){
   norma2 = 0.e+00;
   for(i = k+1; i < M; i++){
      norma2 += A[i][k] * A[i][k];
   }
   norma = sqrt( A[k][k] * A[k][k] + norma2 );
   
   A[k][k] = ( A[k][k] > 0 ) ? ( A[k][k] + norma ) : ( A[k][k] - norma ) ;
   
   work[k] = 2.0 / ( 1.0 + norma2 / ( A[k][k] * A[k][k] ) ) ;
   
   for(i = k+1; i < M; i++){
      A[i][k] /= A[k][k];
   }
   A[k][k]= ( A[k][k] > 0 ) ? ( - norma ) : ( norma ) ;
   
   for(j = k+1; j < N; j++){
      work[j] = A[k][j];
      for(i = k+1; i < M; i++){
            work[j] += A[i][k] * A[i][j];
      }
      work[j] = work[k] * work[j];
      A[k][j] = A[k][j] - work[j];
      for(i = k+1; i < M; i++){
         A[i][j] = A[i][j] - A[i][k] * work[j];
      }
   }
}

for(i = 0; i < N; i++)
   for(j = i; j < N; j++)
      R[i][j] = A[i][j];

for(k = N-1; k > -1; k--){

   for(j = k+1; j < N; j++){

      work[j] = 0.e+00;

         for(i = k+1; i < M; i++){

            work[j] += A[i][k] * A[i][j];

         }

   }

   for(j = k+1; j < N; j++){ 

      work[j] *= work[k];

   }

   A[k][k] = 1.0e+00 - work[k];

   for(j = k+1; j < N; j++){

      A[k][j] = -work[j];

   }
   for(j = k+1; j < N; j++){

      for(i = k+1; i < M; i++){

         A[i][j] -= A[i][k] * work[j];

      }
   }
   
   for(i = k+1; i < M; i++){
      A[i][k] = - A[i][k] * work[k];
   }


}
#pragma endscop
}
