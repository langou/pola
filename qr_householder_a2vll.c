#include <math.h>

void qr_householder_a2vll ( int M, int N, double A[M][N], double tau[N] )
{
  int i, j, k;
  double norma2, norma;

#pragma scop

for(k = 0; k < N; k++){

   for(j = 0; j < k-1; j++){
      tau[j] = A[k][j];
      for(i = k+1; i < M; i++){
            tau[j] += A[i][k] * A[i][j];
      }
      tau[j] = tau[k] * tau[j];
      A[k][j] = A[k][j] - tau[j];
      for(i = k+1; i < M; i++){
         A[i][j] = A[i][j] - A[i][k] * tau[j];
      }
   }

   norma2 = 0.e+00;
   for(i = k+1; i < M; i++){
      norma2 += A[i][k] * A[i][k];
   }
   norma = sqrt( A[k][k] * A[k][k] + norma2 );
   
   A[k][k] = ( A[k][k] > 0 ) ? ( A[k][k] + norma ) : ( A[k][k] - norma ) ;
   
   tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k][k] * A[k][k] ) ) ;
   
   for(i = k+1; i < M; i++){
      A[i][k] /= A[k][k];
   }
   A[k][k]= ( A[k][k] > 0 ) ? ( - norma ) : ( norma ) ;
   

}
#pragma endscop
}
