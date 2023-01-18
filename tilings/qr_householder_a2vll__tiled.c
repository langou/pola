#include <math.h>

void qr_householder_a2vll__tiled ( int M, int N, int B, double A[M][N], double tau[N] )
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
         tmp = A[j][k];
         for(i = j+1; i < M; i++)
            tmp += A[i][j] * A[i][k];
         tmp = tau[j] * tmp;
         A[j][k] = A[j][k] - tmp;
         for(i = j+1; i < M; i++)
            A[i][k] = A[i][k] - A[i][j] * tmp;
      }
//    discard A(j:M,j)
   }

   for (k = k0; ((k < k0+B)&&(k < N)); k++) {
      for(j = k0; j < k; j++){
         tmp = A[j][k];
         for(i = j+1; i < M; i++)
            tmp += A[i][j] * A[i][k];
         tmp = tau[j] * tmp;
         A[j][k] = A[j][k] - tmp;
         for(i = j+1; i < M; i++)
            A[i][k] = A[i][k] - A[i][j] * tmp;
      }

      norma2 = 0.e+00;
      for(i = k+1; i < M; i++)
         norma2 += A[i][k] * A[i][k];
      norma = sqrt( A[k][k] * A[k][k] + norma2 );
   
      A[k][k] = ( A[k][k] > 0 ) ? ( A[k][k] + norma ) : ( A[k][k] - norma ) ;
   
      tau[k] = 2.0 / ( 1.0 + norma2 / ( A[k][k] * A[k][k] ) ) ;
   
      for(i = k+1; i < M; i++)
         A[i][k] /= A[k][k];
      A[k][k]= ( A[k][k] > 0 ) ? ( - norma ) : ( norma ) ;
   }
// write A(1:M,k0:k0+B)
}
#pragma endscop
}
