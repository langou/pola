#include <math.h>

void qr_householder_v0 (int M, int N, int P, double **A )
//void qr_householder_v0 (int M, int N, int P, double A[M][N] )
{
  int i, j, k;
  double norma2, norma, tau, tmp;

#pragma scop

for(k = 0; k < P; k++){
   norma2 = 0.e+00;
   for(i = k+1; i < M; i++){
      norma2 += A[i][k] * A[i][k];
   }
   norma = sqrt( A[k][k] * A[k][k] + norma2 );
   
   A[k][k] = ( A[k][k] > 0 ) ? ( A[k][k] + norma ) : ( A[k][k] - norma ) ;
   
   tau = 2.0 / ( 1.0 + norma2 / ( A[k][k] * A[k][k] ) ) ;
   
   for(i = k+1; i < M; i++){
      A[i][k] /= A[k][k];
   }
   A[k][k]= ( A[k][k] > 0 ) ? ( - norma ) : ( norma ) ;
   
   for(j = k+1; j < N; j++){
      tmp = A[k][j];
      for(i = k+1; i < M; i++){
            tmp += A[i][k] * A[i][j];
      }
      tmp = tau * tmp;
      A[k][j] = A[k][j] - tmp;
      for(i = k+1; i < M; i++){
         A[i][j] = A[i][j] - A[i][k] * tmp;
      }
   }
}
#pragma endscop

}

// => M*P + 1/2*N^2 + M^2*S^(-1/2)*N
// Note: parameter starting with "_" makes Ginac bug
