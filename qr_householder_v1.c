#include <math.h>

// this version has a workarray tmp, this might help with parallelism or data movement
// this version is closer from the LAPACK GEQR2 version

void qr_householder_v1 (int M, int N, int P, double **A, double *tmp )
//void qr_householder_v1 (int M, int N, int P, double A[M][N], tmp[N] )
{
  int i, j, k;
  double norma2, norma, tau;

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
      tmp[j] = A[k][j];
      for(i = k+1; i < M; i++){
            tmp[j] += A[i][k] * A[i][j];
      }
      tmp[j] = tau * tmp[j];
      A[k][j] = A[k][j] - tmp[j];
      for(i = k+1; i < M; i++){
         A[i][j] = A[i][j] - A[i][k] * tmp[j];
      }
   }
}
#pragma endscop

}

// => M*P + 1/2*N^2 + M^2*S^(-1/2)*N
