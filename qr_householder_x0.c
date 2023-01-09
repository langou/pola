#include <math.h>

// this version has a workarray tmp, this might help with parallelism or data movement
// this version is closer from the LAPACK GEQR2 version

void qr_householder_x0 (int M, int N, int P, double **A, double **R, double *tau, double *tmp )
{
  int i, j, k;
  double norma2, norma;

#pragma scop

for(k = 0; k < P; k++){
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
   
   for(j = k+1; j < N; j++){
      tmp[j] = A[k][j];
      for(i = k+1; i < M; i++){
            tmp[j] += A[i][k] * A[i][j];
      }
      tmp[j] = tau[k] * tmp[j];
      A[k][j] = A[k][j] - tmp[j];
      for(i = k+1; i < M; i++){
         A[i][j] = A[i][j] - A[i][k] * tmp[j];
      }
   }
}


   for(i = 0; i < N; i++)
      for(j = 0; j < N; j++)
         R[i][j] = A[i][j];

   for(k = P-1; k > -1; k--){

      for(j = k+1; j < N; j++){

         tmp[j] = 0.e+00;

            for(i = k+1; i < M; i++){

               tmp[j] += A[i][k] * A[i][j];

            }

      }

      tmp[k] = 1.0e+00;

      for(j = k; j < N; j++){ 

         tmp[j] *= tau[k];

      }

      A[k][k] = 1.0e+00 - tmp[k];

      for(j = k+1; j < N; j++){

         A[k][j] = -tmp[j];

      }
      for(j = k+1; j < N; j++){

         for(i = k+1; i < M; i++){

            A[i][j] -= A[i][k] * tmp[j];

         }
      }
      
      for(i = k+1; i < M; i++){
         A[i][k] = - A[i][k] * tmp[k];
      }


   }


#pragma endscop

}
