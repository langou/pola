#include <math.h>

double check_qr_repres(int M, int N, double A[M][N], double Q[M][N], double R[N][N] )
{
   int i, j, k;
   double tmp, normRES, normA;

   normRES = 0.0e+00;
   for(i = 0; i < M; i++){
      for(j = 0; j < N; j++){
         tmp = A[i][j];
         for(k = 0; k <= j; k++){
            tmp -= Q[i][k] * R[k][j];
         }
         normRES += tmp * tmp;
      }
   }
   normRES = sqrt( normRES );

   normA = 0.0e+00;
   for(i = 0; i < M; i++){
      for(j = 0; j < N; j++){
         normA += A[i][j] * A[i][j];
      }
   }
   normA = sqrt( normA );

   return normRES / normA;
}
