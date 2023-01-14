#include <math.h>

double check_orth(int M, int N, double Q[M][N] )
{
   int i, j, k;
   double tmp, normRES;

   normRES = 0.0e+00;
   for(i = 0; i < N; i++){
      for(j = 0; j < N; j++){
         if ( i==j ) tmp = 1.0e+00; else tmp = 0.0e+00; 
         for(k = 0; k < M; k++){
            tmp -= Q[k][i] * Q[k][j];
         }
         normRES += tmp * tmp;
      }
   }
   return sqrt( normRES );

}
