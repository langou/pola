#include <math.h>

double check_orthog(int m, int n, double *Q, int ldq )
{
   int i, j, k;
   double tmp, normRES;

   normRES = 0.0e+00;
   for(i = 0; i < n; i++){
      for(j = 0; j < n; j++){
         if ( i==j ) tmp = 1.0e+00; else tmp = 0.0e+00; 
         for(k = 0; k < m; k++){
            tmp -= Q[k+i*ldq] * Q[k+j*ldq];
         }
         normRES += tmp * tmp;
      }
   }
   return sqrt( normRES );

}
