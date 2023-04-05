#include <stdlib.h>
#include <math.h>
#include <cblas.h>

double check_orthog_blas(int m, int n, double *Q, int ldq )
{
   int i, j;
   double normRES;
   double *RES;

   RES = (double *) malloc( n * n * sizeof(double));
   for(i = 0; i < n; i++) for(j = 0; j < n; j++) if ( i==j ) RES[i+j*n] = 1.0e+00; else RES[i+j*n] = 0.0e+00; 

   cblas_dsyrk(CblasColMajor, CblasLower, CblasTrans, n, m, -1.0, Q, ldq, 1.0, RES, n);

   normRES = 0.0e+00;
   for(j = 0; j < n; j++){
      normRES += RES[i+j*n] * RES[i+j*n];
      for(i = j+1; i < n; i++){
         normRES += 2 * RES[i+j*n] * RES[i+j*n];
      }
   }
   free(RES);
   return sqrt( normRES );


}
