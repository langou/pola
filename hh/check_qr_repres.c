#include <math.h>

double check_qr_repres(int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr )
{
   int i, j, k;
   double tmp, normRES, normA;

   normRES = 0.0e+00;
   for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
         tmp = A[i+j*lda];
         for(k = 0; k <= j; k++){
            tmp -= Q[i+k*ldq] * R[k+j*ldr];
         }
         normRES += tmp * tmp;
      }
   }
   normRES = sqrt( normRES );

   normA = 0.0e+00;
   for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
         normA += A[i+j*lda] * A[i+j*lda];
      }
   }
   normA = sqrt( normA );

   return normRES / normA;
}
