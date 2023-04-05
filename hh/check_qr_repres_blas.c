#include <stdlib.h>
#include <math.h>
#include <cblas.h>

double check_qr_repres_blas(int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr )
{
   int i, j;
   double tmp, normRES, normA;
   double *QxR;

   QxR = (double *) malloc( m * n * sizeof(double));

   for(i = 0; i < m; i++) for(j = 0; j < n; j++) QxR[i+j*m] = Q[i+j*ldq];

   cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, 1.0, R, ldr, QxR, m);

   normRES = 0.0e+00;
   for(i = 0; i < m; i++){
      for(j = 0; j < n; j++){
         tmp = A[i+j*lda] - QxR[i+j*m];
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

   free(QxR);
   return normRES / normA;
}
