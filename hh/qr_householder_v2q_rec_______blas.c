#include <math.h>
#include <cblas.h>

void qr_householder_v2q_rec_______blas ( int m, int n, double *A, int lda, double *tau )
{
   int j, k, n1, n2;


   if( n == 1 ){

      A[0] = 1.0e+00 - tau[0];

      cblas_dscal( (m-1), -tau[0], &(A[1]), 1 );

   } else {

      n2 = n / 2;
      n1 = n - n2;

      qr_householder_v2q_rec_______blas ( m-n1, n2, &(A[n1+n1*lda]), lda, &(tau[n1]) );

      for(k = n1-1; k > -1; k--){

         cblas_dgemv( CblasColMajor, CblasTrans, (m-k-1), n2,
                      1.0, &(A[(k+1)+n1*lda]), lda, &(A[(k+1)+k*lda]), 1,
                      0.0, &(tau[n1]), 1 );

         for(j = n1; j < n; j++) tau[j] *= tau[k];


         for(j = n1; j < n; j++) A[k+lda*j] = -tau[j];

         cblas_dger( CblasColMajor, (m-k-1), n2,
                     -1.0, &(A[(k+1)+k*lda]), 1, &(tau[n1]), 1, &(A[(k+1)+n1*lda]), lda );

      }

      qr_householder_v2q_rec_______blas (m, n1, A, lda, tau);

   }

}
