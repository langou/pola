#include <math.h>
#include <cblas.h>

void qr_householder_a2v_rec_______blas ( int m, int n, double *A, int lda, double *tau )
{
   int j, k, n1, n2;
   double norma2, norma;

   if( n == 1 ){

       norma2 = cblas_ddot( (m-1), &(A[1]), 1, &(A[1]), 1 );
       norma = sqrt( A[0] * A[0] + norma2 );
       
       A[0] = ( A[0] > 0 ) ? ( A[0] + norma ) : ( A[0] - norma ) ;
       
       tau[0] = 2.0 / ( 1.0 + norma2 / ( A[0] * A[0] ) ) ;
       
       cblas_dscal( (m-1), 1.0/A[0], &(A[1]), 1 );
       A[0]= ( A[0] > 0 ) ? ( - norma ) : ( norma ) ;

   } else {

      n1 = n / 2;
      n2 = n - n1;

      qr_householder_a2v_rec_______blas ( m, n1, A, lda, tau );

      for(k = 0; k < n1; k++){

         for(j = n1; j < n; j++) tau[j] = A[k+lda*j];
   
         cblas_dgemv( CblasColMajor, CblasTrans, (m-k-1), n2,
                      1.0, &(A[(k+1)+n1*lda]), lda, &(A[(k+1)+k*lda]), 1,
                      1.0, &(tau[n1]), 1 );

         for(j = n1; j < n; j++) tau[j] = tau[k] * tau[j];
         for(j = n1; j < n; j++) A[k+lda*j] = A[k+lda*j] - tau[j];

         cblas_dger( CblasColMajor, (m-k-1), n2,
                     -1.0, &(A[(k+1)+k*lda]), 1, &(tau[n1]), 1, &(A[(k+1)+n1*lda]), lda );

      }

      qr_householder_a2v_rec_______blas ( m-n1, n2, &(A[n1+n1*lda]), lda, &(tau[n1]) );

   }

}
