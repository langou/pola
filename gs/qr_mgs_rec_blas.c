#include <math.h>
#include <mkl_cblas.h>

// qr_mgs_ll
// this is the left-looking variant of MGS
extern void qr_mgs_ll (int m, int n, double *A, int lda, double *R, int ldr );

void qr_mgs_rec_blas (int m, int n, double *A, int lda, double *R, int ldr )
{
int n1, n2;
int i;

if( n == 1 ){

   (*R) = cblas_ddot( m, A, 1, A, 1 );
   (*R) = sqrt((*R));
   cblas_dscal( m, 1.0/(*R), A, 1 );

} else {

   n1 = n / 2;
   n2 = n - n1;

   qr_mgs_rec_blas ( m, n1, A, lda, R, ldr );


// for (i = 0; i < n1; i++) {
//    for (int j = n1; j < n; j++) {
//       R[i+j*ldr] = 0.0e+00;
//       for (k = 0; k < m; k++)
//          R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];
//       for (k = 0; k < m; k++)           
//          A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];
//    }
// }

// for (i = 0; i < n1; i++) {
//    for (int j = n1; j < n; j++) {
//       R[i+j*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[j*lda]), 1 );
//       cblas_daxpy( m, -R[i+j*ldr], &(A[i*lda]), 1, &(A[j*lda]), 1 );
//    }
// }

   for (i = 0; i < n1; i++) {
      cblas_dgemv( CblasColMajor, CblasTrans, m, n2,
                  1.0, &(A[n1*lda]), lda, &(A[i*lda]), 1,
                  0.0, &(R[i+n1*ldr]), ldr );

      cblas_dger( CblasColMajor, m, n2,
                -1.0, &(A[i*lda]), 1, &(R[i+n1*ldr]), ldr, &(A[n1*lda]), lda );

   }

   qr_mgs_rec_blas ( m, n2, &(A[n1*lda]), lda, &(R[n1+n1*ldr]), ldr );

}

}
