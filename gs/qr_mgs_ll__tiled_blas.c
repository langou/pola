#include <math.h>
#include <mkl_cblas.h>

// qr_mgs_ll__tiled
// this is the left-looking variant of MGS

void qr_mgs_ll__tiled_blas (int m, int n, int b, double *A, int lda, double *R, int ldr )
{
int i, j;
int j0;
for (j0 = 0; j0 < n; j0+=b) {
// read A(1:M,j0:j0+B)
   int l; if ( j0+b < n ) l = b; else l = n-j0;
   for (i = 0; i < j0; i++) {
//    read A(1:M,i)

      cblas_dgemv( CblasColMajor, CblasTrans, m, l,
                  1.0, &(A[j0*lda]), lda, &(A[i*lda]), 1,
                  0.0, &(R[i+j0*ldr]), ldr );

      cblas_dger( CblasColMajor, m, l,
                -1.0, &(A[i*lda]), 1, &(R[i+j0*ldr]), ldr, &(A[j0*lda]), lda );

//    discard A(1:M,i)
   }
   for (j = j0; ((j < j0+b)&&(j < n)); j++) {
      for (i = j0; i < j; i++) {
         R[i+j*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[j*lda]), 1 );
         cblas_daxpy( m, -R[i+j*ldr], &(A[i*lda]), 1, &(A[j*lda]), 1 );
      }
      R[j+j*ldr] = cblas_ddot( m, &(A[j*lda]), 1, &(A[j*lda]), 1 );
      R[j+j*ldr] = sqrt(R[j+j*ldr]);
      cblas_dscal( m, 1.0/R[j+j*ldr], &(A[j*lda]), 1 );
   }
// write A(1:M,j0:j0+B)
}
}
