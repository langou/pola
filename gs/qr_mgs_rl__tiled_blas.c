#include <math.h>
#include <cblas.h>

// qr_mgs_rl

void qr_mgs_rl__tiled_blas (int m, int n, int b, double *A, int lda, double *R, int ldr )
{
int i, j, k;
int i0;

// the coefficients of R are computed, used only once and right after
// their computation, so we write them as soon as we are done with them
// we write the elements of R on the fly, and so 
// (1) the total number of writes due to R is N * (N+1) / 2
// (2) the requirement of cache for R is 1 for the whole computation, so
// we do not take R at all for the cache useage.

#pragma scop
for (i0 = 0; i0 < n; i0+=b) {
////////////////////////////////////////////////////////////
// read A(1:M,i0:i0+B-1)
   for (i = i0; ((i < i0+b)&&(i < n)); i++) {
      R[i+i*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[i*lda]), 1 );
      R[i+i*ldr] = sqrt(R[i+i*ldr]);
      cblas_dscal( m, 1.0/R[i+i*ldr], &(A[i*lda]), 1 );
      for (j = i + 1; ((j < i0+b)&&(j < n)); j++) {
         R[i+j*ldr] = cblas_ddot( m, &(A[i*lda]), 1, &(A[j*lda]), 1 );
         cblas_daxpy( m, -R[i+j*ldr], &(A[i*lda]), 1, &(A[j*lda]), 1 );
      }
   }
////////////////////////////////////////////////////////////
   for (j = i0+b; j < n; j++) {
//    read A(1:M,j)
      for (i = i0; ((i < i0+b)&&(i < n)); i++) {
         R[i+j*ldr] = 0.0e+00;
         for (k = 0; k < m; k++)
            R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];
         for (k = 0; k < m; k++)
            A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];
      }
//    write A(1:M,j)
   }
// write A(1:M,i0:i0+B-1)
}
#pragma endscop
}
