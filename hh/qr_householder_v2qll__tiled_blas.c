#include <math.h>
#include <cblas.h>

void qr_householder_v2qll__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau )
{
   int j, k, k0;
   double tmp;

#pragma scop

for(k0 = N-1; k0 > -1; k0-=B){

   for (k = k0; ((k > k0-B)&&(k > -1)); k--) {

      A[k+lda*k] = 1.0e+00 - tau[k];

      cblas_dscal( (M-k-1), -tau[k], &(A[(k+1)+k*lda]), 1 );

      for(j = k-1; j > k0-B; j--){

         tmp = cblas_ddot( (M-j-1), &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );

         tmp *= tau[j];

         A[j+lda*k] = -tmp;

         cblas_daxpy( (M-j-1), -tmp, &(A[(j+1)+j*lda]), 1, &(A[(j+1)+k*lda]), 1 );

      }

   }

   for(j = k0-B; j > -1; j--){

      int l; if ( k0-B > -1 ) l = B; else l = k0+1;

      cblas_dgemv( CblasColMajor, CblasTrans, (M-j-1), l,
                   1.0, &(A[(j+1)+(k0-l+1)*lda]), lda, &(A[(j+1)+j*lda]), 1,
                   0.0, &(tau[k0-l+1]), 1 );

      for (k = k0; k > k0 - l; k--) tau[k] *= tau[j];

      for (k = k0; k > k0 - l; k--) A[j+lda*k] = -tau[k];

      cblas_dger( CblasColMajor, (M-j-1), l,
                  -1.0, &(A[(j+1)+j*lda]), 1, &(tau[k0-l+1]), 1,
                  &(A[(j+1)+(k0-l+1)*lda]), lda );

   }

}

#pragma endscop
}
