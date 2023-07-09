#include <math.h>
#include <cblas.h>

void qr_householder_v2q_rl__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau )
{
  int i, j, k, k0;
  double tmp;

#pragma scop

for(k0 = N-1; k0 > -1; k0-=B){

   for(j = k0+1; j < N; j++){

      for (k = k0; ((k > k0-B)&&(k > -1)); k--) {

         tmp = cblas_ddot( (M-k-1), &(A[(k+1)+k*lda]), 1, &(A[(k+1)+j*lda]), 1 );

         tmp *= tau[k];

         A[k+lda*j] = -tmp;

         cblas_daxpy( (M-k-1), -tmp, &(A[(k+1)+k*lda]), 1, &(A[(k+1)+j*lda]), 1 );

      }
   
   }

   for (k = k0; ((k > k0-B)&&(k > -1)); k--) {

      cblas_dgemv( CblasColMajor, CblasTrans, (M-k-1), k0-k,
                   1.0, &(A[(k+1)+(k+1)*lda]), lda, &(A[(k+1)+k*lda]), 1,
                   0.0, &(tau[k+1]), 1 );

      for(j = k+1; j < k0+1; j++) tau[j] *= tau[k];

      for(j = k+1; j < k0+1; j++) A[k+lda*j] = -tau[j];

      cblas_dger( CblasColMajor, (M-k-1), k0-k,
                  -1.0, &(A[(k+1)+k*lda]), 1, &(tau[k+1]), 1,
                  &(A[(k+1)+(k+1)*lda]), lda );
   
      A[k+lda*k] = 1.0e+00 - tau[k];

      cblas_dscal( (M-k-1), -tau[k], &(A[(k+1)+k*lda]), 1 );

   }

}
#pragma endscop
}
