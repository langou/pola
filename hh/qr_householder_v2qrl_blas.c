#include <math.h>
#include <cblas.h>

void qr_householder_v2qrl_blas ( int M, int N, double *A, int lda, double *tau )
{
  int j, k;

#pragma scop

for(k = N-1; k > -1; k--){

   cblas_dgemv( CblasColMajor, CblasTrans, (M-k-1), (N-k-1),
                1.0, &(A[(k+1)+(k+1)*lda]), lda, &(A[(k+1)+k*lda]), 1,
                0.0, &(tau[k+1]), 1 );


   for(j = k+1; j < N; j++) tau[j] *= tau[k];


   for(j = k+1; j < N; j++) A[k+lda*j] = -tau[j];

   cblas_dger( CblasColMajor, (M-k-1), (N-k-1),
      -1.0, &(A[(k+1)+k*lda]), 1, &(tau[k+1]), 1, &(A[(k+1)+(k+1)*lda]), lda );

   A[k+lda*k] = 1.0e+00 - tau[k];

   cblas_dscal( (M-k-1), -tau[k], &(A[(k+1)+k*lda]), 1 );

}
#pragma endscop
}
