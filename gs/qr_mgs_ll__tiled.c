#include <math.h>

// qr_mgs_ll__tiled
// this is the left-looking variant of MGS

void qr_mgs_ll__tiled (int m, int n, int b, double *A, int lda, double *R, int ldr )
{
int i, j, k;
int j0;
#pragma scop
for (j0 = 0; j0 < n; j0+=b) {
// read A(1:M,j0:j0+B)
   for (i = 0; i < j0; i++) {
//    read A(1:M,i)
      for (j = j0; ((j < j0+b)&&(j < n)); j++) {
         R[i+j*ldr] = 0.0e+00;
         for (k = 0; k < m; k++)
            R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];
         for (k = 0; k < m; k++)
            A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];
      }
//    discard A(1:M,i)
   }
   for (j = j0; ((j < j0+b)&&(j < n)); j++) {
      for (i = j0; i < j; i++) {
         R[i+j*ldr] = 0.0e+00;
         for (k = 0; k < m; k++)
            R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];
         for (k = 0; k < m; k++)
            A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];
      }
      R[j+j*ldr] = 0.0e+00;
      for (k = 0; k < m; k++)
         R[j+j*ldr] += A[k+j*lda] * A[k+j*lda];
      R[j+j*ldr] = sqrt(R[j+j*ldr]);
      for (k = 0; k < m; k++)
         A[k+j*lda] /= R[j+j*ldr];
   }
// write A(1:M,j0:j0+B)
}
#pragma endscop
}
