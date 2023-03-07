#include <math.h>

// qr_mgs_ll
// this is the left-looking variant of MGS

void qr_mgs_ll (int m, int n, double *A, int lda, double *R, int ldr )
{
int i, j, k;
#pragma scop
for (j = 0; j < n; j++) {
// read A(1:M,j)
   for (i = 0; i < j; i++) {
//    read A(1:M,i)
      R[i+j*ldr] = 0.0e+00;                // LL2  RL6
      for (k = 0; k < m; k++)
         R[i+j*ldr] += A[k+i*lda] * A[k+j*lda];  // LL3  RL7
      for (k = 0; k < m; k++)           
         A[k+j*lda] -= A[k+i*lda] * R[i+j*ldr];  // LL4  RL8
//    write R[i][j]
//    discard A(1:M,i)
   }                                    
   R[j+j*ldr] = 0.0e+00;                   // LL5  RL2   
   for (k = 0; k < m; k++)              
      R[j+j*ldr] += A[k+j*lda] * A[k+j*lda];     // LL6  RL3    
   R[j+j*ldr] = sqrt(R[j+j*ldr]);             // LL7  RL4
   for (k = 0; k < m; k++)
      A[k+j*lda] /= R[j+j*ldr];               // LL8  RL5
// write R[j][j]
// write A(1:M,j)
}
#pragma endscop
}
