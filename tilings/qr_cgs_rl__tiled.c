#include <math.h>

// qr_cgs_rl
//
// this is the right-looking (RL) version of Classical Gram-Schmidt (CGS)
//
// array A is input only
// arrays Q and R do not need to be initialized in input, they are ouput only
// the strictly lower part of R is not referenced

void qr_cgs_rl__tiled (int M, int N, int B, double A[M][N], double Q[M][N], double R[N][N] )
{
int i, j, k;
int iB =B;
int kB1=B;
int kB2=B;
int jB1=B;
int jB2=B;
int i0, j0, k0;
#pragma scop
for(i = 0; i < M; i++) for(j = 0; j < N; j++) Q[i][j] = A[i][j];

// this zeroing out can be done right before the first useage of R[i][j]
// so as to save N^2/2 read and writes. They are moved early on in an
// attempt to make the code more readable.
for(i = 0; i < N; i++) for(j = 0; j < N; j++) R[i][j] = 0e+00;

for (i0 = 0; i0 < N; i0+=iB) {

// this is the panel factorization
// note that there are two matrices here A and Q, not just one
// note that the two panels might not fit in cache but the number of 
// operations is 4 M B^2 so the IO is at most O(MB^2)
// read A(1:M,i0:i0+B-1)
// read Q(1:M,i0:i0+B-1)
   for (i = i0; ((i < i0+iB)&&(i < N)); i++) {

      for (k = 0; k < M; k++)
         R[i][i] += Q[k][i] * Q[k][i];
      R[i][i] = sqrt(R[i][i]);
      for (k = 0; k < M; k++)
         Q[k][i] /= R[i][i];

      for (j = i + 1; ((j < i0+iB)&&(j < N)); j++)
         for (k = 0; k < M; k++)
            R[i][j] += Q[k][i] * A[k][j];

      for (j = i + 1; ((j < i0+iB)&&(j < N)); j++)
         for (k = 0; k < M; k++)
            Q[k][j] -= Q[k][i] * R[i][j];

   }
// discard A(1:M,i0:i0+B-1)
// write Q(1:M,i0:i0+B-1)

   for (j0 = i0+iB; j0 < N; j0+=jB1)
//    ``read`` R[i0:i0+iB-1,j0:j0+jB-1]: size iB x jB1
      for (k0 = 0; k0 < M; k0+=kB1)
         for (k = k0; ((k < k0+kB1)&&(k < M)); k++)
            for (i = i0; ((i < i0+iB)&&(i < N)); i++)
               for (j = j0; ((j < j0+jB1)&&(j < N)); j++)
//                read Q[k,i]
//                read A[k,j]
                  R[i][j] += Q[k][i] * A[k][j];
//    write R[i0:i0+iB-1,j0:j0+jB-1]: size iB x jB1

      for (j0 = i0+iB; j0 < N; j0+=jB2)
         for (k0 = 0; k0 < M; k0+=kB2)
//          read kB2 x jB2
//          read kB2 x iB
//          read iB x jB2
            for (i = i0; ((i < i0+iB)&&(i < N)); i++)
               for (j = j0; ((j < j0+jB2)&&(j < N)); j++)
                  for (k = k0; ((k < k0+kB2)&&(k < M)); k++)
                     Q[k][j] -= Q[k][i] * R[i][j];
//          write kB2 x jB2

}
#pragma endscop
}
