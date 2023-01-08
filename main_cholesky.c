#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SQRT_FUN sqrt

// This is in-place Cholesky factorization. This code takes an n-by-n symmetric
// positive definite matrix A in input and returns its lower triangular factor
// C in output. C is lower triangular such that A = C * C'. In input, array A
// represents the symmetric positive definite matrix. In output, array A
// represents the lower trianglar factor C. We only reference the lower part of
// the n-by-n array A.  The upper part of the n-by-n array is not referenced.
// So we assume symmetry in input and triangularity in output. Anything can be
// stored in the upper part of A.  There is also a numerical check.  The LAPACK
// equivalent subroutine is POTRF. The number of operations is n^3/3 FLOPS.

// There are three variants.
// Variant 1 is the variant in polybench 4.2.1.
// Variant 1 is also called the bordered variant.
// Variant 2 is also called the left looking variant.
// Variant 3 is also called the right looking variant.

int main(int argc, char ** argv) {

	int i, j, k, n, v;
	double **A, **B;
	double normA, normR, tmp;

	srand(0);

    n = 20;
	v = 1;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}
	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-v") == 0) {
			v  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	if (( v != 1 ) && ( v != 2 ) && ( v != 3 ) ) return -1;

	A = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		A[i] = (double *) malloc( n * sizeof(double));
	}

	B = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		B[i] = (double *) malloc( n * sizeof(double));
	}

//	Create a random matrix A.  We create a nonsymmetric matrix but will
//	only reference the lower part and assume symmetry
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Make sure A is positive definite by boosting its diagonal
 	for(i = 0; i < n; i++)
		A[i][i] += (double) n;

//	Save a copy of A in B
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

/*************************************************************/
if ( v == 1 ){
//#pragma scop // variant 1
for ( i = 0 ; i < _PB_N ; i++ ){
  for ( j = 0 ; j < i ; j++ ){
    for ( k = 0 ; k < j ; k++ ){
      A[i][j] -= A[i][k] * A[j][k];
    }
    A[i][j] /= A[j][j];
  }
  for ( j = 0 ; j < i ; j++ ){
    A[i][i] -= A[i][j] * A[i][j];
  }
  A[i][i] = SQRT_FUN(A[i][i]);
}
//#pragma endscop
}
if ( v == 2 ){
//#pragma scop // variant 2
for ( i = 0 ; i < _PB_N ; i++ ){
  for ( j = 0 ; j < i ; j++ ){
    A[i][i] -= A[i][j] * A[i][j];
  }
  A[i][i] = SQRT_FUN(A[i][i]);
  for ( j = 0 ; j < i ; j++ ){
    for ( k = i+1 ; k < _PB_N ; k++ ){
      A[k][i] -= A[k][j] * A[i][j];
    }
  }
  for ( j = i+1 ; j < _PB_N ; j++ ){
    A[j][i] /= A[i][i];
  }
}
//#pragma endscop
}
if ( v == 3 ){
//#pragma scop // variant 3
for ( i = 0 ; i < _PB_N ; i++ ){
  A[i][i] = SQRT_FUN(A[i][i]);
  for ( j = i+1 ; j < _PB_N ; j++ ){
    A[j][i] /= A[i][i];
  }
  for ( j = i+1 ; j < _PB_N ; j++ ){
    for ( k = i+1 ; k < j ; k++ ){
      A[j][k] -= A[j][i] * A[k][i];
    }
    A[j][j] -= A[j][i] * A[j][i];
  }
}
//#pragma endscop
}
/*************************************************************/

//	check ||A - LLt||_F / ||A||_F
	normR = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			tmp = B[i][j];
			for (k = 0; k <= j; k++) {
				tmp -= A[i][k]*A[j][k];
			}
			normR += ( ( i == j ) ? 1.0e+00 : 2.00e+00 ) * tmp * tmp;
		}
	}
	normR = sqrt( normR );
	normA = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			tmp = B[i][j];
			normA += ( ( i == j ) ? 1.0e+00 : 2.00e+00 ) * tmp * tmp;
		}
	}
	normA = sqrt( normA );
	printf("[ POTRF %d ] n = %4d; normR / normA = %6.2e;\n", v, n, normR/normA);

//	Free memory
 	for(i = 0; i < n; i++){
		free( B[i] );
	}
	free( B );

 	for(i = 0; i < n; i++){
		free( A[i] );
	}
	free( A );

	return 0;
}
