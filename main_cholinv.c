#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

// This is in-place Cholesky inversion. This code takes an n-by-n symmetric
// positive definite matrix A in input and returns its inverse. Note that the
// inverse of a symmetric positive definite is also symmetric positive
// definite. We only reference the lower part of the n-by-n array A.  The upper
// part of the n-by-n array is not referenced. So we assume symmetry in input
// and output. Anything can be stored in the upper part of A.  There is also a
// numerical check.  The LAPACK equivalent subroutines are POTRF+TRTRI+LAUUM.
// The number of operations is n^3 FLOPS.

// One can use any variant for POTRF, any variant for TRTRI, and any variant
// for LAUUM. Here we use POTRF_v3, TRTRI_v3 and LAUUM_v1.

int main(int argc, char ** argv) {

	int i, j, k, n;
	double **A, **B;
	double normR, tmp;

	srand(0);

    	n = 20;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	A = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		A[i] = (double *) malloc( n * sizeof(double));
	}

	B = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		B[i] = (double *) malloc( n * sizeof(double));
	}

//	Create a random matrix A
//	We create a dense matrix but will only reference the lower part and assume symmetry
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
//#pragma scop
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
for ( j = 0 ; j < _PB_N ; j++ ){
  for ( i = j+1 ; i < _PB_N ; i++ ){
    A[i][j] = - A[i][j] / A[j][j];
  }
  for ( i = j+1 ; i < _PB_N ; i++ ){
    for ( k = 0 ; k < j ; k++ ){
      A[i][k] += A[i][j] * A[j][k];
    }
  }
  for ( i = 0 ; i < j ; i++ ){
    A[j][i] /= A[j][j];
  }
  A[j][j] = SCALAR_VAL(1.0) / A[j][j];
}
for ( i = 0 ; i < _PB_N ; i++ ){
  for ( j = 0 ; j < i ; j++ ){
    A[j][j] += A[i][j] * A[i][j];
    for ( k = j+1 ; k < i ; k++ ){
      A[k][j] += A[i][k] * A[i][j];
    }
  }
  for ( j = 0 ; j < i ; j++ ){
    A[i][j] *= A[i][i];
  }
  A[i][i] = A[i][i] * A[i][i];
}
//#pragma endscop
/*************************************************************/

//	check ||I - A*inv(A)||_F
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < i; j++) {
			tmp = 0.0e+00;
			for (k = 0; k <= j; k++) {
				tmp -= B[i][k]*A[j][k];
			}
			for (k = j+1; k <= i; k++) {
				tmp -= B[i][k]*A[k][j];
			}
			for (k = i+1; k < n ; k++) {
				tmp -= B[k][i]*A[k][j];
			}
			normR += tmp * tmp;
		}
		tmp = 1.0e+00;
		for (k = 0; k <= i; k++) {
			tmp -= B[i][k]*A[i][k];
		}
		for (k = i+1; k < n ; k++) {
			tmp -= B[k][i]*A[k][i];
		}
		normR += tmp * tmp;
		for (j = i+1; j < n; j++) {
			tmp = 0.0e+00;
			for (k = 0; k <= i; k++) {
				tmp -= B[i][k]*A[j][k];
			}
			for (k = i+1; k <= j; k++) {
				tmp -= B[k][i]*A[j][k];
			}
			for (k = j+1; k < n ; k++) {
				tmp -= B[k][i]*A[k][j];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	printf("[ CHOLINV ] n = %4d; normR = %6.2e;\n", n, normR );

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
