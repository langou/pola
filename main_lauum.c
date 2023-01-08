#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x

// This is in-place matrix-matrix mutiplication of the kind A <- A'*A where in
// input A is lower triangular matrix, and in output A is a symmetric matrix.
// There is also a numerical check. The upper part of the n-by-n array is not
// referenced. (So we assume lower triangularity in input, and symmetry from
// the lower part in output. Anything can be stored in the upper part of A.)
// The LAPACK equivalent subroutine is LAUUM.  The number of operations is
// n^3/3 FLOPS.

// There are three variants.

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

//	Create a random matrix A.  We create a dense matrix but will only
//	reference the lower part and assume triangularity
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Save a copy of A in B
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

/*************************************************************/
if ( v == 1 ){
//#pragma scop // variant 1
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
}
if ( v == 2 ){
//#pragma scop // variant 2
for ( i = 0 ; i < _PB_N ; i++ ){
  for ( j = 0 ; j < i ; j++ ){
    A[i][j] *= A[i][i];
  }
  for ( j = 0 ; j < i ; j++ ){
    for ( k = i+1 ; k < _PB_N ; k++ ){
      A[i][j] += A[k][i] * A[k][j];
    }
  }
  A[i][i] = A[i][i] * A[i][i];
  for ( j = i+1 ; j < _PB_N ; j++ ){
    A[i][i] += A[j][i] * A[j][i];
  }
}
//#pragma endscop
}
if ( v == 3 ){
//#pragma scop // variant 3
for ( i = 0 ; i < _PB_N ; i++ ){
  A[i][i] = A[i][i] * A[i][i];
  for ( j = i+1 ; j < _PB_N ; j++ ){
    A[i][i] += A[j][i] * A[j][i];
  }
  for ( j = i+1 ; j < _PB_N ; j++ ){
    A[j][i] *= A[j][j];
    for ( k = j+1 ; k < _PB_N ; k++ ){
      A[j][i] += A[k][j] * A[k][i];
    }
  }
}
//#pragma endscop
}
/*************************************************************/

//	check ||lauum(A) - AAt||_F / ||AAt||_F
	normR = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = j; i < n; i++) {
			tmp = A[i][j];
			for (k = i; k < n; k++) {
				tmp -= B[k][i]*B[k][j];
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
	printf("[ LAUUM %d ] n = %4d; normR / normA = %6.2e;\n", v, n, normR/normA);

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
