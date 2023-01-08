#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x

// This is in-place lower triangular inversion. In other words, A <- inv(A)
// where A is lower triangular. This code takes an n-by-n invertible lower
// triangular matrix A in input, and returns its inverse (in place). Note that
// the inverse of a lower triangular matrix is also lower triangular. There is
// also a numerical check. The upper part of the n-by-n array is not
// referenced. (So we assume lower triangularity. Anything can be stored in the
// upper part of A.)  The LAPACK equivalent subroutine is TRTRI.  The number of
// operations is n^3/3 FLOPS.

// there are six variants

int main(int argc, char ** argv) {

	int i, j, k, n, v;
	double **A, **B;
	double normR, tmp;

	srand(0);

    	n = 20;
	v = 3;

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

	if (( v != 1 ) && ( v != 2 ) && ( v != 3 ) && ( v != 4 ) && ( v != 5 ) && ( v != 6 )) return -1;

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

//	Make sure A is well conditioned by boosting its diagonal
 	for(i = 0; i < n; i++)
		A[i][i] += (double) n;

//	Save a copy of A in B (so that we can check)
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

/*************************************************************/
if ( v == 1 ){
//#pragma scop // variant 1
for ( i = 0 ; i < _PB_N ; i++ ){
  for ( j = 0 ; j < i ; j++ ){
    A[i][j] = - A[i][j] * A[j][j];
    for ( k = j+1 ; k < i ; k++ ){
      A[i][j] -= A[i][k] * A[k][j];
    }
    A[i][j] /= A[i][i];
  }
  A[i][i] = SCALAR_VAL(1.0) / A[i][i];
}
//#pragma endscop
}
if ( v == 2 ){
//#pragma scop // variant 2
for ( j = 0 ; j < _PB_N ; j++ ){
  for ( i = j+1 ; i < _PB_N ; i++ ){
    for ( k = j+1 ; k < i ; k++ ){
      A[i][j] += A[i][k] * A[k][j];
    }
    A[i][j] = - A[i][j] / A[i][i];
  }
  for ( i = j+1 ; i < _PB_N ; i++ ){
    A[i][j] /= A[j][j];
  }
  A[j][j] = SCALAR_VAL(1.0) / A[j][j];
}
//#pragma endscop
}
if ( v == 3 ){
//#pragma scop // variant 3
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
//#pragma endscop
}
if ( v == 4 ){
//#pragma scop // variant 4
for ( j = _PB_N-1 ; j > -1 ; j-- ){
  for ( i = _PB_N-1 ; i > j ; i-- ){
    A[i][j] *= A[i][i];
    for ( k	 = i-1 ; k > j ; k-- ){
      A[i][j] += A[i][k] * A[k][j];
    }
    A[i][j] = - A[i][j] / A[j][j];
  }
  A[j][j] = SCALAR_VAL(1.0) / A[j][j];
}
//#pragma endscop
}
if ( v == 5 ){
//#pragma scop // variant 5
for ( i = _PB_N-1 ; i > -1 ; i--){
  for ( j = i-1 ; j > -1 ; j-- ){
    for ( k	 = i-1 ; k > j ; k-- ){
      A[i][j] += A[i][k] * A[k][j];
    }
    A[i][j] = - A[i][j] / A[j][j];
  }
  for ( j = i-1 ; j > -1 ; j-- ){
    A[i][j] /= A[i][i];
  }
  A[i][i] = SCALAR_VAL(1.0) / A[i][i];
}
//#pragma endscop
}
if ( v == 6 ){
//#pragma scop // variant 6
for ( i = _PB_N-1 ; i > -1 ; i-- ){
  for ( j = i-1 ; j > -1 ; j-- ){
    A[i][j] = - A[i][j] / A[i][i];
  }
  for ( j = _PB_N-1 ; j > i ; j-- ){
    for ( k = i-1 ; k > -1 ; k-- ){ 
      A[j][k] += A[j][i] * A[i][k];
    }
  }
  for ( j = _PB_N-1 ; j > i ; j-- ){
    A[j][i] /= A[i][i];
  } 
  A[i][i] = SCALAR_VAL(1.0) / A[i][i];
}
//#pragma endscop
}
/*************************************************************/

//	check ||I - A*inv(A)||_F
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j <= i; j++) {
			tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
			for (k = j; k <= i; k++) {
				tmp -= B[i][k]*A[k][j];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	printf("[ TRTRI %d ] n = %4d; normR = %6.2e;\n", v, n, normR );

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
