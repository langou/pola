#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

// This routine takes a matrix pencil (A,B) (e.g., two square n-by-n matrices A
// and B) and reduces it to Hessenberg-Triangular form using orthogonal
// transformations.  The output of the subroutine is an upper Hessenberg matrix
// H, an upper triangular matrix T, an orthogonal matrix Q, and an orthogonal
// matrix Z such that 
//    A = Q * H * Z
//    B = Q * T * Z
//
// This subroutine is the first step to solve the nonsymmetric generalized
// eigenvalue problem. The equivalent LAPACK subroutines is DGGHD2. (The next
// step is the QZ algorithm.)
//
// If only the generalized eigenvalues are needed, (and not the generalized
// eigenvectors,) then only H and T need to be given to the QZ algorithm.
// (I.e. no need to compute the orthogonal factors Q and Z.)
//
// If we want to compute the the generalized eigenvalues and the generalized
// eigenvectors, then we need to give to the QZ algorithm the four matrices H,
// T, Q and Z.
//
// We perform six numerical checks: 
//  (1) || A * Z' - Q * H ||_f / || A ||_f is of the order of machine precision 
//  (2) || B * Z' - Q * T ||_f / || B ||_f is of the order of machine precision 
//  (3) || Z' * Z - I ||_f is of the order of machine precision 
//  (4) || Q' * Q - I ||_f is of the order of machine precision 
//  (5) H is upper Hessenberg 
//  (6) T is upper triangular
//
// This implementation replaces the input A by the output H, and the input B by
// the output T.
//
// Note 1: we have omitted the first step which is to perform a QR
// factorization on B such that B = Q * R.  We start with a B that is already
// triangular and set Q to identity. In a standard example, we cannot expect Q
// to be identity at the start of this algorithm. (It must be the Q-factor of
// the QR factorization of B.) The steps ommited are: [ Q, R ] = qr(B); B = R;
// A = Q'*A; Q = Q';
//
// Note 2: Storing the two pairs (c,s) needed for each (j,i) steps for an
// interval of 'm' (j,i) steps is an option. Then one can apply the stored
// rotations (defined by the pairs (c,s) as a block on Q and Z.  The present
// algorithm updates Q and Z at each (j,i) step of the algorithm.

int main(int argc, char ** argv) {

	int i, j, k, n;
	double **A, **B, **Z, **Q, **AA, **BB;
	double c, s, normA, normR, nrm, tmp;

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

	Z = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		Z[i] = (double *) malloc( n * sizeof(double));
	}

	Q = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		Q[i] = (double *) malloc( n * sizeof(double));
	}

	AA = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		AA[i] = (double *) malloc( n * sizeof(double));
	}

	BB = (double **) malloc( n * sizeof(double*));
 	for(i = 0; i < n; i++){
		BB[i] = (double *) malloc( n * sizeof(double));
	}

//	Create a random n-by-n matrix A.
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Save a copy of A in AA (so that we can check)
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			AA[i][j] = A[i][j];

//	Create a random triangular n-by-n matrix B.
 	for(i = 0; i < n; i++){
 		for(j = 0; j < i; j++)
			B[i][j] = 0.0e+00;
 		for(j = i; j < n; j++)
			B[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;
	}

//	Save a copy of B in BB (so that we can check)
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			BB[i][j] = B[i][j];

//	Set the Q matrix to the identity matrix
 	for(i = 0; i < n; i++){
 		for(j = 0; j < i; j++) Q[i][j] = 0.0e+00;
		Q[i][i] = 1.0e+00;
 		for(j = i+1; j < n; j++) Q[i][j] = 0.0e+00;
	}

//	Set the Z matrix to the identity matrix
 	for(i = 0; i < n; i++){
 		for(j = 0; j < i; j++) Z[i][j] = 0.0e+00;
		Z[i][i] = 1.0e+00;
 		for(j = i+1; j < n; j++) Z[i][j] = 0.0e+00;
	}

/*************************************************************/
//#pragma scop
  for (j = 0; j < _PB_N-2; j++) {
    for (i = _PB_N-2; i > j; i--) {
      nrm = SQRT_FUN ( A[i][j] * A[i][j] + A[i+1][j] * A[i+1][j] );  
      c = A[i][j] / nrm;  
      s = A[i+1][j] / nrm; 
      A[i][j] = nrm;
      A[i+1][j] = SCALAR_VAL(0.0);
      for (k = j+1; k < _PB_N; k++) {
        tmp = c * A[i][k] + s * A[i+1][k];
        A[i+1][k] = - s * A[i][k] + c * A[i+1][k];
        A[i][k] = tmp;
      }
      for (k = i; k < _PB_N; k++) {
        tmp = c * B[i][k] + s * B[i+1][k];
        B[i+1][k] = - s * B[i][k] + c * B[i+1][k];
        B[i][k] = tmp;
      }
      for (k = 0; k < _PB_N; k++) {
        tmp = c * Q[i][k] + s * Q[i+1][k];
        Q[i+1][k] = - s * Q[i][k] + c * Q[i+1][k];
        Q[i][k] = tmp;
      }
      nrm = SQRT_FUN ( B[i+1][i+1] * B[i+1][i+1] + B[i+1][i] * B[i+1][i] );  
      c = B[i+1][i+1] / nrm;  
      s = B[i+1][i] / nrm; 
      B[i+1][i+1] = nrm;
      B[i+1][i] = SCALAR_VAL(0.0);
      for (k = 0; k <= i; k++) {
        tmp = c * B[k][i] - s * B[k][i+1];
        B[k][i+1] = s * B[k][i] + c * B[k][i+1];
        B[k][i] = tmp;
      }
      for (k = 0; k < _PB_N; k++) {
        tmp = c * A[k][i] - s * A[k][i+1];
        A[k][i+1] = s * A[k][i] + c * A[k][i+1];
        A[k][i] = tmp;
      }
      for (k = i-j; k < _PB_N; k++) {
        tmp = c * Z[k][i] - s * Z[k][i+1];
        Z[k][i+1] = s * Z[k][i] + c * Z[k][i+1];
        Z[k][i] = tmp;
      }
    }
  }
//#pragma endscop
/*************************************************************/

//	check || Z' * Z - I ||_f
	printf("[ GGHD2 ] ");
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
			for (k = 0; k < n; k++) {
				tmp -= Z[k][i]*Z[k][j];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	printf("n = %4d; checks = [ %6.2e", n, normR );

//	check || Q' * Q - I ||_f
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp = ( i == j ) ? 1.0e+00 : 0.00e+00;
			for (k = 0; k < n; k++) {
				tmp -= Q[k][i]*Q[k][j];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	printf(" %6.2e", normR );

//	check || A * Z' - Q * H ||_f / || A ||_f
	normR = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp = 0.0e+00;
			for (k = 0; k < n; k++) {
				tmp += Q[i][k]*AA[k][j];
			}
			for (k = 0; k < n; k++) {
				tmp -= A[i][k]*Z[j][k];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	normA = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			normA += AA[i][j] * AA[i][j];
		}
	}
	normA = sqrt( normA );
	printf(" %6.2e", normR/normA );

//	check || B * Z' - Q * T ||_f / || B ||_f
	normR = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp = 0.0e+00;
			for (k = 0; k < n; k++) {
				tmp += Q[i][k]*BB[k][j];
			}
			for (k = 0; k < n; k++) {
				tmp -= B[i][k]*Z[j][k];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	normA = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			normA += BB[i][j] * BB[i][j];
		}
	}
	normA = sqrt( normA );
	printf(" %6.2e", normR/normA );

//	check H is upper Hessenberg 
	k = 0;
	for (i = 2; i < n; i++)
		for (j = 0; j < i-2; j++)
      			if( A[i][j] != 0 ) k++;
	printf(" %d", k );

//	check T is upper triangular
	k = 0;
	for (i = 1; i < n; i++)
		for (j = 0; j < i-1; j++)
      			if( B[i][j] != 0 ) k++;
	printf(" %d", k );

	printf(" ];\n" );

//	Free memory
 	for(i = 0; i < n; i++){
		free( BB[i] );
	}
	free( BB );

 	for(i = 0; i < n; i++){
		free( AA[i] );
	}
	free( AA );

 	for(i = 0; i < n; i++){
		free( Z[i] );
	}
	free( Z );

 	for(i = 0; i < n; i++){
		free( Q[i] );
	}
	free( Q );

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
