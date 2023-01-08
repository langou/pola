#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

// This routine takes an n-by-n matrix A and reduces it to Hessenberg form
// using an orthogonal similarity transformation. The output of the subroutine
// is an upper Hessenberg matrix H and an orthogonal matrix Q such that 
//      A = Q * H * Q^T
//
// This subroutine is the first step to solve the nonsymmetric eigenvalue
// problem. The equivalent LAPACK subroutines is DGEHD2. (The next step is the
// QR algorithm.)
//
// If only the eigenvalues are needed, (and not the eigenvectors,) then only H
// needs to be given to the QR algorithm.  (I.e. no need to compute the
// orthogonal factor Q.)
//
// If we want to compute the the eigenvalues and the eigenvectors, then we need
// to give to the QR algorithm the two matrices H and Q .
//
// We use Householder transformation. This is a level-2 BLAS subroutine.
//
// This implementation requires a workspace of size n. (`tmp` in the code
// below.)
//
// We perform three numerical checks: 
//  (1) || A * Q' - Q * H ||_f / || A ||_f is of the order of machine precision 
//  (2) || Q' * Q - I ||_f is of the order of machine precision 
//  (3) H is upper Hessenberg 
//
// This implementation replaces the input A by the output H and the Householder
// vector matrix V stored beneath H. This is akin to LAPACK DGEHD2.
//
// The matrix Q is not computed by DGEHD2. We do the same and we post-process
// with a sequence of operations similar to DORGHR to explicitly compute Q for
// the check.  (For the check, we could also have applied Q implicitly instead
// of explicitly computing Q.)


int main(int argc, char ** argv) {

	int i, j, k, n;
	double *AAA, *BBB, *QQQ;
	double **A, **B, **Q;
	double *tmp;
	double tau, norma, norma2, normA, normR;

	srand(0);

    n = 20;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	AAA = (double *) malloc( n * n * sizeof(double));
	A = (double **) malloc( n * sizeof(double*) );
 	for(i = 0; i < n; i++){
		A[i] = AAA + i * n ;
	}

	BBB = (double *) malloc( n * n * sizeof(double) );
	B = (double **) malloc( n * sizeof(double*) );
 	for(i = 0; i < n; i++){
		B[i] = BBB + i * n ;
	}

	tmp = (double *) malloc( n * sizeof(double) );

//	Create a random n-by-n matrix A.
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Save a copy of A in B (so that we can check)
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

/*************************************************************/
//#pragma scop
   for ( j = 0; j < _PB_N-2; j++ ) {
      norma2 = SCALAR_VAL(0.0) ;
      for ( i = j+2; i < _PB_N; i++ ) {
         norma2 += A[i][j] * A[i][j] ;
      }
      norma = SQRT_FUN ( A[j+1][j] * A[j+1][j] + norma2 ) ;
      A[j+1][j] = ( A[j+1][j] > 0 ) ? ( A[j+1][j] + norma ) : ( A[j+1][j] - norma ) ;
      tau = SCALAR_VAL(2.0) / ( SCALAR_VAL(1.0) + norma2 / ( A[j+1][j] * A[j+1][j] ) ) ;
      for ( i = j+2; i < _PB_N; i++ ) {
         A[i][j] /= A[j+1][j] ;
      }
      A[j+1][j] = ( A[j+1][j] > 0 ) ? ( - norma ) : ( + norma ) ;
      for (i = j+1; i < _PB_N; i++) {
         tmp[i] = A[j+1][i] ;
         for (k = j+2; k < _PB_N; k++) {
            tmp[i] += A[k][j] * A[k][i];
         }
      }
      for (i = j+1; i < _PB_N; i++) {
         tmp[i] *= tau ;
      }
      for (i = j+1 ; i < _PB_N ; i++) {
         A[j+1][i] -= tmp[i] ;
      }
      for (i = j+2; i < _PB_N; i++) {
         for (k = j+1; k < _PB_N; k++) {
            A[i][k] -= A[i][j] * tmp[k];
         }
      }
      for (i = 0; i < _PB_N; i++) {
         tmp[i] = A[i][j+1];
         for (k = j+2; k < _PB_N; k++) {
            tmp[i] += A[i][k] * A[k][j];
         }
      }
      for (i = 0; i < _PB_N; i++) {
         tmp[i] *= tau ;
      }
      for (i = 0; i < _PB_N; i++) {
         A[i][j+1] -= tmp[i] ;
      }
      for (i = 0; i < _PB_N; i++) {
         for (k = j+2; k < _PB_N; k++) {
            A[i][k] -= tmp[i] * A[k][j] ;
         }
      }
   }
//#pragma endscop
/*************************************************************/

	QQQ = (double *) malloc( n * n * sizeof(double));
	Q = (double **) malloc( n * sizeof(double*) );
 	for(i = 0; i < n; i++){
		Q[i] = QQQ + i * n ;
	}

// the piece of code below is nothing else than LAPACK DORGHD
// it constructs the matrix Q
   for(i = 0; i < n; i++){
      for(j = 0; j < i; j++) Q[i][j] = 0.0e+00;
      Q[i][i] = 1.0e+00;
      for(j = i+1; j < n; j++) Q[i][j] = 0.0e+00;
   }
   for ( j = _PB_N-3; j > -1; j-- ) {
      tau = SCALAR_VAL(1.0) ;
      for (i = j+2; i < _PB_N; i++) {
         tau += A[i][j] * A[i][j] ;
      }
      tau = 2 / tau ;
      for (i = j+1; i < _PB_N; i++) {
         tmp[i] = Q[j+1][i] ;
         for (k = j+2; k < _PB_N; k++) {
            tmp[i] += A[k][j] * Q[k][i] ;
         }
      }
      for (i = j+1; i < _PB_N; i++) {
         tmp[i] *= tau ;
      }
      for (i = j+1; i < _PB_N; i++) {
         Q[j+1][i] -= tmp[i] ;
      }
      for (i = j+2; i < _PB_N; i++) {
         for (k = j+1; k < _PB_N; k++) {
            Q[i][k] -= A[i][j] * tmp[k];
         }
      }
   }

//	check || Q' * Q - I ||_f
	printf("[ GEHD2 ] n = %4d; checks = [ ", n );
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp[0] = ( i == j ) ? 1.0e+00 : 0.00e+00;
			for (k = 0; k < n; k++) {
				tmp[0] -= Q[k][i]*Q[k][j];
			}
			normR += tmp[0] * tmp[0];
		}
	}
	normR = sqrt( normR );
	printf(" %6.2e", normR );

// check || A * Q' - Q * H ||_f / || A ||_f
	normR = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp[0] = 0.0e+00;
			for (k = 0; (k < n)&&(k < j+2); k++) {
				tmp[0] += Q[i][k]*A[k][j];
			}
			for (k = 0; k < n; k++) {
				tmp[0] -= B[i][k]*Q[k][j];
			}
			normR += tmp[0] * tmp[0] ;
		}
	}
	normR = sqrt( normR );
	normA = 0.0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			normA += B[i][j] * B[i][j] ;
		}
	}
	normA = sqrt( normA );
	printf(" %6.2e", normR / normA );

	printf(" ];\n");

//	Free memory
	free( Q );
	free( QQQ );

	free( tmp );

	free( B );
	free( BBB );

	free( A );
	free( AAA );

	return 0;
}
