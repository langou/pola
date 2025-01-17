#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_M m
#define _PB_N n
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

// We must have m >= n.
//
// This routine takes an m-by-n matrix A and reduces it to bidiagonal form
// using an orthogonal transformation. The output of the subroutine
// is a bidiagonal matrix B and orthogonal matrices P and Q such that 
//      A = P * H * Q^T
//
// This subroutine is the first step to solve the singular value problem. The
// equivalent LAPACK subroutines is DGEBD2. (The next step is for example
// DBDSQR, or DBDSDC, or DBDSVDX.)
//
// If only the singular values are needed, (and not the singular vectors,) then
// only B needs to be given to DBDSQR.  (I.e. no need to compute the orthogonal
// factors P and Q.)
//
// If we want to compute the the eigenvalues and the eigenvectors, then we need
// to give to the QR algorithm the two matrices H and Q .
//
// We use Householder transformation. This is a level-2 BLAS subroutine.
//
// This implementation requires a workspace of size n. (`tmp` in the code
// below.)
//
// We perform three numerical checks: (1) || A * Q' - Q * H ||_f / || A ||_f is
// of the order of machine precision (2) || Q' * Q - I ||_f is of the order of
// machine precision (3) H is upper Hessenberg 
//
// This implementation replaces the input A by the output H and the Householder
// vector matrix V stored beneath H. This is akin to LAPACK DGEHD2.
//
// The matrix Q is not computed by DGEHD2. We do the same and we post-process
// with a sequence of operations similar to DORGHR to explicitly compute Q for
// the check.  (For the check, we could also have applied Q implicitly instead
// of explicitly computing Q.)

extern void dgebd2_( int *m, int *n, double *A, int *lda,
		double *d, double *e, double *tauq, double *taup, double *work, int *info );

extern void dorgbr_( char *vect, int *m, int *n, int *k, double *A, int *lda,
		double *tau, double *work, int *lwork, int *info );

extern void dlacpy_( char *uplo, int *m, int *n, double *a, int *lda, double *b, int *ldb );

int main(int argc, char ** argv) {

	int i, j, k, m, n;
	double *AAA, *BBB, *QQQ;
	double **A, **B, **Q;
	double *tmp;
	double tau, norma, norma2, normA, normR;
	double *d, *e, *taup, *tauq, *work;
	int info;
	double *A_r, *Q_r;

	srand(0);

	m = 30;
	n = 20;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-m") == 0) {
			m  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
	}

	AAA = (double *) malloc( m * n * sizeof(double));
	A = (double **) malloc( m * sizeof(double*) );
 	for(i = 0; i < m; i++){
		A[i] = AAA + i * n ;
	}

	BBB = (double *) malloc( m * n * sizeof(double) );
	B = (double **) malloc( m * sizeof(double*) );
 	for(i = 0; i < m; i++){
		B[i] = BBB + i * n ;
	}

	tmp = (double *) malloc( n * sizeof(double) );

//	Create a random n-by-n matrix A.
 	for(i = 0; i < m; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Save a copy of A in B (so that we can check)
 	for(i = 0; i < m; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

	QQQ = (double *) malloc( m * n * sizeof(double));
	Q = (double **) malloc( m * sizeof(double*) );
 	for(i = 0; i < m; i++){
		Q[i] = QQQ + i * n ;
	}

	printf("A = [");
	for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { printf("%+8.6f, ", A[i][j] ); } printf(";" ); }
	printf("];\n");

	A_r = (double *) malloc( m * n * sizeof(double));
	Q_r = (double *) malloc( m * n * sizeof(double));
	for (i = 0; i < m; i++) for (j = 0; j < n; j++) A_r[ i + j * m ] = A[i][j];


	d = (double *) malloc( n * sizeof(double));
	e = (double *) malloc( (n-1) * sizeof(double));
	taup = (double *) malloc( n * sizeof(double));
	tauq = (double *) malloc( n * sizeof(double));
	work = (double *) malloc( m * sizeof(double));

	dgebd2_( &m, &n, A_r, &m, d, e, tauq, taup, work, &info );
	free( work );

	printf("%% [ dgebd2     ] info = %d\n", info );
	printf("Af = [");
	for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { printf("%+8.6f, ", A_r[i+j*m] ); } printf(";" ); }
	printf("];\n");
	printf("Bf = triu(tril(Af,1))\n");

	dlacpy_( "A", &m, &n, A_r, &m, Q_r, &m );

	work = (double *) malloc( n * sizeof(double));
	dorgbr_( "Q", &m, &n, &n, Q_r, &m, tauq, work, &n, &info );
	free( work );

	for (i = 0; i < m; i++) for (j = 0; j < n; j++) Q[i][j] = Q_r[ i + j * m ];

	printf("d=[" ); for (i = 0; i < n; i++) { printf("%+8.6f ", d[i]); } printf("];\n" ); 
	printf("e=[" ); for (i = 0; i < n-1; i++) { printf("%+8.6f ", e[i]); } printf("];\n" ); 
	printf("%% [ dorgbr (1) ] info = %d\n", info );
	for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { printf("%+8.6f ", Q[i][j] ); } printf("\n" ); }

	free( A_r );
	free( Q_r );

 	for(i = 0; i < m; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];


#pragma scop

int M = m;
int N = n;
double ttmp;

for(k = 0; k < N; k++){
// on the left
   norma2 = 0.e+00;
   for(i = k+1; i < M; i++){
      norma2 += B[i][k] * B[i][k];
   }
   norma = sqrt( B[k][k] * B[k][k] + norma2 );
   
   B[k][k] = ( B[k][k] > 0 ) ? ( B[k][k] + norma ) : ( B[k][k] - norma ) ;
   
   tauq[k] = 2.0 / ( 1.0 + norma2 / ( B[k][k] * B[k][k] ) ) ;
   
   for(i = k+1; i < M; i++){
      B[i][k] /= B[k][k];
   }
   B[k][k]= ( B[k][k] > 0 ) ? ( - norma ) : ( norma ) ;
   
   for(j = k+1; j < N; j++){
      ttmp = B[k][j];
      for(i = k+1; i < M; i++){
         ttmp += B[i][k] * B[i][j];
      }
      ttmp = tauq[k] * ttmp;
      B[k][j] = B[k][j] - ttmp;
      for(i = k+1; i < M; i++){
         B[i][j] = B[i][j] - B[i][k] * ttmp;
      }
   }

// on the right
   norma2 = 0.e+00;
   for(j = k+2; j < N; j++){
      norma2 += B[k][j] * B[k][j];
   }
   norma = sqrt( B[k][k+1] * B[k][k+1] + norma2 );
   
   B[k][k+1] = ( B[k][k+1] > 0 ) ? ( B[k][k+1] + norma ) : ( B[k][k+1] - norma ) ;
   
   taup[k] = 2.0 / ( 1.0 + norma2 / ( B[k][k+1] * B[k][k+1] ) ) ;
   
   for(j = k+2; j < N; j++){
      B[k][j] /= B[k][k+1];
   }
   B[k][k+1]= ( B[k][k+1] > 0 ) ? ( - norma ) : ( norma ) ;
   
   for(i = k+1; i < M; i++){
      ttmp = B[i][k+1];
      for(j = k+2; j < N; j++){
         ttmp += B[i][j] * B[k][j];
      }
      ttmp = ttmp * taup[k];
      B[i][k+1] = B[i][k+1] - ttmp;
      for(j = k+2; j < N; j++){
         B[i][j] = B[i][j] - ttmp * B[k][j] ;
      }
   }

}

	printf("B = [");
	for (i = 0; i < m; i++) { for (j = 0; j < n; j++) { printf("%+8.6f, ", B[i][j] ); } printf(";" ); }
	printf("];"); 
	//printf("]; for j=3:n,B(1,j)=0;end;for i=2:m, B(i,1)=0; end;\n");
	printf("Bf = triu(tril(B,1))\n");


#pragma endscop









	return 0;
	printf("%% [ GEBD2 ] m = %4d; n = %4d; checks = [ ", m, n );

//	check || Q' * Q - I ||_f
	normR = 0e+00;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			tmp[0] = ( i == j ) ? 1.0e+00 : 0.00e+00;
			for (k = 0; k < m; k++) {
				tmp[0] -= Q[k][i]*Q[k][j];
			}
			normR += tmp[0] * tmp[0];
		}
	}
	normR = sqrt( normR );
	printf(" %6.2e", normR );


	free( Q );
	free( QQQ );

	free( tauq );
	free( taup );
	free( e );
	free( d );


/*************************************************************/
//#pragma scop
// for ( j = 0; j < _PB_N-2; j++ ) {
//    norma2 = SCALAR_VAL(0.0) ;
//    for ( i = j+1; i < _PB_M; i++ ) {
//       norma2 += A[i][j] * A[i][j] ;
//    }
//    norma = SQRT_FUN ( A[j][j] * A[j][j] + norma2 ) ;
//    A[j][j] = ( A[j][j] > 0 ) ? ( A[j][j] + norma ) : ( A[j][j] - norma ) ;
//    tau = SCALAR_VAL(2.0) / ( SCALAR_VAL(1.0) + norma2 / ( A[j][j] * A[j][j] ) ) ;
//    for ( i = j+1; i < _PB_M; i++ ) {
//       A[i][j] /= A[j][j] ;
//    }
//    A[j][j] = ( A[j][j] > 0 ) ? ( - norma ) : ( + norma ) ;
//    for (i = j+1; i < _PB_N; i++) {
//       tmp[i] = A[j][i] ;
//       for (k = j+1; k < _PB_M; k++) {
//          tmp[i] += A[k][j] * A[k][i];
//       }
//    }
//    for (i = j+1; i < _PB_N; i++) {
//       tmp[i] *= tau ;
//    }
//    for (i = j+1 ; i < _PB_N ; i++) {
//       A[j][i] -= tmp[i] ;
//    }
//    for (i = j+1; i < _PB_M; i++) {
//       for (k = j+1; k < _PB_N; k++) {
//          A[i][k] -= A[i][j] * tmp[k];
//       }
//    }
//     for (i = 0; i < _PB_N; i++) {
//        tmp[i] = A[i][j+1];
//        for (k = j+2; k < _PB_N; k++) {
//           tmp[i] += A[i][k] * A[k][j];
//        }
//     }
//     for (i = 0; i < _PB_N; i++) {
//        tmp[i] *= tau ;
//     }
//     for (i = 0; i < _PB_N; i++) {
//        A[i][j+1] -= tmp[i] ;
//     }
//     for (i = 0; i < _PB_N; i++) {
//        for (k = j+2; k < _PB_N; k++) {
//           A[i][k] -= tmp[i] * A[k][j] ;
//        }
//     }
// }
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
	printf("[ GEHD2 ] m = %4d; n = %4d; checks = [ ", m, n );
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
//	this needs to be a loop no?
	free( QQQ );

	free( tmp );

	free( B );
	free( BBB );

	free( A );
	free( AAA );

	return 0;
}
