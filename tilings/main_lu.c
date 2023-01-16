#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_N n
#define SQRT_FUN sqrt

int main(int argc, char ** argv) {

	int i, j, k, n, nb;
	double **A, **B;
	double normA, normR, tmp;

	srand(0);

    	n = 22;
	nb = 4;

	for(i = 1; i < argc; i++){
		if( strcmp( *(argv + i), "-n") == 0) {
			n  = atoi( *(argv + i + 1) );
			i++;
		}
		if( strcmp( *(argv + i), "-nb") == 0) {
			nb  = atoi( *(argv + i + 1) );
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

//	Create a random matrix A. 
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

//	Save a copy of A in B
 	for(i = 0; i < n; i++)
 		for(j = 0; j < n; j++)
			B[i][j] = A[i][j];

/*************************************************************/

	long int number_read, number_write, cache_useage, max_cache_useage;

	number_read = 0;
	number_write = 0;
	cache_useage = 0;
	max_cache_useage = 0;

	int jb, kb;
	int ii, jj, kk;

 	for( j = 0; j < n; j+=nb ){

		jb = (( nb < n-j ) ? ( nb ) : ( n-j ));

 		for( k= j+jb; k < n; k+=nb ){

			kb = (( nb < n-k ) ? ( nb ) : ( n-k ));

//			read A[ k:k+kb, j:j+jb ]
			number_read += ( kb * jb );
			cache_useage += ( kb * jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 			for( kk= 0; kk < j; kk++ ){

//				read A[ k:k+kb, kk ]
//				read A[ kk, j:j+jb ]
				number_read += ( kb + jb );
				cache_useage += ( kb + jb );
				max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 				for( ii= k; ii < k+kb; ii++ )
 					for( jj= j; jj < j+jb; jj++ )
						A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

//				erase A[ k:k+kb, kk ]
//				erase A[ kk, j:j+jb ]
				cache_useage -= ( kb + jb );
				max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

			}

//			write A[ k:k+kb, j:j+jb ]
			number_write += ( kb * jb );
			cache_useage -= ( kb * jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

		}

 		for( k= j+jb; k < n; k+=nb ){

			kb = (( nb < n-k ) ? ( nb ) : ( n-k ));

//			read A[ j:j+jb, k:k+kb ]
			number_read += ( jb * kb );
			cache_useage += ( jb * kb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 			for( kk= 0; kk < j; kk++ ){

//				read A[ j:j+jb, kk ]
//				read A[ kk, k:k+kb ]
				number_read += ( jb + kb );
				cache_useage += ( jb + kb );
				max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 				for( ii= j; ii < j+jb; ii++ )
 					for( jj= k; jj < k+kb; jj++ )
						A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

//				erase A[ j:j+jb, kk ]
//				erase A[ kk, k:k+kb ]
				cache_useage -= ( jb + kb );
				max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

			}

//			write A[ j:j+jb, k:k+kb ]
			number_write += ( jb * kb );
			cache_useage -= ( jb * kb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

		}

//		read A[ 1:jb, 1:jb ]
		number_read += ( jb * jb );
		cache_useage += ( jb * jb );
		max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 		for( kk= 0; kk < j; kk++ ){

//			read A[ 1:jb, kk ]
//			read A[ kk, 1:jb ]
			number_read += ( jb + jb );
			cache_useage += ( jb + jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 			for( ii= j; ii < j+jb; ii++ )
 				for( jj= j; jj < j+jb; jj++ )
					A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

//			erase A[ 1:jb, kk ]
//			erase A[ kk, 1:jb ]
			cache_useage -= ( jb + jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

		}


 		for( jj = j; jj < j+jb; jj++ ){

			for(ii = j; ii < jj; ii++)
 				for(kk = j; kk < ii; kk++)
					A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

			for(ii = jj; ii < j+jb; ii++)
 				for(kk = j; kk < jj; kk++)
					A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

 			for( ii = jj+1; ii < j+jb; ii++ )
				A[ ii ][ jj ] /= A[ jj ][ jj ];

		}


		for(ii = j+jb; ii < n; ii++){

			// read A[ ii, 1:jb ]
			number_read += ( jb );
			cache_useage += ( jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 			for( jj = j; jj < j+jb; jj++ ){

 				for(kk = j; kk < jj; kk++)
					A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];

				A[ ii ][ jj ] /= A[ jj ][ jj ];
			}

			// write A[ ii, 1:jb ]
			number_write += ( jb );
			cache_useage -= ( jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

		}


		for(jj = j+jb; jj < n; jj++){

			// read A[ 1:jb, jj ]
			number_read += ( jb );
			cache_useage += ( jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

 			for( ii = j; ii < j+jb; ii++ ){

 				for(kk = j; kk < ii; kk++){
					A[ ii ][ jj ] -= A[ ii ][ kk ] * A[ kk ][ jj ];
				}

			}

			// write A[ 1:jb, jj ]
			number_write += ( jb );
			cache_useage -= ( jb );
			max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));

		}

//		write A[ 1:jb, 1:jb ]
		number_write += ( jb * jb );
		cache_useage -= ( jb * jb );
		max_cache_useage = (( max_cache_useage > cache_useage ) ? ( max_cache_useage ) : ( cache_useage ));


	}

//	check
	normR = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = 0; i <= j; i++) {
			tmp = B[i][j];
			for (k = 0; k < i ; k++) {
				tmp -= A[i][k]*A[k][j];
			}
			tmp -= A[i][j];
			normR += tmp * tmp;
		}
		for (i = j+1; i < n; i++) {
			tmp = B[i][j];
			for (k = 0; k <= j ; k++) {
				tmp -= A[i][k]*A[k][j];
			}
			normR += tmp * tmp;
		}
	}
	normR = sqrt( normR );
	normA = 0e+00;
	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {
			tmp = B[i][j];
			normA += tmp * tmp;
		}
	}
	normA = sqrt( normA );
	printf("[ GETRF ] n = %4d; nb = %4d; normR / normA = %6.2e;\n", n, nb, normR/normA);

	long int number_read_, number_write_, max_cache_useage_;

	max_cache_useage_ = nb * nb + 2 * nb ;
	number_read_ = ( ( 2 * n * n * n / nb + 3 * n * n - 2 * nb * n ) / 3 ) ;
	number_write_ = 2 * n * n - n * nb ;

//	printf(" %4d %4d ", n, nb);
//	printf(" %10ld", number_read );
//	printf(" %10ld", number_write );
//	printf(" %10ld", max_cache_useage );
//	printf("\n");

//	printf("number_read      = %10ld;\n", number_read );
//	printf("number_write     = %10ld;\n", number_write );
//	printf("max_cache_useage = %10ld;\n", max_cache_useage );

	printf("number_read      = %10ld %10ld %2ld;\n", number_read, number_read_, labs(number_read - number_read_));
	printf("number_write     = %10ld %10ld %2ld;\n", number_write, number_write_, labs(number_write - number_write_));
	printf("max_cache_useage = %10ld %10ld %2ld;\n", max_cache_useage, max_cache_useage_, labs(max_cache_useage - max_cache_useage_));

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
