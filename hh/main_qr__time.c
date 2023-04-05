#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "lapacke.h"

#define HH_A2VLL                8
#define HH_A2VLL_BLAS         108
#define HH_A2VLL__TILED       208
#define HH_A2VLL__TILED_BLAS  308
#define UNKNOWN               999

extern void qr_householder_a2vll ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2vll_blas ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2vll__tiled ( int M, int N, int B, double *A, int lda, double *tau );
extern void qr_householder_a2vll__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau );
extern void qr_householder_v2q ( int M, int N, double *A, int lda, double *tau );

extern double check_qr_repres( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_qr_repres_blas( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_orthog( int m, int n, double *Q, int ldq );
extern double check_orthog_blas( int m, int n, double *Q, int ldq );


int main(int argc, char ** argv) {

   int b, i, j, m, n;
   int lda, ldq, ldr;
   int method;

   m = 10;
   n = 8;
   b = 3;
   method = HH_A2VLL;

   for(i = 1; i < argc; i++){
      if( strcmp( *(argv + i), "-m") == 0) {
         m = atoi( *(argv + i + 1) );
         i++;
      }
      if( strcmp( *(argv + i), "-n") == 0) {
         n = atoi( *(argv + i + 1) );
         i++;
      }
      if( strcmp( *(argv + i), "-b") == 0) {
         b = atoi( *(argv + i + 1) );
         i++;
      }
      if( strcmp( *(argv + i), "-method") == 0) {
         if( strcmp( *(argv + i + 1), "hh_a2vll") == 0)
           method = HH_A2VLL;
         else if( strcmp( *(argv + i + 1), "hh_a2vll_blas") == 0)
           method = HH_A2VLL_BLAS;
         else if( strcmp( *(argv + i + 1), "hh_a2vll__tiled") == 0)
           method = HH_A2VLL__TILED;
         else if( strcmp( *(argv + i + 1), "hh_a2vll__tiled_blas") == 0)
           method = HH_A2VLL__TILED_BLAS;
	 else 
           method = UNKNOWN;
         i++;
      }
   }

   if ( m < n ) { printf("m < n is not OK -- m = %d and n= %d -- quick return.\n", m, n ); return -1; }

   if ( method == UNKNOWN ) { printf("method is UNKNOWN -- quick return.\n"); return -1; }

   lda = m; 
   ldq = m;
   ldr = n;

   double *A;
   A = (double *) malloc( lda * n * sizeof(double));
   double *Q;
   Q = (double *) malloc( ldq * n * sizeof(double));
   double *R;
   R = (double *) malloc( ldr * n * sizeof(double));

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i+j*lda] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

/*************************************************************/

   if ( method == HH_A2VLL )             printf("%%%% [ HH_A2VLL              ] m = %4d; n = %4d;           ",m,n);
   if ( method == HH_A2VLL_BLAS )        printf("%%%% [ HH_A2VLL_BLAS         ] m = %4d; n = %4d;           ",m,n);
   if ( method == HH_A2VLL__TILED )      printf("%%%% [ HH_A2VLL__TILED       ] m = %4d; n = %4d; b = %4d; ",m,n,b);
   if ( method == HH_A2VLL__TILED_BLAS ) printf("%%%% [ HH_A2VLL__TILED_BLAS  ] m = %4d; n = %4d; b = %4d; ",m,n,b);

   double(*tau) = malloc(sizeof(double[n]));
   for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];

   if ( method == HH_A2VLL )             qr_householder_a2vll (m, n, Q, ldq, tau);
   if ( method == HH_A2VLL_BLAS )        qr_householder_a2vll_blas (m, n, Q, ldq, tau);
   if ( method == HH_A2VLL__TILED )      qr_householder_a2vll__tiled (m, n, b, Q, ldq, tau);
   if ( method == HH_A2VLL__TILED_BLAS ) qr_householder_a2vll__tiled_blas (m, n, b, Q, ldq, tau);

   for(i = 0; i < n; i++) for(j = i; j < n; j++) R[i+j*ldr] = Q[i+j*ldq];;
   qr_householder_v2q (m, n, Q, ldq, tau);
   //LAPACKE_dorgqr( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau );

   free(tau);

/*************************************************************/

   printf("repres = %8.1e; ", check_qr_repres_blas( m, n, A, lda, Q, ldq, R, ldr ));

   //printf("orth = %8.1e;\n", check_orthog( m, n, Q, ldq ));
   printf("orth = %8.1e;\n", check_orthog_blas( m, n, Q, ldq ));

   free( R );
   free( Q );
   free( A );

   return 0;

}
