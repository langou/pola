#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

extern void qr_householder_a2vll__tiled ( int M, int N, int B, double A[M][N], double tau[N] );
extern void qr_householder_v2q ( int M, int N, double A[M][N], double tau[N] );

extern double check_qr_repres( int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern double check_orth( int M, int N, double Q[M][N] );


int main(int argc, char ** argv) {

   int i, j, m, n, b;

   m = 10;
   n = 4;
   b = 2;

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
   }

   if ( m < n ) { printf("m < n is not OK. m = %d and n= %d. quick return.\n", m, n ); return -1; }

   double(*A)[n] = malloc(sizeof(double[m][n]));
   double(*Q)[n] = malloc(sizeof(double[m][n]));
   double(*R)[n] = malloc(sizeof(double[n][n]));

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

/*************************************************************/

   printf("%%%% [ HH_A2VLL_V2Q  ] m = %4d; n = %4d; b = %4d; ",m,n,b);
   double(*tau) = malloc(sizeof(double[n]));
   for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
   qr_householder_a2vll__tiled (m, n, b, Q, tau);
   for(i = 0; i < n; i++) for(j = i; j < n; j++) R[i][j] = Q[i][j];
   qr_householder_v2q (m, n, Q, tau);
   free(tau);

/*************************************************************/

   printf("repres = %17.10e; ", check_qr_repres( m, n, A, Q, R ));

   printf("orth = %17.10e;\n", check_orth( m, n, Q ));

   free( R );
   free( Q );
   free( A );

   return 0;

}
