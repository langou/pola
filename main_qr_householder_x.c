#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

extern void qr_householder_x0 ( int M, int N, int P, double A[M][N], double R[N][N], double work[N] );

extern double check_qr_repres(int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern double check_orth(int M, int N, double Q[M][N] );

int main(int argc, char ** argv) {

   int i, j, k, m, n, p, l;
   double norma, norma2, normR, normA;

   m = 10;
   n = 4;

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

   double(*A)[n] = malloc(sizeof(double[m][n]));
   double(*Q)[n] = malloc(sizeof(double[m][n]));
   double(*R)[n] = malloc(sizeof(double[n][n]));
   
//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         Q[i][j] = A[i][j];

   if ( m > n ) p = n; else p = m-1;

/*************************************************************/

   printf("%%%% [ HH_X0         ] m = %4d; n = %4d; ",m,n);
   double(*work) = malloc(sizeof(double[n]));
   qr_householder_x0 (m, n, p, Q, R, work);
   free(work);

/*************************************************************/

   printf("repres = %8.1e; ", check_qr_repres( m, n, A, Q, R ));

   printf("orth = %8.1e;\n", check_orth( m, n, Q ));


//Free memory
   free( R );
   free( Q );
   free( A );

   return 0;

}
