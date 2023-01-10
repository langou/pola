#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define CGS_V1    1
#define MGS_V0    2
#define MGS_V1    3
#define CGS2_V1   4
#define UNKNOWN 999

extern void qr_mgs_v0 (int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern void qr_mgs_v1 (int M, int N, double Q[M][N], double R[N][N] );
extern void qr_cgs_v1 (int M, int N, double Q[M][N], double R[N][N] );
extern void qr_cgs2_v1 (int M, int N, double Q[M][N], double R[N][N], double *tmp );

extern double check_qr_repres(int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern double check_orth(int M, int N, double Q[M][N] );


int main(int argc, char ** argv) {

   int i, j, k, m, n, p, l;
   double norma, norma2, normR, normA, tmp;
   int method;

   m = 10;
   n = 4;
   method = MGS_V1;

   for(i = 1; i < argc; i++){
      if( strcmp( *(argv + i), "-m") == 0) {
         m = atoi( *(argv + i + 1) );
         i++;
      }
      if( strcmp( *(argv + i), "-n") == 0) {
         n = atoi( *(argv + i + 1) );
         i++;
      }
      if( strcmp( *(argv + i), "-method") == 0) {
         if( strcmp( *(argv + i + 1), "cgs_v1") == 0)
           method = CGS_V1;
	 else if( strcmp( *(argv + i + 1), "mgs_v0") == 0)
           method = MGS_V0;
	 else if( strcmp( *(argv + i + 1), "mgs_v1") == 0)
           method = MGS_V1;
	 else if( strcmp( *(argv + i + 1), "cgs2_v1") == 0)
           method = CGS2_V1;
	 else 
           method = UNKNOWN;
         i++;
      }
   }

   if ( m < n ) { printf("m < n is not OK. m = %d and n= %d. quick return.\n", m, n ); return -1; }

   if ( method == UNKNOWN ) { printf("method is UNKNOWN. quick return.\n"); return -1; }

   double(*A)[n] = malloc(sizeof(double[m][n]));
   double(*Q)[n] = malloc(sizeof(double[m][n]));
   double(*R)[n] = malloc(sizeof(double[n][n]));

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

/*************************************************************/

   if ( method == MGS_V0 ) { 
      printf("%%%% [ MGS_V0        ] m = %4d; n = %4d; ",m,n);
      double(*B)[n] = malloc(sizeof(double[m][n]));
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) B[i][j] = A[i][j];
      qr_mgs_v0 (m, n, B, Q, R);
      free(B);
   }

   if ( method == MGS_V1 ) { 
      printf("%%%% [ MGS_V1        ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_mgs_v1 (m, n, Q, R);
   }

   if ( method == CGS_V1 ) { 
      printf("%%%% [ CGS_V1        ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_cgs_v1 (m, n, Q, R);
   }

   if ( method == CGS2_V1 ) { 
      printf("%%%% [ CGS2_V1       ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      double(*work) = malloc(sizeof(double[n]));
      qr_cgs2_v1 (m, n, Q, R, work);
      free(work);
   }

/*************************************************************/

   printf("repres = %8.1e; ", check_qr_repres( m, n, A, Q, R ));

   printf("orth = %8.1e;\n", check_orth( m, n, Q ));

   free( R );
   free( Q );
   free( A );

   return 0;

}
