#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define CGS_LL         1
#define CGS_RL         2
#define MGS_PB         3
#define MGS_LL         4
#define MGS_RL         5
#define CGS2_LL        6
#define CGS2_RL        7
#define HH_A2VLL_V2Q   8
#define HH_A2VRL_V2Q   9
#define HH_A2Q        10
#define UNKNOWN      999

extern void qr_mgs_pb (int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern void qr_mgs_ll (int M, int N, double Q[M][N], double R[N][N] );
extern void qr_mgs_rl (int M, int N, double Q[M][N], double R[N][N] );
extern void qr_cgs_ll (int M, int N, double Q[M][N], double R[N][N] );
extern void qr_cgs_rl (int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern void qr_cgs2_ll (int M, int N, double Q[M][N], double R[N][N], double *tmp );
extern void qr_householder_a2vll ( int M, int N, double A[M][N], double tau[N] );
extern void qr_householder_a2vrl ( int M, int N, double A[M][N], double tau[N] );
extern void qr_householder_v2q ( int M, int N, double A[M][N], double tau[N] );
extern void qr_householder_a2q ( int M, int N, double A[M][N], double R[N][N], double work[N] );

extern double check_qr_repres(int M, int N, double A[M][N], double Q[M][N], double R[N][N] );
extern double check_orth(int M, int N, double Q[M][N] );


int main(int argc, char ** argv) {

   int i, j, m, n;
   int method;

   m = 10;
   n = 4;
   method = MGS_RL;

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
         if( strcmp( *(argv + i + 1), "cgs_rl") == 0)
           method = CGS_RL;
	 else if( strcmp( *(argv + i + 1), "cgs_ll") == 0)
           method = CGS_LL;
	 else if( strcmp( *(argv + i + 1), "mgs_pb") == 0)
           method = MGS_PB;
	 else if( strcmp( *(argv + i + 1), "mgs_rl") == 0)
           method = MGS_RL;
	 else if( strcmp( *(argv + i + 1), "mgs_ll") == 0)
           method = MGS_LL;
	 else if( strcmp( *(argv + i + 1), "cgs2_ll") == 0)
           method = CGS2_LL;
	 else if( strcmp( *(argv + i + 1), "hh_a2vll_v2q") == 0)
           method = HH_A2VLL_V2Q;
	 else if( strcmp( *(argv + i + 1), "hh_a2vrl_v2q") == 0)
           method = HH_A2VRL_V2Q;
	 else if( strcmp( *(argv + i + 1), "hh_a2q") == 0)
           method = HH_A2Q;
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

   if ( method == MGS_PB ) { 
      printf("%%%% [ MGS_PB        ] m = %4d; n = %4d; ",m,n);
      double(*B)[n] = malloc(sizeof(double[m][n]));
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) B[i][j] = A[i][j];
      qr_mgs_pb (m, n, B, Q, R);
      free(B);
   }

   if ( method == MGS_RL ) { 
      printf("%%%% [ MGS_RL        ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_mgs_rl (m, n, Q, R);
   }

   if ( method == MGS_LL ) { 
      printf("%%%% [ MGS_LL        ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_mgs_ll (m, n, Q, R);
   }

   if ( method == CGS_RL ) { 
      printf("%%%% [ CGS_RL        ] m = %4d; n = %4d; ",m,n);
      qr_cgs_rl (m, n, A, Q, R);
   }
   
   if ( method == CGS_LL ) { 
      printf("%%%% [ CGS_LL        ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_cgs_ll (m, n, Q, R);
   }

   if ( method == CGS2_LL ) { 
      printf("%%%% [ CGS2_LL       ] m = %4d; n = %4d; ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      double(*work) = malloc(sizeof(double[n]));
      qr_cgs2_ll (m, n, Q, R, work);
      free(work);
   }

   if ( method == HH_A2VLL_V2Q ) { 
      printf("%%%% [ HH_A2VLL_V2Q  ] m = %4d; n = %4d; ",m,n);
      double(*tau) = malloc(sizeof(double[n]));
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_householder_a2vll (m, n, Q, tau);
      for(i = 0; i < n; i++) for(j = i; j < n; j++) R[i][j] = Q[i][j];
      qr_householder_v2q (m, n, Q, tau);
      free(tau);
   }

   if ( method == HH_A2VRL_V2Q ) { 
      printf("%%%% [ HH_A2VRL_V2Q  ] m = %4d; n = %4d; ",m,n);
      double(*tau) = malloc(sizeof(double[n]));
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_householder_a2vrl (m, n, Q, tau);
      for(i = 0; i < n; i++) for(j = i; j < n; j++) R[i][j] = Q[i][j];
      qr_householder_v2q (m, n, Q, tau);
      free(tau);
   }

   if ( method == HH_A2Q ) { 
      printf("%%%% [ HH_A2Q        ] m = %4d; n = %4d; ",m,n);
      double(*work) = malloc(sizeof(double[n]));
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
      qr_householder_a2q (m, n, Q, R, work);
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
