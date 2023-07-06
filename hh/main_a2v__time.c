#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "lapacke.h"

#define HH_A2V_LL____________   8
#define HH_A2V_RL____________   9
#define HH_A2V_LL________BLAS 108
#define HH_A2V_RL________BLAS 109
#define HH_A2V_REC_______BLAS 110
#define HH_A2V_LL__TILED_____ 208
#define HH_A2V_RL__TILED_____ 209
#define HH_A2V_LL__TILED_BLAS 308
#define HH_A2V_RL__TILED_BLAS 309
#define GEQR2                 901
#define GEQRF                 902
#define UNKNOWN               999

extern void qr_householder_a2v_rl____________ ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2v_rec_______blas ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2v_rl________blas ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2v_ll____________ ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2v_ll________blas ( int M, int N, double *A, int lda, double *tau );
extern void qr_householder_a2v_ll__tiled_____ ( int M, int N, int B, double *A, int lda, double *tau );
extern void qr_householder_a2v_ll__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau );
extern void qr_householder_a2v_rl__tiled_____ ( int M, int N, int B, double *A, int lda, double *tau );
extern void qr_householder_a2v_rl__tiled_blas ( int M, int N, int B, double *A, int lda, double *tau );

extern double check_qr_repres( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_qr_repres_blas( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_orthog( int m, int n, double *Q, int ldq );
extern double check_orthog_blas( int m, int n, double *Q, int ldq );

int main(int argc, char ** argv) {

   int b, i, j, lwork, m, n;
   int human_readable;
   int lda, ldq, ldr;
   int method;
   double *A, *Q, *R, *tau, *work;

   m = 10;
   n = 8;
   b = 3;
   human_readable = 0;
   method = HH_A2V_LL____________ ;

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
      if( strcmp( *(argv + i), "-h") == 0) {
	human_readable = 1;
      }
     if( strcmp( *(argv + i), "-method") == 0) {
         if( strcmp( *(argv + i + 1), "hh_a2v_rl") == 0)
           method = HH_A2V_RL____________;
         else if( strcmp( *(argv + i + 1), "hh_a2v_rec_blas") == 0)
           method = HH_A2V_REC_______BLAS;
         else if( strcmp( *(argv + i + 1), "hh_a2v_rl_blas") == 0)
           method = HH_A2V_RL________BLAS;
         else if( strcmp( *(argv + i + 1), "hh_a2v_ll") == 0)
           method = HH_A2V_LL____________;
         else if( strcmp( *(argv + i + 1), "hh_a2v_ll_blas") == 0)
           method = HH_A2V_LL________BLAS;
         else if( strcmp( *(argv + i + 1), "hh_a2v_ll_tiled") == 0)
           method = HH_A2V_LL__TILED_____;
         else if( strcmp( *(argv + i + 1), "hh_a2v_ll_tiled_blas") == 0)
           method = HH_A2V_LL__TILED_BLAS;
         else if( strcmp( *(argv + i + 1), "hh_a2v_rl_tiled") == 0)
           method = HH_A2V_RL__TILED_____;
         else if( strcmp( *(argv + i + 1), "hh_a2v_rl_tiled_blas") == 0)
           method = HH_A2V_RL__TILED_BLAS;
         else if( strcmp( *(argv + i + 1), "geqr2") == 0)
           method = GEQR2;
         else if( strcmp( *(argv + i + 1), "geqrf") == 0)
           method = GEQRF;
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

   A = (double *) malloc( lda * n * sizeof(double));
   Q = (double *) malloc( ldq * n * sizeof(double));
   R = (double *) malloc( ldr * n * sizeof(double));

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i+j*lda] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

   tau = (double *) malloc(sizeof(double[n]));

   if ( method == GEQRF ){
      lwork = -1;
      work = (double *) malloc( 1 * sizeof(double));
      LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, lwork );
      lwork = ( (int) work[0] );
      free( work );
      work = (double *) malloc( lwork * sizeof(double));
   }

   if ( method == GEQR2 ){
      work = (double *) malloc( n * sizeof(double));
   }

   if (human_readable) {
     if ( method == HH_A2V_REC_______BLAS ) printf("%%%% [ HH_A2V_REC_______BLAS ] m = %4d; n = %4d;           ",m,n);
     if ( method == HH_A2V_RL____________ ) printf("%%%% [ HH_A2V_RL____________ ] m = %4d; n = %4d;           ",m,n);
     if ( method == HH_A2V_RL________BLAS ) printf("%%%% [ HH_A2V_RL________BLAS ] m = %4d; n = %4d;           ",m,n);
     if ( method == HH_A2V_LL____________ ) printf("%%%% [ HH_A2V_LL____________ ] m = %4d; n = %4d;           ",m,n);
     if ( method == HH_A2V_LL________BLAS ) printf("%%%% [ HH_A2V_LL________BLAS ] m = %4d; n = %4d;           ",m,n);
     if ( method == HH_A2V_LL__TILED_____ ) printf("%%%% [ HH_A2V_LL__TILED_____ ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     if ( method == HH_A2V_LL__TILED_BLAS ) printf("%%%% [ HH_A2V_LL__TILED_BLAS ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     if ( method == HH_A2V_RL__TILED_____ ) printf("%%%% [ HH_A2V_RL__TILED_____ ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     if ( method == HH_A2V_RL__TILED_BLAS ) printf("%%%% [ HH_A2V_RL__TILED_BLAS ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     if ( method == GEQR2         )         printf("%%%% [ GEQR2                 ] m = %4d; n = %4d;           ",m,n);
     if ( method == GEQRF         )         printf("%%%% [ GEQRF                 ] m = %4d; n = %4d;           ",m,n);
   } else {
     if ( method == HH_A2V_REC_______BLAS ) printf("HH_A2V_REC_BLAS       %4d %4d N/A ",m,n);
     if ( method == HH_A2V_RL____________ ) printf("HH_A2V_RL             %4d %4d N/A ",m,n);
     if ( method == HH_A2V_RL________BLAS ) printf("HH_A2V_RL_BLAS        %4d %4d N/A ",m,n);
     if ( method == HH_A2V_LL____________ ) printf("HH_A2V_LL             %4d %4d N/A ",m,n);
     if ( method == HH_A2V_LL________BLAS ) printf("HH_A2V_LL_BLAS        %4d %4d N/A ",m,n);
     if ( method == HH_A2V_LL__TILED_____ ) printf("HH_A2V_LL_TILED       %4d %4d %4d ",m,n,b);
     if ( method == HH_A2V_LL__TILED_BLAS ) printf("HH_A2V_LL_TILED_BLAS  %4d %4d %4d ",m,n,b);
     if ( method == HH_A2V_RL__TILED_____ ) printf("HH_A2V_RL_TILED       %4d %4d %4d ",m,n,b);
     if ( method == HH_A2V_RL__TILED_BLAS ) printf("HH_A2V_RL_TILED_BLAS  %4d %4d %4d ",m,n,b);
     if ( method == GEQR2         )         printf("GEQR2                 %4d %4d N/A ",m,n);
     if ( method == GEQRF         )         printf("GEQRF                 %4d %4d N/A ",m,n);
   }

   for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];

/*************************************************************/

   if ( method == HH_A2V_REC_______BLAS ) qr_householder_a2v_rec_______blas (m, n, Q, ldq, tau);
   if ( method == HH_A2V_RL____________ ) qr_householder_a2v_rl____________ (m, n, Q, ldq, tau);
   if ( method == HH_A2V_RL________BLAS ) qr_householder_a2v_rl________blas (m, n, Q, ldq, tau);
   if ( method == HH_A2V_LL____________ ) qr_householder_a2v_ll____________ (m, n, Q, ldq, tau);
   if ( method == HH_A2V_LL________BLAS ) qr_householder_a2v_ll________blas (m, n, Q, ldq, tau);
   if ( method == HH_A2V_LL__TILED_____ ) qr_householder_a2v_ll__tiled_____ (m, n, b, Q, ldq, tau);
   if ( method == HH_A2V_LL__TILED_BLAS ) qr_householder_a2v_ll__tiled_blas (m, n, b, Q, ldq, tau);
   if ( method == HH_A2V_RL__TILED_____ ) qr_householder_a2v_rl__tiled_____ (m, n, b, Q, ldq, tau);
   if ( method == HH_A2V_RL__TILED_BLAS ) qr_householder_a2v_rl__tiled_blas (m, n, b, Q, ldq, tau);
   if ( method == GEQR2 )                 LAPACKE_dgeqr2_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work );
   if ( method == GEQRF )                 LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, Q, ldq, tau, work, lwork );

/*************************************************************/

   for(i = 0; i < n; i++) for(j = i; j < n; j++) R[i+j*ldr] = Q[i+j*ldq];

   LAPACKE_dorgqr( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau );

   if ( method == GEQR2 ) free(work);
   if ( method == GEQRF ) free(work);

   free(tau);

   if (human_readable) {
     printf("repres = %8.1e; ", check_qr_repres_blas( m, n, A, lda, Q, ldq, R, ldr ));

     printf("orth = %8.1e;", check_orthog_blas( m, n, Q, ldq ));
   }

   printf("\n");

   free( R );
   free( Q );
   free( A );

   return 0;

}
