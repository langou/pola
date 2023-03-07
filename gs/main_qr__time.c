#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MGS_LL                4
#define MGS_LL_BLAS         104
#define MGS_LL__TILED       204
#define MGS_LL__TILED_BLAS  304
#define MGS_RL                5
#define UNKNOWN             999

extern void qr_mgs_ll             ( int m, int n, double *Q, int ldq, double *R, int ldr );
extern void qr_mgs_ll_blas        ( int m, int n, double *Q, int ldq, double *R, int ldr );
extern void qr_mgs_ll__tiled      (int m, int n, int b, double *A, int lda, double *R, int ldr );
extern void qr_mgs_ll__tiled_blas (int m, int n, int b, double *A, int lda, double *R, int ldr );


extern double check_qr_repres( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_orthog( int m, int n, double *Q, int ldq );


int main(int argc, char ** argv) {

   int b, i, j, m, n;
   int lda, ldq, ldr;
   int method;

   m = 10;
   n = 8;
   b = 3;
   method = MGS_LL;

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
         if( strcmp( *(argv + i + 1), "mgs_ll") == 0)
           method = MGS_LL;
         else if( strcmp( *(argv + i + 1), "mgs_ll_blas") == 0)
           method = MGS_LL_BLAS;
         else if( strcmp( *(argv + i + 1), "mgs_ll__tiled") == 0)
           method = MGS_LL__TILED;
         else if( strcmp( *(argv + i + 1), "mgs_ll__tiled_blas") == 0)
           method = MGS_LL__TILED_BLAS;
	 else 
           method = UNKNOWN;
         i++;
      }
   }

   if ( m < n ) { printf("m < n is not OK. m = %d and n= %d. quick return.\n", m, n ); return -1; }

   if ( method == UNKNOWN ) { printf("method is UNKNOWN. quick return.\n"); return -1; }

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

   if ( method == MGS_LL ) { 
      printf("%%%% [ MGS_LL               ] m = %4d; n = %4d;          ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];
      qr_mgs_ll (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL_BLAS ) { 
      printf("%%%% [ MGS_LL_BLAS          ] m = %4d; n = %4d;          ",m,n);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];
      qr_mgs_ll_blas (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL__TILED ) { 
      printf("%%%% [ MGS_LL__TILED        ] m = %4d; n = %4d; b = %4d; ",m,n,b);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];
      qr_mgs_ll__tiled (m, n, b, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL__TILED_BLAS ) { 
      printf("%%%% [ MGS_LL__TILED_BLAS   ] m = %4d; n = %4d; b = %4d; ",m,n,b);
      for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];
      qr_mgs_ll__tiled_blas (m, n, b, Q, ldq, R, ldr);
   }

/*************************************************************/

   printf("repres = %8.1e; ", check_qr_repres( m, n, A, lda, Q, ldq, R, ldr ));

   printf("orth = %8.1e;\n", check_orthog( m, n, Q, ldq ));

   free( R );
   free( Q );
   free( A );

   return 0;

}
