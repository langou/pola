#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <papi/timing.h>

#define MGS_LL                4
#define MGS_LL_BLAS         104
#define MGS_LL__TILED       204
#define MGS_LL__TILED_BLAS  304
#define MGS_RL                5
#define MGS_RL_BLAS         105
#define MGS_RL__TILED       205
#define MGS_RL__TILED_BLAS  305
#define MGS_REC_BLAS        106
#define UNKNOWN             999

extern void qr_mgs_ll             ( int m, int n, double *A, int lda, double *R, int ldr );
extern void qr_mgs_rl             ( int m, int n, double *A, int lda, double *R, int ldr );
extern void qr_mgs_ll_blas        ( int m, int n, double *A, int lda, double *R, int ldr );
extern void qr_mgs_rl_blas        ( int m, int n, double *A, int lda, double *R, int ldr );
extern void qr_mgs_ll__tiled      ( int m, int n, int b, double *A, int lda, double *R, int ldr );
extern void qr_mgs_rl__tiled      ( int m, int n, int b, double *A, int lda, double *R, int ldr );
extern void qr_mgs_ll__tiled_blas ( int m, int n, int b, double *A, int lda, double *R, int ldr );
extern void qr_mgs_rl__tiled_blas ( int m, int n, int b, double *A, int lda, double *R, int ldr );
extern void qr_mgs_rec_blas       ( int m, int n, double *A, int lda, double *R, int ldr );


extern double check_qr_repres( int m, int n, double *A, int lda, double *Q, int ldq, double *R, int ldr );
extern double check_orthog( int m, int n, double *Q, int ldq );


int main(int argc, char ** argv) {

   int b, i, j, m, n;
   int human_readable;
   int use_papi;
   int lda, ldq, ldr;
   int method;

   m = 10;
   n = 8;
   b = 3;
   human_readable = 0;
   use_papi = 1;
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
      if( strcmp( *(argv + i), "-h") == 0) {
	human_readable = 1;
      }
      if( strcmp( *(argv + i), "-np") == 0) {
	use_papi = 0;
      }
      if( strcmp( *(argv + i), "-method") == 0) {
         if( strcmp( *(argv + i + 1), "mgs_ll") == 0)
           method = MGS_LL;
         else if( strcmp( *(argv + i + 1), "mgs_rl") == 0)
           method = MGS_RL;
         else if( strcmp( *(argv + i + 1), "mgs_ll_blas") == 0)
           method = MGS_LL_BLAS;
         else if( strcmp( *(argv + i + 1), "mgs_rl_blas") == 0)
           method = MGS_RL_BLAS;
         else if( strcmp( *(argv + i + 1), "mgs_ll__tiled") == 0)
           method = MGS_LL__TILED;
         else if( strcmp( *(argv + i + 1), "mgs_rl__tiled") == 0)
           method = MGS_RL__TILED;
         else if( strcmp( *(argv + i + 1), "mgs_ll__tiled_blas") == 0)
           method = MGS_LL__TILED_BLAS;
         else if( strcmp( *(argv + i + 1), "mgs_rl__tiled_blas") == 0)
           method = MGS_RL__TILED_BLAS;
         else if( strcmp( *(argv + i + 1), "mgs_rec_blas") == 0)
           method = MGS_REC_BLAS;
	 else 
           method = UNKNOWN;
         i++;
      }
   }

   if ( m < n ) { printf("m < n is not OK. m = %d and n= %d. quick return.\n", m, n ); return -1; }

   if ( method == UNKNOWN ) { printf("method is UNKNOWN. quick return.\n"); return -1; }

   // Initializing PAPI and preparing the data structure for the result
   papi_info_t papi_info;
   
   // To store the results of the counter
   size_t num_events = 0;

   if(use_papi) {
     init_papi();
     papi_info = build_papi_info();
     num_events = papi_info.num_events;
   }

   long long results [num_events];

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

   for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i+j*ldq] = A[i+j*lda];

   if ( method == MGS_LL ) { 
     if(human_readable)
       printf("%%%% [ MGS_LL               ] m = %4d; n = %4d;          ",m,n);
     else
       printf("MGS_LL %4d %4d N/A",m,n);
   }

   if ( method == MGS_RL ) { 
     if(human_readable)
       printf("%%%% [ MGS_RL               ] m = %4d; n = %4d;          ",m,n);
     else
       printf("MGS_RL %4d %4d N/A",m,n);
   }

   if ( method == MGS_LL_BLAS ) { 
     if(human_readable)
       printf("%%%% [ MGS_LL_BLAS          ] m = %4d; n = %4d;          ",m,n);
     else
       printf("MGS_LL_BLAS %4d %4d N/A",m,n);
   }

   if ( method == MGS_RL_BLAS ) { 
     if(human_readable)
       printf("%%%% [ MGS_RL_BLAS          ] m = %4d; n = %4d;          ",m,n);
     else
       printf("MGS_RL_BLAS %4d %4d N/A",m,n);
   }

   if ( method == MGS_LL__TILED ) { 
     if(human_readable)
       printf("%%%% [ MGS_LL__TILED        ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     else
       printf("MGS_LL_TILED %4d %4d %4d",m,n,b);
   }

   if ( method == MGS_RL__TILED ) { 
     if(human_readable)
       printf("%%%% [ MGS_RL__TILED        ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     else
       printf("MGS_RL_TILED %4d %4d %4d",m,n,b);
   }

   if ( method == MGS_LL__TILED_BLAS ) { 
     if(human_readable)
       printf("%%%% [ MGS_LL_TILED_BLAS   ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     else
       printf("MGS_LL_TILED_BLAS %4d %4d %4d",m,n,b);
   }

   if ( method == MGS_RL__TILED_BLAS ) { 
     if(human_readable)
       printf("%%%% [ MGS_RL__TILED_BLAS   ] m = %4d; n = %4d; b = %4d; ",m,n,b);
     else
       printf("MGS_RL_TILED_BLAS %4d %4d %4d",m,n,b);
   }

   if ( method == MGS_REC_BLAS ) { 
     if(human_readable)
       printf("%%%% [ MGS_REC_BLAS         ] m = %4d; n = %4d;          ",m,n);
     else
       printf("MGS_REC_BLAS %4d %4d N/A",m,n);
   }

/*************************************************************/

   struct timespec start, end;
   clock_gettime(CLOCK_MONOTONIC, &start);
   if(use_papi) {
     record_events(papi_info);
   }

   if ( method == MGS_LL ) { 
      qr_mgs_ll (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_RL ) { 
     qr_mgs_rl (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL_BLAS ) { 
     qr_mgs_ll_blas (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_RL_BLAS ) { 
     qr_mgs_rl_blas (m, n, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL__TILED ) { 
     qr_mgs_ll__tiled (m, n, b, Q, ldq, R, ldr);
   }

   if ( method == MGS_RL__TILED ) { 
     qr_mgs_rl__tiled (m, n, b, Q, ldq, R, ldr);
   }

   if ( method == MGS_LL__TILED_BLAS ) { 
     qr_mgs_ll__tiled_blas (m, n, b, Q, ldq, R, ldr);
   }

   if ( method == MGS_RL__TILED_BLAS ) { 
     qr_mgs_rl__tiled_blas (m, n, b, Q, ldq, R, ldr);
   }

   if ( method == MGS_REC_BLAS ) { 
      qr_mgs_rec_blas (m, n, Q, ldq, R, ldr);
   }

   if(use_papi)
     retrieve_results(papi_info, results);
   clock_gettime(CLOCK_MONOTONIC, &end);

/*************************************************************/

   double cycles, l1miss, l2miss = 0;
   if (use_papi) {
     cycles = (double) results[CYCLES];
     l1miss = results[CACHE_MISS_L1];
     l2miss = results[CACHE_MISS_L2];
   }

   double time_taken;
   time_taken = end.tv_sec - start.tv_sec;
   time_taken = time_taken + (end.tv_nsec - start.tv_nsec) * 1e-9;
   
   if(human_readable) {
     printf("repres = %8.1e; ", check_qr_repres( m, n, A, lda, Q, ldq, R, ldr ));
     printf("orth = %8.1e; ", check_orthog( m, n, Q, ldq ));

     printf("time = %.5g; ", time_taken);
     if (use_papi) 
       printf("cycles = %.5g; L1 miss = %.5g; L2 miss = %.5g ", cycles, l1miss, l2miss);
   }
   else {
     printf(" %.5g ", time_taken);
     if (use_papi) 
       printf(" %.5g %.5g %.5g ", cycles, l1miss, l2miss);
     else 
       printf(" N/A N/A N/A ");
   }


   printf("\n");
   free( R );
   free( Q );
   free( A );

   return 0;

}
