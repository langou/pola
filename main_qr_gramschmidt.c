#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

extern void qr_mgs_v0 (int M, int N, double **A, double **Q, double **R );
extern void qr_mgs_v1 (int M, int N, double **Q, double **R );
extern void qr_cgs_v1 (int M, int N, double **Q, double **R );

int main(int argc, char ** argv) {

   int i, j, k, m, n, p, l;
   double *AAA, *BBB, *QQQ, *RRR;
   double **A, **Q, **B, **R;
   double norma, norma2, normR, normA, tmp;

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

   AAA = (double *) malloc( m * n * sizeof(double) );
   A = (double **) malloc( m * sizeof(double*) );
   for(i = 0; i < m; i++) A[i] = AAA + i * n ;

   QQQ = (double *) malloc( m * n * sizeof(double) );
   Q = (double **) malloc( m * sizeof(double*) );
   for(i = 0; i < m; i++) Q[i] = QQQ + i * n ;

   BBB = (double *) malloc( m * n * sizeof(double) );
   B = (double **) malloc( m * sizeof(double*) );
   for(i = 0; i < m; i++) B[i] = BBB + i * n ;

   RRR = (double *) malloc( n * n * sizeof(double) );
   R = (double **) malloc( n * sizeof(double*) );
   for(i = 0; i < n; i++) R[i] = RRR + i * n ;

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;


//   Save a copy of A in B (so that we can later check)
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         B[i][j] = A[i][j];

   if ( m < n ) printf("m < n is not OK. m = %d and n= %d\n", m, n );

/*************************************************************/

   //qr_mgs_v0 (m, n, A, Q, R);


   //for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
   //qr_mgs_v1 (m, n, Q, R);

   for(i = 0; i < m; i++) for(j = 0; j < n; j++) Q[i][j] = A[i][j];
   qr_cgs_v1 (m, n, Q, R);


/*************************************************************/

//printf("A=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",A[i][j]);} printf("\n");} printf("];\n");

//printf("B=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",B[i][j]);} printf("\n");} printf("];\n");

//printf("Q=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",Q[i][j]);} printf("\n");} printf("];\n");

   printf("[ MGS   ] m = %4d; n = %4d;",m,n);

   normR = 0.0e+00;

   for(i = 0; i < m; i++){

      for(j = 0; j < n; j++){

         tmp = B[i][j];

         for(k = 0; k <= j; k++){

            tmp -= Q[i][k] * R[k][j];

         }

         normR += tmp * tmp;

      }

   }

   normR = SQRT_FUN( normR );

   normA = 0.0e+00;

   for(i = 0; i < m; i++){

      for(j = 0; j < n; j++){

         normA += B[i][j] * B[i][j];

      }

   }

   normA = SQRT_FUN( normA );
 
   printf(" %e;", normR/normA );

   normR = 0.0e+00;

   for(i = 0; i < n; i++){

      for(j = 0; j < n; j++){

         if ( i==j ) tmp = 1.0e+00; else tmp = 0.0e+00; 

         for(k = 0; k < m; k++){

            tmp -= Q[k][i] * Q[k][j];

         }

         normR += tmp * tmp;

      }

   }

   normR = SQRT_FUN( normR );
 
   printf(" %e;", normR );

   printf("\n" );

//Free memory
   free( R );
   free( RRR );
   free( Q );
   free( QQQ );
   free( B );
   free( BBB );
   free( A );
   free( AAA );

   return 0;

}
