#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#define _PB_P p
#define _PB_M m
#define _PB_N n
#define SCALAR_VAL(x) x
#define SQRT_FUN sqrt

extern void qr_householder_v0 ( int M, int N, int P, double **A, double *tmp );


int main(int argc, char ** argv) {

   int i, j, k, m, n, p, l;
   double *AAA, *BBB, *QQQ;
   double **A, **Q, **B;
   double *tmp;
   double norma, norma2, tau, normR, normA;

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
   
   BBB = (double *) malloc( m * n * sizeof(double) );
   B = (double **) malloc( m * sizeof(double*) );
   for(i = 0; i < m; i++) B[i] = BBB + i * n ;
   
   //tmp = (double *) malloc( 1 * sizeof(double) );
   //tmp = (double *) malloc( n * sizeof(double) );

//   Create a random m-by-n matrix A.
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         A[i][j] = (double)rand() / (double)(RAND_MAX) - 0.5e+00;


//   Save a copy of A in B (so that we can later check)
   for(i = 0; i < m; i++)
      for(j = 0; j < n; j++)
         B[i][j] = A[i][j];

   if ( m > n ) p = n; else p = m-1;

/*************************************************************/

   qr_householder_v0 (m, n, p, A, tmp);

/*************************************************************/
//free(tmp);

//printf("A=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",A[i][j]);} printf("\n");} printf("];\n");

//printf("B=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",B[i][j]);} printf("\n");} printf("];\n");

   tmp = (double *) malloc( m * sizeof(double) );

   printf("[ GEQR2 ] m=%4d; n = %4d;",m,n);

   if ( m > n ) p = n; else p = m-1;
   QQQ = (double *) malloc( m * n * sizeof(double));
   Q = (double **) malloc( m * sizeof(double*) );
   for(i = 0; i < m; i++) Q[i] = QQQ + i * n ;

   for(k = p-1; k > -1; k--){

      norma2 = 0.e+00;

      for(i = k+1; i < m; i++){

         norma2 += A[i][k] * A[i][k];

      }

      tau = 2.0e+00 / ( 1.0e+00 + norma2 );

      for(j = k+1; j < n; j++){

         tmp[j] = 0.e+00;

            for(i = k+1; i < m; i++){

               tmp[j] += A[i][k] * Q[i][j];

            }

      }

      tmp[k] = 1.0e+00;

      for(j = k; j < n; j++){ 

         tmp[j] *= tau;

      }

      Q[k][k] = 1.0e+00 - tmp[k];

      for(j = k+1; j < n; j++){

         Q[k][j] = -tmp[j];

      }

      for(j = k; j < n; j++){

         for(i = k+1; i < m; i++){

            Q[i][j] -= A[i][k] * tmp[j];

         }
      }
   }

//printf("Q=[\n");
//for(i = 0; i < m; i++){
//for(j = 0; j < n; j++){
//printf("%+e ",Q[i][j]);} printf("\n");} printf("];\n");


   normR = 0.0e+00;

   for(i = 0; i < m; i++){

      for(j = 0; j < n; j++){

         tmp[0] = B[i][j];

         for(k = 0; k <= j; k++){

            tmp[0] -= Q[i][k] * A[k][j];

         }

         normR += tmp[0] * tmp[0];

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

         if ( i==j ) tmp[0] = 1.0e+00; else tmp[0] = 0.0e+00; 

         for(k = 0; k < m; k++){

            tmp[0] -= Q[k][i] * Q[k][j];

         }

         normR += tmp[0] * tmp[0];

      }

   }

   normR = SQRT_FUN( normR );
 
   printf(" %e;", normR );

   printf("\n" );



//Free memory
   free( Q );
   free( QQQ );
   free( tmp );
   free( B );
   free( BBB );
   free( A );
   free( AAA );

   return 0;

}
