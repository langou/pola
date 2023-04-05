CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas
LDFLAGS = -lm
CBLAS_INCLUDE = /Users/langou/Documents/repositories/lapack.git/CBLAS/include/
CBLAS_LIB = /Users/langou/Documents/repositories/lapack.git/libcblas.a /Users/langou/Documents/repositories/lapack.git/librefblas.a
LAPACKE_INCLUDE = /Users/langou/Documents/repositories/lapack.git/LAPACKE/include/
LAPACKE_LIB = /Users/langou/Documents/repositories/lapack.git/liblapacke.a /Users/langou/Documents/repositories/lapack.git/liblapack.a


all: main_qr__time.exe

qr_householder_a2vll.o: qr_householder_a2vll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2vll.c

qr_householder_a2vll_blas.o: qr_householder_a2vll_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2vll_blas.c

qr_householder_v2q.o: qr_householder_v2q.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q.c

qr_householder_a2vll__tiled.o: qr_householder_a2vll__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2vll__tiled.c

qr_householder_a2vll__tiled_blas.o: qr_householder_a2vll__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2vll__tiled_blas.c

main_qr__time.o: main_qr__time.c
	$(CC) -c $(CCFLAGS) -I$(LAPACKE_INCLUDE) -I$(CBLAS_INCLUDE) -o $@ main_qr__time.c

main_qr__time.exe: main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled_blas.o qr_householder_a2vll__tiled.o check_qr_repres.o check_orthog.o
	gcc -o $@ main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled_blas.o check_qr_repres.o check_orthog.o $(CBLAS_LIB) $(LDFLAGS)
	#gfortran -o $@ main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled_blas.o check_qr_repres.o check_orthog.o $(LAPACKE_LIB) $(CBLAS_LIB) $(LDFLAGS)

clean:
	rm -f *.exe *.o

