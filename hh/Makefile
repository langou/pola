CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas -static-libgcc -lgcc 
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

check_orthog.o: check_orthog.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_orthog.c

check_orthog_blas.o: check_orthog_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_orthog_blas.c

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_qr_repres.c

check_qr_repres_blas.o: check_qr_repres_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_qr_repres_blas.c

main_qr__time.o: main_qr__time.c
	$(CC) -c $(CCFLAGS) -I$(LAPACKE_INCLUDE) -I$(CBLAS_INCLUDE) -o $@ main_qr__time.c

main_qr__time.exe: main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled_blas.o qr_householder_a2vll__tiled.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o
	gcc -o $@ main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled_blas.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o $(LAPACKE_LIB) $(CBLAS_LIB) /usr/local/Cellar/gcc/12.2.0/lib/gcc/12/libgfortran.a $(LDFLAGS) 

clean:
	rm -f *.exe *.o

