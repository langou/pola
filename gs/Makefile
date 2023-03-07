CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas
LDFLAGS = -lm
CBLAS_INCLUDE = /Users/langou/Documents/repositories/lapack.git/CBLAS/include/
CBLAS_LIB = /Users/langou/Documents/repositories/lapack.git/libcblas.a /Users/langou/Documents/repositories/lapack.git/librefblas.a


all: main_qr__time.exe

qr_mgs_ll.o: qr_mgs_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_ll.c

qr_mgs_rl.o: qr_mgs_rl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_rl.c

qr_mgs_ll_blas.o: qr_mgs_ll_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_ll_blas.c 

qr_mgs_rl_blas.o: qr_mgs_rl_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_rl_blas.c 

qr_mgs_ll__tiled.o: qr_mgs_ll__tiled.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_ll__tiled.c 

qr_mgs_rl__tiled.o: qr_mgs_rl__tiled.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_rl__tiled.c 

qr_mgs_ll__tiled_blas.o: qr_mgs_ll__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_ll__tiled_blas.c 

qr_mgs_rl__tiled_blas.o: qr_mgs_rl__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_mgs_rl__tiled_blas.c 

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -o $@ check_qr_repres.c

check_orthog.o: check_orthog.c
	$(CC) -c $(CCFLAGS) -o $@ check_orthog.c

main_qr__time.o: main_qr__time.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr__time.c

main_qr__time.exe: main_qr__time.o qr_mgs_ll.o qr_mgs_rl.o qr_mgs_ll_blas.o qr_mgs_rl_blas.o qr_mgs_ll__tiled.o qr_mgs_rl__tiled.o qr_mgs_ll__tiled_blas.o qr_mgs_rl__tiled_blas.o check_qr_repres.o check_orthog.o
	$(CC) -o $@ main_qr__time.o qr_mgs_ll.o qr_mgs_rl.o qr_mgs_ll_blas.o qr_mgs_rl_blas.o qr_mgs_ll__tiled.o qr_mgs_rl__tiled.o qr_mgs_ll__tiled_blas.o qr_mgs_rl__tiled_blas.o check_qr_repres.o check_orthog.o $(CBLAS_LIB) $(LDFLAGS) 

clean:
	rm -f *.exe *.o

