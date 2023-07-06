CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas -static-libgcc -lgcc -O3 -march=native
LDFLAGS = -lm -lgfortran
CBLAS_INCLUDE = ${MKLROOT}/include/
CBLAS_LIB = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl
# module tools/lapacke/3.9.1 adds the correct directories to C_INCLUDE_PATH, so no need to provide more.
# LAPACKE_INCLUDE = /Users/langou/Documents/repositories/lapack.git/LAPACKE/include/
# LAPACKE_LIB = /Users/langou/Documents/repositories/lapack.git/liblapacke.a /Users/langou/Documents/repositories/lapack.git/liblapack.a

all: main_a2v__time.exe main_v2q__time.exe

qr_householder_a2v_rec_______blas.o: qr_householder_a2v_rec_______blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2v_rec_______blas.c

qr_householder_a2v_rl____________.o: qr_householder_a2v_rl____________.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2v_rl____________.c

qr_householder_a2v_rl________blas.o: qr_householder_a2v_rl________blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2v_rl________blas.c

qr_householder_a2v_ll____________.o: qr_householder_a2v_ll____________.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2v_ll____________.c

qr_householder_a2v_ll________blas.o: qr_householder_a2v_ll________blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2v_ll________blas.c

qr_householder_a2v_ll__tiled_____.o: qr_householder_a2v_ll__tiled_____.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2v_ll__tiled_____.c

qr_householder_a2v_ll__tiled_blas.o: qr_householder_a2v_ll__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2v_ll__tiled_blas.c

qr_householder_a2v_rl__tiled_____.o: qr_householder_a2v_rl__tiled_____.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2v_rl__tiled_____.c

qr_householder_a2v_rl__tiled_blas.o: qr_householder_a2v_rl__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_a2v_rl__tiled_blas.c

qr_householder_v2q_rl____________.o: qr_householder_v2q_rl____________.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q_rl____________.c

qr_householder_v2q_rl________blas.o: qr_householder_v2q_rl________blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_v2q_rl________blas.c

qr_householder_v2q_rec_______blas.o: qr_householder_v2q_rec_______blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_v2q_rec_______blas.c

qr_householder_v2q_ll____________.o: qr_householder_v2q_ll____________.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q_ll____________.c

qr_householder_v2q_ll________blas.o: qr_householder_v2q_ll________blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_v2q_ll________blas.c

qr_householder_v2q_ll__tiled_____.o: qr_householder_v2q_ll__tiled_____.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q_ll__tiled_____.c

qr_householder_v2q_ll__tiled_blas.o: qr_householder_v2q_ll__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_v2q_ll__tiled_blas.c

qr_householder_v2q_rl__tiled_____.o: qr_householder_v2q_rl__tiled_____.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q_rl__tiled_____.c

qr_householder_v2q_rl__tiled_blas.o: qr_householder_v2q_rl__tiled_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ qr_householder_v2q_rl__tiled_blas.c

check_orthog.o: check_orthog.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_orthog.c

check_orthog_blas.o: check_orthog_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_orthog_blas.c

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_qr_repres.c

check_qr_repres_blas.o: check_qr_repres_blas.c
	$(CC) -c $(CCFLAGS) -I$(CBLAS_INCLUDE) -o $@ check_qr_repres_blas.c

main_a2v__time.o: main_a2v__time.c
	$(CC) -c $(CCFLAGS) -I$(LAPACKE_INCLUDE) -I$(CBLAS_INCLUDE) -o $@ main_a2v__time.c

main_v2q__time.o: main_v2q__time.c
	$(CC) -c $(CCFLAGS) -I$(LAPACKE_INCLUDE) -I$(CBLAS_INCLUDE) -o $@ main_v2q__time.c

main_a2v__time.exe: main_a2v__time.o qr_householder_a2v_rl__tiled_blas.o qr_householder_a2v_rl__tiled_____.o qr_householder_a2v_rec_______blas.o qr_householder_a2v_rl____________.o qr_householder_a2v_rl________blas.o qr_householder_a2v_ll____________.o qr_householder_a2v_ll________blas.o qr_householder_a2v_ll__tiled_blas.o qr_householder_a2v_ll__tiled_____.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o
	gfortran -o $@ main_a2v__time.o qr_householder_a2v_rl__tiled_blas.o qr_householder_a2v_rl__tiled_____.o qr_householder_a2v_rec_______blas.o qr_householder_a2v_rl____________.o qr_householder_a2v_rl________blas.o qr_householder_a2v_ll____________.o qr_householder_a2v_ll________blas.o qr_householder_a2v_ll__tiled_____.o qr_householder_a2v_ll__tiled_blas.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o $(LAPACKE_LIB) $(CBLAS_LIB) $(LDFLAGS) 

main_v2q__time.exe: main_v2q__time.o qr_householder_v2q_rl__tiled_blas.o qr_householder_v2q_rl__tiled_____.o qr_householder_v2q_rl____________.o qr_householder_v2q_rl________blas.o qr_householder_v2q_rec_______blas.o qr_householder_v2q_ll____________.o qr_householder_v2q_ll________blas.o qr_householder_v2q_ll__tiled_____.o qr_householder_v2q_ll__tiled_blas.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o
	gfortran -o $@ main_v2q__time.o qr_householder_v2q_rl__tiled_blas.o qr_householder_v2q_rl__tiled_____.o qr_householder_v2q_rl____________.o qr_householder_v2q_rl________blas.o qr_householder_v2q_rec_______blas.o qr_householder_v2q_ll____________.o qr_householder_v2q_ll________blas.o qr_householder_v2q_ll__tiled_____.o qr_householder_v2q_ll__tiled_blas.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o $(LAPACKE_LIB) $(CBLAS_LIB) $(LDFLAGS) 
# main_qr__time.exe: main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled_blas.o qr_householder_a2vll__tiled.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o
# 	gcc -o $@ main_qr__time.o qr_householder_a2vll.o qr_householder_a2vll_blas.o qr_householder_v2q.o qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled_blas.o check_qr_repres.o check_qr_repres_blas.o check_orthog.o check_orthog_blas.o $(LAPACKE_LIB) $(CBLAS_LIB) $(LDFLAGS) 

clean:
	rm -f *.exe *.o

