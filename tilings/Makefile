CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas
LDFLAGS = -lm

all: main_qr__check_qr_mgs_rl__tiled.exe main_qr__check_qr_mgs_ll__tiled.exe main_qr__check_qr_cgs_rl__tiled.exe main_qr__check_qr_householder_a2vll__tiled.exe

qr_mgs_rl__tiled.o: qr_mgs_rl__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_rl__tiled.c

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -o $@ check_qr_repres.c

check_orth.o: check_orth.c
	$(CC) -c $(CCFLAGS) -o $@ check_orth.c

main_qr__check_qr_mgs_rl__tiled.o: main_qr__check_qr_mgs_rl__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr__check_qr_mgs_rl__tiled.c

main_qr__check_qr_mgs_rl__tiled.exe: main_qr__check_qr_mgs_rl__tiled.o qr_mgs_rl__tiled.o check_qr_repres.o check_orth.o
	$(CC) -o $@ main_qr__check_qr_mgs_rl__tiled.o qr_mgs_rl__tiled.o check_qr_repres.o check_orth.o $(LDFLAGS)

qr_mgs_ll__tiled.o: qr_mgs_ll__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_ll__tiled.c

main_qr__check_qr_mgs_ll__tiled.o: main_qr__check_qr_mgs_ll__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr__check_qr_mgs_ll__tiled.c

main_qr__check_qr_mgs_ll__tiled.exe: main_qr__check_qr_mgs_ll__tiled.o qr_mgs_ll__tiled.o check_qr_repres.o check_orth.o
	$(CC) -o $@ main_qr__check_qr_mgs_ll__tiled.o qr_mgs_ll__tiled.o check_qr_repres.o check_orth.o $(LDFLAGS)

qr_cgs_rl__tiled.o: qr_cgs_rl__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs_rl__tiled.c

main_qr__check_qr_cgs_rl__tiled.o: main_qr__check_qr_cgs_rl__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr__check_qr_cgs_rl__tiled.c

main_qr__check_qr_cgs_rl__tiled.exe: main_qr__check_qr_cgs_rl__tiled.o qr_cgs_rl__tiled.o check_qr_repres.o check_orth.o
	$(CC) -o $@ main_qr__check_qr_cgs_rl__tiled.o qr_cgs_rl__tiled.o check_qr_repres.o check_orth.o $(LDFLAGS)

qr_householder_a2vll__tiled.o: qr_householder_a2vll__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2vll__tiled.c

qr_householder_v2q.o: qr_householder_v2q.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q.c

main_qr__check_qr_householder_a2vll__tiled.o: main_qr__check_qr_householder_a2vll__tiled.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr__check_qr_householder_a2vll__tiled.c

main_qr__check_qr_householder_a2vll__tiled.exe: main_qr__check_qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled.o check_qr_repres.o check_orth.o qr_householder_v2q.o
	$(CC) -o $@ main_qr__check_qr_householder_a2vll__tiled.o qr_householder_a2vll__tiled.o check_qr_repres.o check_orth.o qr_householder_v2q.o $(LDFLAGS)



clean:
	rm -f *.exe *.o

