CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas
LDFLAGS = -lm

all: main_cholesky.exe main_trtri.exe main_lauum.exe main_cholinv.exe main_cholinv_onesweep.exe main_gghd2.exe main_gehd2_householder.exe main_qr_householder.exe main_qr_gramschmidt.exe main_qr_householder_x.exe

main_gehd2_householder.exe: main_gehd2_householder.c
	$(CC) -o $@ main_gehd2_householder.c $(LDFLAGS)

main_gghd2.exe: main_gghd2.c
	$(CC) -o $@ main_gghd2.c  $(LDFLAGS)

main_cholesky.exe: main_cholesky.c
	$(CC) -o $@ main_cholesky.c  $(LDFLAGS)

main_cholinv.exe: main_cholinv.c
	$(CC) -o $@ main_cholinv.c  $(LDFLAGS)

main_cholinv_onesweep.exe: main_cholinv_onesweep.c
	$(CC) -o $@ main_cholinv_onesweep.c  $(LDFLAGS)

main_trtri.exe: main_trtri.c
	$(CC) -o $@ main_trtri.c  $(LDFLAGS)

main_lauum.exe: main_lauum.c
	$(CC) -o $@ main_lauum.c  $(LDFLAGS)

qr_householder_v0.o: qr_householder_v0.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v0.c

qr_householder_v1.o: qr_householder_v1.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v1.c

main_qr_householder.exe: main_qr_householder.c qr_householder_v0.o qr_householder_v1.o
	$(CC) -o $@ main_qr_householder.c qr_householder_v0.o qr_householder_v1.o $(LDFLAGS)

qr_householder_x0.o: qr_householder_x0.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_x0.c

main_qr_householder_x.exe: main_qr_householder_x.c qr_householder_x0.o
	$(CC) -o $@ main_qr_householder_x.c qr_householder_x0.o $(LDFLAGS)

qr_mgs_v0.o: qr_mgs_v0.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_v0.c

qr_mgs_v1_rl.o: qr_mgs_v1_rl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_v1_rl.c

qr_mgs_v1_ll.o: qr_mgs_v1_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_v1_ll.c

qr_cgs_v1_ll.o: qr_cgs_v1_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs_v1_ll.c

qr_cgs_v1_rl.o: qr_cgs_v1_rl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs_v1_rl.c

qr_cgs2_v1.o: qr_cgs2_v1.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs2_v1.c

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -o $@ check_qr_repres.c

check_orth.o: check_orth.c
	$(CC) -c $(CCFLAGS) -o $@ check_orth.c

main_qr_gramschmidt.exe: main_qr_gramschmidt.c qr_mgs_v0.o qr_mgs_v1_rl.o qr_mgs_v1_ll.o qr_cgs_v1_ll.o qr_cgs_v1_rl.o qr_cgs2_v1.o check_qr_repres.o check_orth.o
	$(CC) -o $@ main_qr_gramschmidt.c qr_mgs_v0.o qr_mgs_v1_rl.o qr_mgs_v1_ll.o qr_cgs_v1_ll.o qr_cgs_v1_rl.o qr_cgs2_v1.o check_qr_repres.o check_orth.o $(LDFLAGS)

clean:
	rm -f *.exe *.o

