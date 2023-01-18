CC = gcc
CCFLAGS = -Wall -Wno-unknown-pragmas
LDFLAGS = -lm

all: main_cholesky.exe main_trtri.exe main_lauum.exe main_cholinv.exe main_cholinv_onesweep.exe main_gghd2.exe main_gehd2_householder.exe main_qr.exe

main_gehd2_householder.exe: main_gehd2_householder.o
	$(CC) -o $@ main_gehd2_householder.o $(LDFLAGS)

main_gghd2.exe: main_gghd2.o
	$(CC) -o $@ main_gghd2.o  $(LDFLAGS)

main_cholesky.exe: main_cholesky.o
	$(CC) -o $@ main_cholesky.o  $(LDFLAGS)

main_cholinv.exe: main_cholinv.o
	$(CC) -o $@ main_cholinv.o  $(LDFLAGS)

main_cholinv_onesweep.exe: main_cholinv_onesweep.o
	$(CC) -o $@ main_cholinv_onesweep.o  $(LDFLAGS)

main_trtri.exe: main_trtri.o
	$(CC) -o $@ main_trtri.o  $(LDFLAGS)

main_lauum.exe: main_lauum.o
	$(CC) -o $@ main_lauum.o  $(LDFLAGS)

main_gehd2_householder.o: main_gehd2_householder.c
	$(CC) -c $(CCFLAGS) -o $@ main_gehd2_householder.c

main_gghd2.o: main_gghd2.c
	$(CC) -c $(CCFLAGS) -o $@ main_gghd2.c

main_cholesky.o: main_cholesky.c
	$(CC) -c $(CCFLAGS) -o $@ main_cholesky.c

main_cholinv.o: main_cholinv.c
	$(CC) -c $(CCFLAGS) -o $@ main_cholinv.c

main_cholinv_onesweep.o: main_cholinv_onesweep.c
	$(CC) -c $(CCFLAGS) -o $@ main_cholinv_onesweep.c

main_trtri.o: main_trtri.c
	$(CC) -c $(CCFLAGS) -o $@ main_trtri.c

main_lauum.o: main_lauum.c
	$(CC) -c $(CCFLAGS) -o $@ main_lauum.c

qr_householder_a2vrl.o: qr_householder_a2vrl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2vrl.c

qr_householder_a2v_ll.o: qr_householder_a2vll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2vll.c

qr_householder_v2q.o: qr_householder_v2q.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_v2q.c

qr_householder_a2q.o: qr_householder_a2q.c
	$(CC) -c $(CCFLAGS) -o $@ qr_householder_a2q.c

qr_mgs_pb.o: qr_mgs_pb.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_pb.c

qr_mgs_rl.o: qr_mgs_rl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_rl.c

qr_mgs_ll.o: qr_mgs_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_mgs_ll.c

qr_cgs_ll.o: qr_cgs_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs_ll.c

qr_cgs_rl.o: qr_cgs_rl.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs_rl.c

qr_cgs2_ll.o: qr_cgs2_ll.c
	$(CC) -c $(CCFLAGS) -o $@ qr_cgs2_ll.c

check_qr_repres.o: check_qr_repres.c
	$(CC) -c $(CCFLAGS) -o $@ check_qr_repres.c

check_orth.o: check_orth.c
	$(CC) -c $(CCFLAGS) -o $@ check_orth.c

main_qr.o: main_qr.c
	$(CC) -c $(CCFLAGS) -o $@ main_qr.c

main_qr.exe: main_qr.o qr_householder_a2vll.o qr_householder_a2vrl.o qr_householder_v2q.o qr_householder_a2q.o qr_mgs_pb.o qr_mgs_rl.o qr_mgs_ll.o qr_cgs_ll.o qr_cgs_rl.o qr_cgs2_ll.o check_qr_repres.o check_orth.o
	$(CC) -o $@ main_qr.o qr_householder_a2vll.o qr_householder_a2vrl.o qr_householder_v2q.o qr_householder_a2q.o qr_mgs_pb.o qr_mgs_rl.o qr_mgs_ll.o qr_cgs_ll.o qr_cgs_rl.o qr_cgs2_ll.o check_qr_repres.o check_orth.o $(LDFLAGS)

clean:
	rm -f *.exe *.o

