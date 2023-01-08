CC = gcc
CCFLAGS = 
LDFLAGS = -lm

all: main_cholesky.exe main_trtri.exe main_lauum.exe main_cholinv.exe main_cholinv_onesweep.exe main_gghd2.exe main_gehd2_householder.exe main_qr_householder.exe

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

main_qr_householder.exe: main_qr_householder.c qr_householder_v0.o
	$(CC) -o $@ main_qr_householder.c qr_householder_v0.o $(LDFLAGS)

clean:
	rm -f *.exe *.o

