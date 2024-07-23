# test/Makefile
# build three executables for testing bam

DIR = ../../par/test
DIR = .

all: $(DIR)/bamtest $(DIR)/bamcompare $(DIR)/floatdiff $(DIR)/mpitest

$(DIR)/bamtest: bamtest.c
	$(CC) -O2 -o $(DIR)/bamtest bamtest.c $(libsys)

$(DIR)/bamcompare: bamcompare.c
	$(CC) -O2 -o $(DIR)/bamcompare bamcompare.c $(libsys)


$(DIR)/floatdiff: floatdiff.c wolfio.c wolfio.h
	$(CC) -o $(DIR)/floatdiff  floatdiff.c wolfio.c $(libsys)

$(DIR)/mpitest: mpitest.c
	$(CC) -O2 $(MPIDIRI) -o $(DIR)/mpitest mpitest.c $(libsys) $(MPIDIRL) $(MPILIBS)