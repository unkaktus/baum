# amr/Makefile

NAME := amr

OBJS := bam_amr.o box.o check.o flag.o grid.o inject.o move.o print.o points.o \
	regrid.o shells.o storage.o


include $(TOP)/Makefile.subdirs

bam_amr.o: bam_amr_loops_box.h bam_amr_advection.h
