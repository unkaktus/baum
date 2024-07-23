# hydroanalysis/Makefile

NAME := hydroanalysis

OBJS := bam_$(NAME).o hydroanalysis.o compute_hydrovars.o compute_magvars.o

include $(TOP)/Makefile.subdirs

compute_hydrovars.c: compute_hydrovars.m
	math < compute_hydrovars.m > /dev/null