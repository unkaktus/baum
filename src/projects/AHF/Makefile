# AHF/Makefile

NAME := AHF

OBJS := bam_$(NAME).o AHF.o AHF_output.o \
	AHF_1d_1.o AHF_1d_2.o

include $(TOP)/Makefile.subdirs




AHF_1d_1.c: AHF_1d_1.m
	math < AHF_1d_1.m > /dev/null

AHF_1d_2.c: AHF_1d_2.m
	math < AHF_1d_2.m > /dev/null