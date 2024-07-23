# adm/Makefile

NAME := adm
OBJS := bam_$(NAME).o $(NAME).o \
	adm_constraints_N.o \
	adm_data.o \
	adm_utility.o

include $(TOP)/Makefile.subdirs

adm_constraints_N.c: adm_constraints_N.m
	math < adm_constraints_N.m > /dev/null