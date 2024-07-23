# bssn/Makefile

NAME := bssn
OBJS := bam_$(NAME).o $(NAME).o \
	bssn_init.o \
	bssn_eta.o \
	bssn_rhs_movpunc_N.o \
	bssn_algcon.o \
	bssn_rhs_add_source.o \
	bssn_boundary.o 


include $(TOP)/Makefile.subdirs

bssn_init.c: bssn_init.m
	math < bssn_init.m > /dev/null

bssn_rhs_movpunc_N.c: bssn_rhs_movpunc_N.m
	math < bssn_rhs_movpunc_N.m > /dev/null

bssn_algcon.c: bssn_algcon.m
	math < bssn_algcon.m > /dev/null

bssn_rhs_add_source.c: bssn_rhs_add_source.m
	math < bssn_rhs_add_source.m > /dev/null

bssn_boundary.c: bssn_boundary.m 
	math < bssn_boundary.m > /dev/null


