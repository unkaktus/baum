# Invariants/Makefile

NAME := Invariants
OBJS := bam_$(NAME).o $(NAME).o \
	curvature_invariants_N.o \
	dphi_int2_psi4_N.o 
	

include $(TOP)/Makefile.subdirs

curvature_invariants_N.c: curvature_invariants_N.m
	math < curvature_invariants_N.m > /dev/null
    
dphi_int2_psi4_N.c: dphi_int2_psi4_N.m
	math < dphi_int2_psi4_N.m > /dev/null
