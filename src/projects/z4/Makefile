# Z4c/Makefile

NAME := z4
OBJS := bam_$(NAME).o $(NAME).o \
	z4_boundary.o z4_boundary_box.o z4_boundary_box_rad.o z4_boundary_shell.o \
	z4_init.o \
	z4_algcon.o \
	z4_rhs_movpunc_N.o \
	z4_rhs_add_source.o 

# z4_boundary.o \

include $(TOP)/Makefile.subdirs

z4_algcon.c: z4_algcon.m
	math < z4_algcon.m > /dev/null

z4_boundary.c: z4_boundary.m
	math < z4_boundary.m > /dev/null

z4_boundary_box.c: z4_boundary_box.m
	math < z4_boundary_box.m > /dev/null

z4_boundary_shell.c: z4_boundary_shell.m
	math < z4_boundary_shell.m > /dev/null

z4_boundary_box_rad.c: z4_boundary_box_rad.m
	math < z4_boundary_box_rad.m > /dev/null

z4_init.c: z4_init.m
	math < z4_init.m > /dev/null

z4_rhs_add_source.c: z4_rhs_add_source.m
	math < z4_rhs_add_source.m > /dev/null

z4_rhs_movpunc_N.c: z4_rhs_movpunc_N.m
	math < z4_rhs_movpunc_N.m > /dev/null
