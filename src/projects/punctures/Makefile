# punctures/Makefile

NAME := punctures
OBJS := bam_$(NAME).o ConformalFactor.o ExtrinsicCurvature.o FlagRegrid.o \
 	Lapse.o Mass.o Punctures.o PuncturesPS.o PunctureSolve.o Shift.o \
	VectorLieFlat.o VectorLaplaceFlat.o VectorLaplaceFlatGS.o OutsideBox.o \
	FullBox.o FuncsForOutsideBox.o Sequences.o \
	forMovingPunctures.o PunctureData_with_given_M_ADM.o

include $(TOP)/Makefile.subdirs

VectorLieFlat.c: VectorLieFlat.m
	math < VectorLieFlat.m > /dev/null
VectorLaplaceFlat.c: VectorLaplaceFlat.m
	math < VectorLaplaceFlat.m > /dev/null
VectorLaplaceFlatGS.c: VectorLaplaceFlatGS.m
	math < VectorLaplaceFlatGS.m > /dev/null
