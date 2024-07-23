# boundary/Makefile

NAME := boundary

OBJS := bam_$(NAME).o boundary.o excision.o extrapolate.o inversion.o periodic.o plane.o radiative.o robin.o symmetry.o tube.o normals.o

deps := $(shell find $(TOP)/src/math -name "*.m") $(shell find $(TOP)/src/utility -name "*.h") 


include $(TOP)/Makefile.subdirs
