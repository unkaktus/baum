# Gauge/Makefile

NAME := Gauge
OBJS := bam_$(NAME).o \
	moving_punctures.o moving_punctures_integrate.o \
        puncture_properties.o  puncture_properties_integrand.o \
	
#	maximal.o maximal_L.o maximal_init.o \
#	commotion.o corotation.o scale_shift.o \
#	moving_punctures.o moving_punctures_integrate.o \
#	alpha_center.o \
#	generic_Gauge.o

include $(TOP)/Makefile.subdirs

maximal_L.c: maximal_L.m
	math < maximal_L.m > /dev/null

maximal_init.c: maximal_init.m
	math < maximal_init.m > /dev/null

puncture_properties_integrand.c: puncture_properties_integrand.m
	math < puncture_properties_integrand.m > /dev/null
