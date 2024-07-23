# NumericUtils/Makefile

NAME := NumericUtils
OBJS := bam_$(NAME).o integral.o integral3D.o newton_lnsrch.o nrutil.o \
	f_to_minimize.o fd_jacobian.o ludcmp.o lubksb.o lnsrch_double.o \
	rtsafe_itsP.o attenuation.o Integral_over_sphere.o \
	Integral_over_volume.o Integral_over_volume_4.o \
	path_integral.o \
	fact.o \
	SphericalHarmonics.o SphericalHarmonics_precomp.o\

include $(TOP)/Makefile.subdirs




SphericalHarmonics_precomp.c: SphericalHarmonics_precomp.m
	math < SphericalHarmonics_precomp.m > /dev/null








