# ADM_mass/Makefile

NAME := ADM_mass
OBJS := bam_$(NAME).o ADM_mass.o \
	adm_mass_integrand.o \
	adm_mass_integrand_chi.o \
	adm_pj_integrand.o \

include $(TOP)/Makefile.subdirs

adm_mass_integrand.c: adm_mass_integrand.m
	math < adm_mass_integrand.m > /dev/null

adm_pj_integrand.c: adm_pj_integrand.m
	math < adm_pj_integrand.m > /dev/null

adm_mass_integrand_chi.c: adm_mass_integrand_chi.m
	math < adm_mass_integrand_chi.m > /dev/null

