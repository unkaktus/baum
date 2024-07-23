# matter/Makefile

NAME := grhd
OBJS := bam_$(NAME).o $(NAME).o \
	grhd_atm.o \
	grhd_Vac.o \
	grhd_c2p.o \
	grhd_c2p_v2.o \
	grhd_c2p_Vac.o \
	grhd_c2p_hroot.o \
	grhd_eig.o \
	grhd_flx1d.o \
	grhd_flx1d_Vac.o \
	grhd_flx1d_turb.o \
	grhd_ppm.o \
	grhd_flx1d_HLLC_PPM.o \
	grhd_p2c.o \
	grhd_sources.o \
	grhd_sources_Vac.o \
	grhd_sources_turb.o \
	grhd_excision.o	
	

#	grhd_p2a.o \
#	grhd_a2p.o \

include $(TOP)/Makefile.subdirs

grhd_sources.c: grhd_sources.m
	math < grhd_sources.m > /dev/null

grhd_sources_Vac.c: grhd_sources_Vac.m
	math < grhd_sources_Vac.m > /dev/null
	
grhd_sources_turb.c: grhd_sources_turb.m
	math < grhd_sources_turb.m > /dev/null

#grhd_a2p.c: grhd_a2p.m
#	math < grhd_a2p.m > /dev/null

#grhd_p2a.c: grhd_p2a.m
#	math < grhd_p2a.m > /dev/null
