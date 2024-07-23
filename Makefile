# Makefile
# Bernd Bruegmann, 12/99, 5/02
# Builds the bam executable based on the file MyConfig

# where am I
UNAME := $(shell uname)
TOP   := $(shell pwd)
# to shorten the output: about to be removed for non-standard directories
# TOP   := ../../..

# name the fruit of our labor
EXEC = bam
EXECDIR = $(TOP)/exe
PROJECTDIR = $(TOP)/src/projects

# variables common to all setups
CC = # uses e.g. gcc or icc, set in MyConfig
DFLAGS =
OFLAGS =
WARN = -Wno-unused-result # -Wall

INCS = -I$(TOP)/src/main/main
LIBS = -L$(TOP)/lib
SPECIALLIBS =
libsys = -lm


# --------------------------------------------------------------------------
# some libraries are currently required
libpaths  = src/main/amr 
libpaths += src/utility/AwA
libpaths += src/utility/boundary 
libpaths += src/utility/checkpoint
libpaths += src/utility/elliptic/iterative
libpaths += src/utility/elliptic/multigrid
libpaths += src/utility/evolve
libpaths += src/utility/interpolate
libpaths += src/utility/output 
libpaths += src/utility/NumericUtils
libpaths += src/utility/units

# --------------------------------------------------------------------------
# the user choses the libraries and some options in the file MyConfig

projects =#
include MyConfig
libpaths += $(projects)
projectnamesnocore = $(notdir $(projects))
projectnames = "../.." $(projectnamesnocore)


# --------------------------------------------------------------------------
# manage how the bam sources are compiled

# you probably want to run with MPI (actually required for now)
libpaths += src/main/bampi

# note that the order matters, e.g.
# main has to go last since it has to be compiled last
libpaths += src/main/main

# extract the list of directory names
libdirs = $(dir $(libpaths))

# extract list of names
libnames = $(notdir $(libpaths))

# make the list of libraries
liblist := $(foreach libname,$(libnames),-l$(libname))

# remove -lmain from that list
liblist := $(subst -lmain,,$(liblist))

# make final list of libraries that is passed to the linker
# uses standard hack to resolve interdependencies by repeating the libraries
# system libraries go in the end
LIBS += $(MPIDIRL) $(liblist) $(liblist)
LIBS += $(SPECIALLIBS) $(libsys) $(MPILIBS)

# make the list of include files that will be automatically included for each
# module
libincludes := $(foreach libpath,$(libpaths),\
	$(libpath)/bam_$(notdir $(libpath)))

# define the automatic configuration files
autoinclude = src/main/main/bam_automatic_include.h
autoinitial = src/main/main/bam_automatic_initialize.c
autofinal   = src/main/main/bam_automatic_finalize.c
autotext    = \/\* automatically generated from MyConfig \*\/

# --------------------------------------------------------------------------
# some of the above variables are meant to be global, so pass them on
# to the shell 
CFLAGS = $(DFLAGS) $(OFLAGS) $(INCS) $(MPIDIRI) $(WARN)
MFILES = $(shell find $(TOP)/src/math/MathToC -name "*.m") 
export

# Automatic includes compile before anything else, thus not parallel.
.bamauto: $(autoinclude) $(autofinal) $(autoinitial)
.NOTPARALLEL: .bamauto

# --------------------------------------------------------------------------
# --------------------------------------------------------------------------
# default target
bam: .bamauto
	mkdir -p $(EXECDIR);
	mkdir -p $(PROJECTDIR); 
	for X in $(libnames); do mkdir -p lib/obj/$$X; done
	for X in $(libpaths); do $(MAKE) -C $$X; done

# --------------------------------------------------------------------------
# other targets 

# automatic configuration files
$(autoinclude): MyConfig
	echo $(autotext) > $(autoinclude) 
	for X in $(libincludes); do \
	  echo \#include \"$(TOP)/$$X.h\" >> $(autoinclude); \
	done

$(autoinitial): MyConfig
	echo $(autotext) > $(autoinitial) 
	for X in $(libnames); do \
	  echo bam\_$$X\(\)\; >> $(autoinitial); \
	  echo void bam\_$$X\(\)\; >> $(autoinclude); \
	done

$(autofinal): MyConfig
	echo $(autotext) > $(autofinal) 
	for X in $(libincludes); do \
	if grep '_\final' $$X.c > /dev/null ; then \
		echo $$(basename bam\_$$X\_final\(\)\;) >> $(autofinal); \
		echo void $$(basename bam\_$$X\_final\(\)\;) >> $(autoinclude); \
	fi \
	done

# take a fresh look at things
clean:
	-rm -r lib 
	-rm $(autoinclude)
	-rm $(autoinitial)
	-rm $(autofinal)

# make clean and recompile
new: clean bam