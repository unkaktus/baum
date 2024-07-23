# Makefile
# Bernd Bruegmann, 12/99, 5/02
# Builds the bam executable based on the file MyConfig
# See http://www.gnu.org/software/make/manual for the manual of GNU make

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
bam: .bamauto bamtest
	mkdir -p $(EXECDIR);
	mkdir -p $(PROJECTDIR); 
	for X in $(libnames); do mkdir -p lib/obj/$$X; done
	for X in $(libpaths); do $(MAKE) -C $$X; done

# --------------------------------------------------------------------------
# other targets 

# if there is no MyConfig file, use the example provided in doc
MyConfig:
	-if test ! -f MyConfig; then cp doc/MyConfig.example MyConfig; fi 


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

# build executables for testing
bamtest:
	$(MAKE) -C src/test

# copy test suites into standard location
# note that we remove all old tests first (which may remove other files, too)
# note that the update option in "cp -u" is therefore redundant 
testnew: testclean
	find ./src -name "*_correct" -exec cp -pru {} ./par/test/ \;
	-rm -rf ./par/test/*_correct/.git
	cd ./par/test; cp -pu *_correct/*.par .

testclean:
	find ./par/test/* -not \( \
	-name .git -o -name Readme.txt \
	-o -name bamtest -o -name bamcompare -o -name floatdiff \) \
	-exec rm -rf {} \;

# run test suites
test:
	cd par/test; ./bamtest -2 -rt 1e-12 -at 8e-15 *.par

test1:
	cd par/test; ./bamtest -1 -rt 1e-12 -at 8e-15 *1proc.par

test2:
	cd par/test; ./bamtest -2 -rt 1e-12 -at 8e-15 *2proc.par

# create tar file
tar:
	cd ..; tar czf bam.tgz --exclude lib --exclude exe --exclude .git ./bam
tarbzip:
	cd ..; tar cjf bam.tbz --exclude lib --exclude exe --exclude .git ./bam
tarbox:
	cd ..; tar cjf bambox.tbz --exclude lib --exclude exe --exclude par --exclude .git ./bambox

# take a fresh look at things
clean: #testclean
	-rm -r lib 
	-rm $(autoinclude)
	-rm $(autoinitial)
	-rm $(autofinal)

# make clean and recompile
new: clean bam

# remove emacs backup files
tildeclean:
	find . -name "*~" -exec rm {} \;

# git 
git_clone:
	mkdir -p $(PROJECTDIR)
	for X in $(projectnames); do mkdir -p $(PROJECTDIR)/$$X; git clone git@git.tpi.uni-jena.de:bamdev/$$X.git $(PROJECTDIR)/$$X; done

git_clone_https:
	mkdir -p $(PROJECTDIR)
	for X in $(projectnames); do mkdir -p $(PROJECTDIR)/$$X; git clone https://git.tpi.uni-jena.de/bamdev/$$X.git $(PROJECTDIR)/$$X; done

git_pull:
	for X in $(projectnames); do if [ -d "$(PROJECTDIR)/$$X" ]; then echo $$X; cd $(PROJECTDIR)/$$X; git pull; fi done

git_clone_bb:
	mkdir -p $(PROJECTDIR)
	for X in $(projectnames); do git clone https://bruegmann@git.tpi.uni-jena.de/bamdev/$$X.git $(PROJECTDIR)/$$X; done

branch_BAM23:
	for X in $(projectnames); do cd $(PROJECTDIR)/$$X; git checkout BAM23; done


# svn for projects --- old
# svn_commit: 	svn_ci
# svn_checkout: 	svn_co
# svn_diff: 	svn_di
# svn_status: 	svn_st
# svn_stat: 	svn_st
# svn_update: 	svn_up

# svn_ci:	
# 	for X in $(projectnames); do svn commit src/projects/$$X; done 
# svn_cl:	
# 	for X in $(projectnames); do svn cleanup src/projects/$$X; done 
# svn_co:
# 	for X in $(projectnamesnocore); do svn checkout https://svn.tpi.uni-jena.de/numrel/group/bam/14.07/projects/$$X src/projects/$$X; done 
# svn_co_munich:
# 	for X in $(projectnamesnocore); do svn checkout https://localhost:10011/numrel/group/bam/12.04/projects/$$X src/projects/$$X; done 
# svn_di:
# 	for X in $(projectnames); do svn diff -rHEAD src/projects/$$X; done 
# svn_st:
# 	for X in $(projectnames); do svn status src/projects/$$X; done 
# svn_sts:
# 	@for X in $(projectnames); do svn status src/projects/$$X | grep -v ?; done
# svn_up:
# 	for X in $(projectnames); do svn update src/projects/$$X; done






