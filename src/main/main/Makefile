# main/Makefile

NAME := main

OBJS := bam_main.o main.o parameters.o physics.o skeleton.o system.o \
	tensors.o utilities.o variables.o


# main is special, so we do not use Makefile.subdirs here

# directory of the libraries
LIBDIR  = $(TOP)/lib

# directory for the object files of this library
OBJDIR  = $(LIBDIR)/obj/$(NAME)
OBJSWITHPATH = $(addprefix $(OBJDIR)/,$(OBJS))


# main target is the executable
# note that $@ is the matched target
$(EXECDIR)/$(EXEC): $(OBJSWITHPATH) $(wildcard $(LIBDIR)/*)
	$(CC) $(CFLAGS) -o $@ $(OBJSWITHPATH) $(LIBS)

# secondary target: build object files
# note that $@ is the matched target, 
# $* is the stem of the target defined by %
$(OBJDIR)/%.o: %.c bam.h bam_$(NAME).h $(NAME).h Makefile
	$(CC) $(CFLAGS) -o $@ -c $*.c -DREVISION="$$(git rev-parse  HEAD)"


# dependencies for automatic files
bam.h: bam_automatic_include.h
main.o: bam_automatic_initialize.c bam_automatic_finalize.c
