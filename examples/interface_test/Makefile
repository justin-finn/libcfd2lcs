###############
#Example Makefile for compiling with libcfd2lcs.

#Change location of CFD2LCS_HOME on your system:
CFD2LCS_HOME = /home/finnj/work/CFD/cfd2lcs/code/libcfd2lcs
##############

include $(CFD2LCS_HOME)/Makefile.in
CFD2LCS_LIB = $(CFD2LCS_HOME)/lib
CFD2LCS_MOD = $(CFD2LCS_HOME)/mod

.SUFFIXES:
.SUFFIXES: .f90 .o

default: CFD2LCS

USER_OBJS =
CFD2LCS_OBJS =  $(USER_OBJS) cfd2lcs.o
CFD2LCS_LIBS = -L$(CFD2LCS_LIB) -lcfd2lcs $(LAPACK_LIBS)

CFD2LCS: $(CFD2LCS_OBJS)
	$(F90) -o $@ $(CFD2LCS_OBJS) $(CFD2LCS_LIBS)

#Rule for objects
.f90.o:
	$(F90) $(FFLAGS) -$(INC)$(CFD2LCS_MOD) $*.f90

clean:
	rm -f *.o *.mod CFD2LCS
