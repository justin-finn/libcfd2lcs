###############
CFD2LCS_HOME = /home/finnj/work/CFD/cfd2lcs/code/libcfd2lcs

#Example Makefile for compiling an application with libcfd2lcs.
#Everything below here should remain unchanged:
##############

include $(CFD2LCS_HOME)/Makefile.in

CFD2LCS_LIBS = -L $(CFD2LCS_HOME)/lib/ -lcfd2lcs

.SUFFIXES:
.SUFFIXES: .f90 .o

default: CFD2LCS_TEST

OBJS =  interface_test.o
LIBS = $(CFD2LCS_LIBS) $(THIRD_PARTY_LIBS)

CFD2LCS_TEST: $(OBJS)
	$(F90) -o $@ $(OBJS) $(LIBS)
	mkdir -p dump

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $(THIRD_PARTY_INC) $*.f90

clean:
	rm -f *.o *.mod CFD2LCS_TEST
