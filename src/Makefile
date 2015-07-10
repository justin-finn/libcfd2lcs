CFD2LCS_LIB = /home/finnj/CFD/cfd2lcs/code/libcfd2lcs/
CFD2LCS_LIB = $(CFD2LCS_HOME)/lib/
CFD2LCS_MOD = $(CFD2LCS_HOME)/mod/

include ./../Makefile.in

.SUFFIXES:
.SUFFIXES: .f90 .o

default:	libcfd2lcs.a
	mkdir -p $(CFD2LCS_LIB)
	mkdir -p $(CFD2LCS_MOD)
	cp libcfd2lcs.a $(CFD2LCS_LIB)
	cp *.mod $(CFD2LCS_MOD)

OBJS = precision_m.o data_m.o

libcfd2lcs.a:	$(OBJS)
	$(AR) $@ $(OBJS)

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $*.f90

clean:
	rm -f *.o *.mod libcfd2lcs.a
