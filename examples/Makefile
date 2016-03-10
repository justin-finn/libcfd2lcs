###############
#Example Makefile for a c or fortran program linking to libcfd2lcs
#All dependencies and include paths are stored in the variables
#CFD2LCS_LIBS and CFD2LCS_INC in the Makefile.in of the top-level
#cfd2lcs directory.  Modify the line here to point to it on your system.
CFD2LCS_HOME = /home/finnj/work/CFD/cfd2lcs/code/libcfd2lcs
include $(CFD2LCS_HOME)/Makefile.in
##############


.SUFFIXES:
.SUFFIXES: .c .o
.SUFFIXES: .f90 .o

default:  DOUBLE_GYRE_F90 DOUBLE_GYRE_C ABC_FLOW_F90 ABC_FLOW_C
	$(info ************  COMPILATION SUCCESSFUL ************)

DOUBLE_GYRE_F90: double_gyre_f.o
	$(F90) -o $@ double_gyre_f.o $(CFD2LCS_SP_LIBS)

DOUBLE_GYRE_C: double_gyre_c.o
	$(CC) -o $@ double_gyre_c.o $(CFD2LCS_SP_LIBS) $(CLIBS)

ABC_FLOW_F90: abc_flow_f.o
	$(F90) -o $@ abc_flow_f.o $(CFD2LCS_SP_LIBS)

ABC_FLOW_C: abc_flow_c.o
	$(CC) -o $@ abc_flow_c.o $(CFD2LCS_SP_LIBS) $(CLIBS)


.c.o:
	$(CC) $(CFLAGS) -$(INC). $(CFD2LCS_INC) $*.c

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $(CFD2LCS_INC) $*.f90

clean:
	rm -f *.o *.mod DOUBLE_GYRE_F90 DOUBLE_GYRE_C ABC_FLOW_F90 ABC_FLOW_C

dataclean:
	rm -rf cfd2lcs_output cfd2lcs_tmp
