###############
CFD2LCS_HOME = /home/finnj/work/CFD/cfd2lcs/code/libcfd2lcs
include $(CFD2LCS_HOME)/Makefile.in

NETCDF_LIBS = -L/usr/lib/ -lnetcdff
NETCDF_INC = -I/usr/include/
##############

.SUFFIXES:
.SUFFIXES: .c .o
.SUFFIXES: .f90 .o

default:  ROMS_LCS
	$(info ************  COMPILATION SUCCESSFUL ************)


OBJS = roms_lcs.o

ROMS_LCS: $(OBJS)
	$(F90) -o $@ $(OBJS) $(NETCDF_LIBS) $(CFD2LCS_SP_LIBS)

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $(CFD2LCS_INC) $(NETCDF_INC) $*.f90


clean:
	rm -f *.o *.mod ROMS_LCS

dataclean:
	rm -rf cfd2lcs_output cfd2lcs_tmp
