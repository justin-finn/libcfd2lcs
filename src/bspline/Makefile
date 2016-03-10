include ./../../Makefile.in

#Define the headers:
FHEADER_DP = ./../cfd2lcs_inc_double.f90
FHEADER_SP = ./../cfd2lcs_inc_single.f90

#Objects:
OBJS = bspline_sub_module.o bspline_oo_module.o bspline_module.o

.SUFFIXES:
.SUFFIXES: .f90 .o

default:
	$(info ERROR: Compile with either:  1. make sp  2. make dp)

sp: headers_sp libbspline_sp.a
	mkdir -p ./../../lib
	mv libbspline_sp.a ./../../lib
	$(info ********  libbspline_sp.a COMPILATION SUCCESSFUL ***********)

dp: headers_dp libbspline_dp.a
	mkdir -p ./../../lib
	mv libbspline_dp.a ./../../lib
	$(info ********  libbspline_dp.a COMPILATION SUCCESSFUL ***********)

headers_sp:
	$(info ************  COMPILING libbspline_sp.a ********************)
	ln -fs $(FHEADER_SP) cfd2lcs_inc.f90

headers_dp:
	$(info ************  COMPILING libbspline_dp.a ********************)
	ln -fs $(FHEADER_DP) cfd2lcs_inc.f90

libbspline_sp.a libbspline_dp.a:	$(OBJS)
	$(AR) $@ $(OBJS)

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $*.f90

clean:
	rm -f *.o *.mod *.a
