#Include any platform specific settings
include ./../Makefile.in

#Set the correct headers depending
#on the precision set in the Makefile.in
ifeq ($(PRECISION),DOUBLE)
FHEADER = cfd2lcs_inc_double.f90
CHEADER = cfd2lcs_inc_double.h
else
ifeq ($(PRECISION),SINGLE)
FHEADER = cfd2lcs_inc_single.f90
CHEADER = cfd2lcs_inc_single.h
endif
endif

#Set the correct hdf5_file
ifeq ($(HDF5_SUPPORT),"TRUE")
IO = io_hdf5_m.f90
else
IO = io_nohdf5_m.f90
endif

.SUFFIXES:
.SUFFIXES: .f90 .o

default: headers libbspline.a libcfd2lcs.a
	mkdir -p ./../lib
	cp libcfd2lcs.a ./../lib
	$(info ************ LIBCFD2LCS COMPILATION SUCCESSFUL *******************)

OBJS = \
	data_m.o \
	structured_m.o \
	unstructured_m.o \
	lp_m.o \
	comms_m.o \
	sgrid_m.o \
	io_m.o \
	lp_tracking_m.o \
	flowmap_m.o \
	lp_motion_m.o \
	lcs_m.o \
	cfd2lcs.o

headers:
	mkdir -p ./../include
	$(info *******Compiling libcfd2lcs for $(PRECISION) precision************)
	cp $(FHEADER) ./../include/cfd2lcs_inc.f90
	cp $(CHEADER) ./../include/cfd2lcs_inc.h
	ln -fs $(IO) io_m.f90

libcfd2lcs.a:	$(OBJS)
	$(AR) $@ $(OBJS)

libbspline.a:
	cd bspline; make

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $(CFD2LCS_INC) $*.f90

clean:
	cd bspline; make clean
	rm -f *.o *.mod libcfd2lcs.a
