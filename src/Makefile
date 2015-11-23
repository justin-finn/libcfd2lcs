include ./../Makefile.in

.SUFFIXES:
.SUFFIXES: .f90 .o

default: libbspline.a libcfd2lcs.a
	mkdir -p ./../lib
	mkdir -p ./../include
	cp libcfd2lcs.a ./../lib
	cp cfd2lcs_inc* ./../include
	$(info ************ LIBCFD2LCS COMPILATION SUCCESSFUL ************)

OBJS = \
	data_m.o \
	structured_m.o \
	unstructured_m.o \
	lp_m.o \
	comms_m.o \
	sgrid_m.o \
	io_m.o \
	lp_tracking_m.o \
	lp_motion_m.o \
	lcs_m.o \
	cfd2lcs.o

libcfd2lcs.a:	$(OBJS)
	$(AR) $@ $(OBJS)

libbspline.a:
	cd bspline; make

.f90.o:
	$(F90) $(FFLAGS) -$(INC). -$(INC)./../include $(THIRD_PARTY_INC) $*.f90

clean:
	cd bspline; make clean
	rm -f *.o *.mod libcfd2lcs.a

libclean:
	cd bspline; make clean
	rm -f *.o *.mod libcfd2lcs.a
	rm -f ./../lib/*.a
	rm -f ./../include/cfd2lcs_inc*
	rm -f ./../include/*.mod
