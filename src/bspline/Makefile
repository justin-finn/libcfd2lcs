include ./../../Makefile.in

.SUFFIXES:
.SUFFIXES: .f90 .o

default:	libbspline.a
	$(info ************  COMPILING libbspline.a ********************)
	mkdir -p ./../../lib
	mkdir -p ./../../include
	cp libbspline.a ./../../lib
	cp *.mod ./../../include
	$(info ********  libbspline.a COMPILATION SUCCESSFUL ***********)

OBJS = bspline_sub_module.o bspline_oo_module.o bspline_module.o

libbspline.a:	$(OBJS)
	$(AR) $@ $(OBJS)

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $*.f90

clean:
	rm -f *.o *.mod libbspline.a
