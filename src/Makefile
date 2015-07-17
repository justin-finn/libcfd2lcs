include ./../Makefile.in

.SUFFIXES:
.SUFFIXES: .f90 .o

default:	libcfd2lcs.a
	mkdir -p ./../lib
	mkdir -p ./../mod
	cp libcfd2lcs.a ./../lib
	cp *.mod ./../mod

OBJS = precision_m.o mpi_m.o data_m.o io_m.o cfd2lcs.o

libcfd2lcs.a:	$(OBJS)
	$(AR) $@ $(OBJS)

.f90.o:
	$(F90) $(FFLAGS) -$(INC). $(HDF5_INC) $(MPI_INC) $*.f90

clean:
	rm -f *.o *.mod libcfd2lcs.a

libclean:
	rm -f *.o *.mod libcfd2lcs.a
	rm -f ./../lib/*.a
	rm -f ./../mod/*.mod
