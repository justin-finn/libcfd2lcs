#########################################################
# Top level Makefile for cfd2lcs
# - builds cfd2lcs libraries and some examples.
#
# Sample usage:
# % make PLATFORM
#
# Where PLATFORM is one of the supported machines.
#
# To add your own PLATFORM, do the following:
#
# 1. Create a new Makefile.PLATFORM.in in the subdirectory
#    ./makefiles. This file gets included in both the src
#    and  examples Makefiles, and contains compiler/linker/archiver
#    names, compiler options, and library locations for the
#    various dependencies (lapack, MPI, etc).
#
# 2. Add the relevant lines below to link and make your new
#    PLATFORM.
#
#########################################################

# by default, assume a Makefile.in has been provided...
default:
	(cd src ; make)

# AWESOMO4000: HPZ620 Linux Workstation at UoL...
AWESOMO4000:
	ln -fs makefiles/Makefile.AWESOMO4000.in Makefile.in
	(cd src ; make)

clean:
	(cd src ; make clean)

libclean:
	(cd src ; make libclean)
