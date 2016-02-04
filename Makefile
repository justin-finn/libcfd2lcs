#########################################################
# Top level Makefile for cfd2lcs
# - builds cfd2lcs libraries and some examples.
#
# Sample usage:
# % make $(PLATFORM)
# % make EXAMPLES
#
# Where $(PLATFORM) is one of the supported machines.
#
# To add your own $(PLATFORM), do the following:
#
# 1. Create a new Makefile.$(PLATFORM).in in the subdirectory
#    ./makefiles. This file gets included in both the src
#    and  examples Makefiles, and contains compiler/linker/archiver
#    names, compiler options, and library locations for the
#    various dependencies (lapack, MPI, HDF5, etc).
#
# 2. Add the relevant lines below to link and make your new
#    $(PLATFORM).
#
#########################################################

# by default, assume a Makefile.in has been provided...
default:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(cd src ; make)

# Examples
EXAMPLES:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(cd examples ; make)

# Documentation
DOC:
	(cd doc; pdflatex libcfd2lcs_manual.tex)





# AWESOMO4000: HPZ620 Linux Workstation at UoL...
AWESOMO4000:
	(ln -fs makefiles/Makefile.AWESOMO4000.in Makefile.in)
	(sed -i '/PREFIX =/c\PREFIX = '"$$PWD"'' ./makefiles/Makefile.AWESOMO4000.in)
	(make libclean)
	(cd src ; make)

# AWESOMO4000-PROFILE: HPZ620 Linux Workstation at UoL...  Parallel profiling with Scalasca/Score-p
AWESOMO4000-PROFILE:
	(ln -fs makefiles/Makefile.AWESOMO4000-PROFILE.in Makefile.in)
	(sed -i '/PREFIX =/c\PREFIX = '"$$PWD"'' ./makefiles/Makefile.AWESOMO4000-PROFILE.in)
	(make libclean)
	(cd src ; make)

# LAPPY386: Linux Mint laptop...
LAPPY386:
	(ln -fs makefiles/Makefile.LAPPY386.in Makefile.in)
	(sed -i '/PREFIX =/c\PREFIX = '"$$PWD"'' ./makefiles/Makefile.LAPPY386.in)
	(make libclean)
	(cd src ; make)

# ARCHER:  CRAY XE-6
ARCHER:
	(ln -fs makefiles/Makefile.ARCHER.in Makefile.in)
	(sed -i '/PREFIX =/c\PREFIX = '"$$PWD"'' ./makefiles/Makefile.ARCHER.in)
	(make libclean)
	(cd src ; make)






#clean:
clean:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(cd src ; make clean)
	(cd examples ; make clean)

libclean:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(cd src ; make clean)
	(cd examples ; make clean)
	rm -f ./lib/*.a
	rm -f ./include/cfd2lcs_inc*
	rm -f ./include/*.mod

distclean:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(cd src ; make clean)
	(cd examples ; make clean)
	(cd doc ; rm -f *.pdf *.aux *.log *.backup)
	rm -f ./lib/*.a
	rm -f ./include/cfd2lcs_inc*
	rm -f ./include/*.mod


