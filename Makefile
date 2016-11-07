#
#Copyright (C) 2015-2016, Justin R. Finn.  All rights reserved.
#libcfd2lcs is distributed is under the terms of the GNU General Public License
#
# Top level Makefile for libcfd2lcs
# - builds libcfd2lcs libraries and some example programs.
#
# Sample usage:
# % make $(PLATFORM)
# % make EXAMPLES
# % make DOC
#
# Where $(PLATFORM) is one of the supported machines.
#########################################################
date=$(shell date +%y%m%d)
version=1.0

#Default:  Tell the user what to do:
default:
	$(info ----------------- libcfd2lcs build instructions-- --------------------------)
	$(info A. Sample usage of this makefile:                                           )
	$(info >>    make PLATFORM             [Builds library for machine called PLATFORM])
	$(info >>    make EXAMPLES             [Builds example programs]                   )
	$(info >>    make DOC                  [Builds documentation, requires pdflatex]   )
	$(info                                                                             )
	$(info B. To build libcfd2lcs on a currently supported platform:                   )
	$(info >>    make AWESOMO4000          [UoL workstation]                           )
	$(info >>    make AWESOMO4000-PROFILE  [UoL workstation with Scalasca profiler]    )
	$(info >>    make LAPPY386             [Linux Mint laptop]                         )
	$(info >>    make ARCHER               [Cray XE-6]                                 )
	$(info >>    make BLUECRYSTAL          [U of Bristol HPC]                          )
	$(info                                                                             )
	$(info C. To Compile on a new platform:  Edit the variables in the file            )
	$(info ./makefiles/Makefile.YOUR_NEW_PLATFORM.in to match your system config. Then:)
	$(info >>   make YOUR_NEW_PLATFORM                                                 )
	$(info ----------------------------------------------------------------------------)

# Examples (Include some which may or may not be distributed)
EXAMPLES:
	(sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/Makefile)
	(if [ -d "./examples/mobile" ]; then sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/mobile/mobile_src/Makefile    ; fi)
	(if [ -d "./examples/roms" ]; then sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/roms/Makefile    ; fi)
	(if [ -d "./examples/ttrack3D" ]; then sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/ttrack3D/src/Makefile    ; fi)
	(if [ -d "./examples/cgs_dem/settling" ]; then sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/cgs_dem/settling/Makefile    ; fi)
	(if [ -d "./examples/cgs_dem/tgv" ]; then sed -i '/CFD2LCS_HOME =/c\CFD2LCS_HOME = '"$$PWD"'' ./examples/cgs_dem/tgv/Makefile    ; fi)
	(cd examples ; make)

# Documentation
DOC:
	(cd doc/src; pdflatex libcfd2lcs_manual.tex; bibtex libcfd2lcs_manual; pdflatex libcfd2lcs_manual.tex; pdflatex libcfd2lcs_manual.tex; mv libcfd2lcs_manual.pdf ./../)


###############################
#CURRENTLY SUPPORTED PLATFORMS
###############################

# AWESOMO4000: HPZ620 Linux Workstation at UoL...
AWESOMO4000:
	(ln -fs makefiles/Makefile.AWESOMO4000.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.AWESOMO4000.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)

# AWESOMO4000-PROFILE: HPZ620 Linux Workstation at UoL...  Parallel profiling with Scalasca/Score-p
AWESOMO4000-PROFILE:
	(ln -fs makefiles/Makefile.AWESOMO4000-PROFILE.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.AWESOMO4000-PROFILE.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)

# LAPPY386: Linux Mint laptop...
LAPPY386:
	(ln -fs makefiles/Makefile.LAPPY386.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.LAPPY386.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)

# ARCHER:  CRAY XE-6
ARCHER:
	(ln -fs makefiles/Makefile.ARCHER.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.ARCHER.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)

# BLUECRYSTAL
BLUECRYSTAL:
	(ln -fs makefiles/Makefile.BLUECRYSTAL.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.BLUECRYSTAL.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)

# YOUR_NEW_PLATFORM
YOUR_NEW_PLATFORM:
	(ln -fs makefiles/Makefile.YOUR_NEW_PLATFORM.in Makefile.in)
	(sed -i '/CFD2LCS_PREFIX =/c\CFD2LCS_PREFIX = '"$$PWD"'' ./makefiles/Makefile.YOUR_NEW_PLATFORM.in)
	(make libclean)
	(cd src ; make clean; make sp; make clean; make dp)



########
#cleanup:
########
tarball:
	(cd ./../; tar -czvf libcfd2lcs_$(date)_$(version).tar.gz libcfd2lcs)

tarball-release:
	(cd ./../; tar -zcvf libcfd2lcs_$(date)_$(version)_release.tar.gz libcfd2lcs --exclude=libcfd2lcs/examples/mobile --exclude=libcfd2lcs/examples/roms/inputData --exclude=libcfd2lcs/examples/ttrack3D --exclude=libcfd2lcs/examples/cgs_dem)

clean:
	(cd src ; make clean)
	(cd examples ; make clean)

libclean:
	(cd src ; make clean)
	rm -f ./lib/*.a
	rm -f ./include/*.f90
	rm -f ./include/*.h

distclean:
	(cd src ; make clean)
	(cd examples ; make clean; make dataclean)
	(cd examples/roms ; make clean; make dataclean)
	(cd examples/mobile ; make clean; make dataclean)
	(cd examples/ttrack3D ; make clean; make dataclean)
	(cd examples/cgs_dem ; make clean; make dataclean)
	(cd doc/src ; rm -f *.pdf *.aux *.log *.backup *.bak *.bbl *.blg *.out)
	rm -f ./lib/*.a
	rm -f ./include/cfd2lcs_inc*
