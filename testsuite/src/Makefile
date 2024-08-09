# $Id: Makefile 771 2009-11-24 16:06:08Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

.SUFFIXES: .F90 .o 

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../../make.arch/$(PDAF_ARCH).h

######################################################
info:
	@echo "Makefile to build PDAF examples";
	@echo  "------------------------------------------------------------------";
	@echo  "Use as   make TARGET   where TARGET can be:";
	@echo  "  ----- BUILD -----------------------------------------------------";
	@echo  "  info               - This text";
	@echo  "  dummymodel         - Dummy model without PDAF";
	@echo  "                       (Do not use -DUSE_PDAF for this target!)"
	@echo  "  pdaf_dummy_online  - PDAF online configuration with dummy model";
	@echo  "  pdaf_dummy_offline - PDAF in offline configuration";
	@echo  "  pdaf_dummy_si      - PDAF online with dummy model & simplified interface";
	@echo  "  pdaf_dummy_snglprec- PDAF online with dummy model in single precision";
	@echo  "  ----- CLEAN -----------------------------------------------------";
	@echo  "  clean              - Clean up objects in test suite";
	@echo  "  cleanall           - Clean up test suite";
	@echo  "  cleanmodels        - Clean up model directories";
	@echo  "  cleanpdaf          - Clean up directory of PDAF library";
	@echo  "  distclean          - Combined cleanpdaf and clean (cleans all)";
	@echo  "  ----- TEST ------------------------------------------------------";
	@echo  "  test_pdaf_online   - Run serial tests for PDAF in online config.";
	@echo  "  test_pdaf_offline  - Run serial tests for PDAF in offline config.";
	@echo  "  test_pdaf_smoother - Run serial tests for PDAF smoother in online config.";
	@echo  "  NOTE: Tests are executed in the directory ../bin";
	@echo  "";
	@echo  "Depending on whether PDAF should be active, one has to define";
	@echo  "USE_PDAF in the include file in make.arch/ file for preprocessing.";
	@echo  "List available settings for PDAF_ARCH with 'make listarch'";
	@echo  "Use machine-specific definitions for PDAF_ARCH: " $(PDAF_ARCH);
	@echo  "------------------------------------------------------------------";

######################################################

driver :
	@echo "++++++ Build driver ++++++"
	@cd main; make driver.a

driver_offline :
	@echo "++++++ Build off line driver ++++++"
	@cd main; make driver_offline.a

cleandriver:
	@echo "+++ Clean up directory main"
	@cd main; make clean

######################################################

dummymodel : driver dummymodel_1D.a
	@echo "++++++ Linking ++++++"
	$(LD) $(OPT_LNK) -o ../bin/$@ \
	main/main.o dummymodel_1D/dummymodel_1D.a main/driver.a \
	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

dummymodel_1D.a:
	@echo "++++++ Make 1D dummymodel ++++++"
	@cd dummymodel_1D; make $@

cleandummymodel:
	@echo "+++ Clean up directory dummymodel_1D"
	@cd dummymodel_1D; make clean

######################################################

pdaf_dummy_online : libpdaf-d.a driver dummymodel_1D_pdaf.a
	@echo "++++++ Linking ++++++"
	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
	main/main.o dummymodel_1D/dummymodel_1D_pdaf.a main/driver.a \
	-L../../lib -lpdaf-d \
	$(NC_LIB) $(LINK_LIBS)
	@echo "++++++ Done ++++++"

dummymodel_1D_pdaf.a:
	@echo "++++++ Make 1D dummymodel and user routines for PDAF ++++++"
	@cd dummymodel_1D; make $@

######################################################

pdaf_dummy_snglprec : libpdaf-d.a driver dummymodel_1D_sngl_pdaf.a
	@echo "++++++ Linking ++++++"
	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
	main/main.o dummymodel_1D_snglprec/dummymodel_1D_sngl_pdaf.a main/driver.a \
	-L../../lib -lpdaf-d \
	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

dummymodel_1D_sngl_pdaf.a:
	@echo "++++++ Make 1D dummymodel and user routines for PDAF ++++++"
	@cd dummymodel_1D_snglprec; make $@

cleandummysngl:
	@echo "+++ Clean up directory dummymodel_1D_snglprec"
	@cd dummymodel_1D_snglprec; make clean

######################################################

pdaf_dummy_offline : libpdaf-d.a driver_offline dummymodel_1D_pdaf-ol.a
	@echo "++++++ Linking ++++++"
	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
	main/main_offline.o \
	offline_1D/dummymodel_1D_pdaf-ol.a main/driver_offline.a \
	-L../../lib -lpdaf-d \
	$(LINK_LIBS)
	@echo "++++++ Done ++++++"

dummymodel_1D_pdaf-ol.a:
	@echo "++++++ Make 1D dummymodel and user routines for PDAF offline mode ++++++"
	@cd offline_1D; make $@

cleanoffline:
	@echo "+++ Clean up directory offline_1D"
	@cd offline_1D; make clean

######################################################

pdaf_dummy_nobs : libpdaf-d.a driver dummymodel_1D_pdaf_nobs.a
	@echo "++++++ Linking ++++++"
	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
	main/main.o dummymodel_1D_nobs/dummymodel_1D_pdaf_nobs.a main/driver.a \
	-L../../lib -lpdaf-d \
	$(NC_LIB) $(LINK_LIBS)
	@echo "++++++ Done ++++++"

dummymodel_1D_pdaf_nobs.a:
	@echo "++++++ Make 1D dummymodel and user routines for PDAF nobs case ++++++"
	@cd dummymodel_1D_nobs; make $@

cleandummynobs:
	@echo "+++ Clean up directory dummymodel_1D_nobs"
	@cd dummymodel_1D_nobs; make clean

######################################################

pdaf_dummy_si : libpdaf-d.a driver dummymodel_1D_pdaf_si.a
# Here, we don't link the PDAF library as a library, but as an archive file. This was
# required for the Cray compiler CCE, at least until version 8.3.7, because
# is didn't find the call-back routines which are decleared 'external' in the 
# PDAF core routine.
	@echo "++++++ Linking ++++++"
	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
	main/main.o dummymodel_1D_si/dummymodel_1D_pdaf_si.a main/driver.a \
	../../lib/libpdaf-d.a \
	dummymodel_1D_si/dummymodel_1D_pdaf_si.a ../../lib/libpdaf-d.a \
	$(LINK_LIBS)
# Linking libpdaf-d.a explicitly as a library is more consistent, but didn't work
# with the Cray compiler CCE, at least until version 8.3.7.
#	$(LD)  $(OPT_LNK)  -o ../bin/$@  \
#	main/main.o dummymodel_1D_si/dummymodel_1D_pdaf_si.a main/driver.a \
#	-L../../lib -lpdaf-d \
#	dummymodel_1D_si/dummymodel_1D_pdaf_si.a -lpdaf-d \
#	$(LINK_LIBS)
	echo "++++++ Done ++++++"

dummymodel_1D_pdaf_si.a:
	@echo "++++++ Make 1D dummymodel and user routines for PDAF ++++++"
	@cd dummymodel_1D_si; make $@

cleandummy1Dsi:
	@echo "+++ Clean up directory dummymodel_1D_si"
	@cd dummymodel_1D_si; make clean

######################################################

libpdaf-d.a: 
	@echo "++++++ Generate Filter library ++++++"
	@cd ../..; make;

.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) $(CPP_PDAF2) -Ioffline_1D -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

######################################################
# Cleans

clean: cleandriver cleanmodels

cleanall: cleandriver cleanbin cleanmodels

cleanbin:
	@echo "+++ Clean up bin directory"
	@cd ../bin; rm -f dummymodel pdaf_dummy_online pdaf_dummy_offline \
	  pdaf_dummy_si pdaf_dummy_snglprec

cleanmodels: cleandummymodel cleanoffline cleandummy1Dsi cleandummysngl

cleanpdaf:
	@echo "+++ Clean up PDAF directory"
	cd ../..; make clean

distclean: cleanpdaf cleanall

######################################################
# List arch files

listarch:
	@echo Available architecture-specific input files for PDAF_ARCH
	@echo ---------------------------------------------------------
	@ls -1 ../../make.arch | cut -d"." -f1

######################################################
# Tests

test_pdaf_online:
	@echo "Running tests with PDAF in online configuration in ../bin"
	@echo "The output files can be compared with outputs in testsuite/tests_dummy1D"
	cd ../bin; ../tests_dummy1D/runtests_1pe.sh

test_pdaf_offline:
	@echo "Running tests with PDAF in offline configuration in ../bin"
	@echo "The output files can be compared with outputs in testsuite/tests_dummy1D"
	cd ../bin; ../tests_dummy1D/runtests_offline_1pe.sh 

test_pdaf_smoother:
	@echo "Running tests with PDAF smoother in online configuration in ../bin"
	@echo "The output files can be compared with outputs in testsuite/tests_dummy1D"
	cd ../bin; ../tests_dummy1D/runtests_smoother_1pe.sh


