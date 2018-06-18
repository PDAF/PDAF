# $Id: Makefile 746 2009-08-04 12:16:28Z lnerger $

#######################################################
# Generic Makefile for to build PDAF with dummy model #
# To choose the architecture set $PDAF_ARCH           #
#######################################################

.SUFFIXES: .F90 .o 

# Modules used for the model part
MODULES =  	mod_parallel.o \
		mod_modeltime.o \
		timer.o \
		mod_memcount.o

MOD_PARSER =	parser_mpi.o

# Main driver routine
OBJ_DRIVER =  	main.o

# Main driver routine for PDAF offline mode
OBJ_OFFLINE =  	main_offline.o

######################################################

# Include machine-specific definitions
# For available include files see directory make.arch
# To choose a file, set PDAF_ARCH either here or by an
# environment variable.
include ../../../make.arch/$(PDAF_ARCH).h

######################################################

driver.a : $(MODULES) $(MOD_PARSER) $(OBJ_MPI) $(OBJ_DRIVER)
	@echo "++++++ Generate driver archive for PDAF testsuite ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_PARSER) \
		$(OBJ_MPI)
	@echo "++++++ Done ++++++"

driver_offline.a : $(MODULES) $(MOD_PARSER) $(OBJ_MPI) $(OBJ_OFFLINE)
	@echo "++++++ Generate driver archive for PDAF testsuite ++++++"
	$(AR) -r $(AR_SPEC) $@ $(MODULES) $(MOD_PARSER) $(OBJ_MPI)
	@echo "++++++ Done ++++++"

.F90.o :
	$(FC) $(OPT) $(MPI_INC) $(CPP_DEFS) -c $*.F90

# For older compilers one might need to separate the
# preprocessing from the compilation as defined below:
#.F90.o :
#	$(CPP) -P -C $(MPI_INC) $(CPP_DEFS) $*.F90 $*.f90
#	$(FC) $(OPT) $(MPI_INC) -c $*.F90
#	@rm -f  $*.f90

clean :
	rm -f *.o *.mod *.a

$(OBJ_DRIVER) : $(MODULES) $(OBJ_MPI)