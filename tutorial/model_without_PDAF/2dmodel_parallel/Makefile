#######################################################
# Simple Makefile to build the 2D tutotrial model     #
#######################################################

# Include machine-specific definitions
include compile_settings.txt

#######################################################
# List of sources

SRC_MODEL =  	model_parallel_mod.F90 \
		model_mod.F90 \
		model_io_mod.F90 \
		model_init_mod.F90 \
		model_step_mod.F90 \
		model_post_mod.F90 \
		main.F90

#######################################################
# Define object files - to be create in build/
OBJ_MODEL := $(SRC_MODEL:%.F90=build/%.o)

######################################################

.phony: info
info:
	@echo "Makefile to build advanced 2d tutorial model of PDAF";
	@echo "------------------------------------------------------------------";
	@echo "Use as   make TARGET   where TARGET can be:";
	@echo "  model         - the 2D tutorial model";
	@echo "  clean         - clean up compilation";
	@echo "  cleandata     - remove output files";
	@echo "------------------------------------------------------------------";


model : $(OBJ_MODEL)
	$(LD)  $(FFLAGS)  -o $@  $(OBJ_MODEL) $(LIBS)

######################################################
# use pattern rule to create rules for all object files

$(OBJ_MODEL): build/%.o: src/%.F90
	$(FC) $(FFLAGS) $(INCLUDES) -c -o $@ $(MODULEOPT) build/ $<

######################################################
# Cleans

.phony: cleanall clean cleandata
cleanall : clean cleandata

clean :
	rm -f build/*.o build/*.mod model

cleandata:
	rm -f field*.txt


######################################################
# Dependencies
# (created using 'mkdepends' in PDAF's directory external/)

build/main.o: ./src/main.F90 build/model_parallel_mod.o build/model_init_mod.o build/model_step_mod.o build/model_post_mod.o
build/model_init_mod.o: ./src/model_init_mod.F90 build/model_mod.o build/model_parallel_mod.o build/model_io_mod.o
build/model_io_mod.o: ./src/model_io_mod.F90 build/model_parallel_mod.o build/model_mod.o
build/model_mod.o: ./src/model_mod.F90
build/model_parallel_mod.o: ./src/model_parallel_mod.F90
build/model_post_mod.o: ./src/model_post_mod.F90 build/model_mod.o build/model_parallel_mod.o
build/model_step_mod.o: ./src/model_step_mod.F90 build/model_mod.o build/model_parallel_mod.o build/model_io_mod.o
