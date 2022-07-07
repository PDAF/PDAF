cmake_minimum_required(VERSION 3.4)

project(echam6 C Fortran)

set(PDAF_ACTIVE OFF CACHE BOOL "compile fesom with PDAF-binding active")

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB sources_Fortran ${src_home}/*.f90)
file(GLOB sources_C ${src_home}/version.c)
if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
	set(ECHAM_FFLAGS -r8 -ansi_alias -fp-model strict -no-fast-transcendentals)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/../../cmake/FindNETCDF.cmake)

add_subdirectory(../support ${PROJECT_BINARY_DIR}/support)
add_subdirectory(../cdi/src ${PROJECT_BINARY_DIR}/cdi)
add_subdirectory(../lapack ${PROJECT_BINARY_DIR}/lapack)
add_subdirectory(../blas ${PROJECT_BINARY_DIR}/blas)

add_library(${PROJECT_NAME}_C ${sources_C})
target_compile_definitions(${PROJECT_NAME}_C PRIVATE pgiFortran __NOUSERNAME__)

# create our binary (set its name to name of this project)
add_executable(${PROJECT_NAME} ${sources_Fortran})
target_include_directories(${PROJECT_NAME} PRIVATE ${src_home}/../include ${PROJECT_BINARY_DIR} ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${PROJECT_BINARY_DIR}/cdi/yaxt)
target_link_libraries(${PROJECT_NAME} ${PROJECT_NAME}_C ${NETCDF_Fortran_LIBRARIES} ${NETCDF_C_LIBRARIES} support cdi lapack blas psmile yaxt)
target_compile_definitions(${PROJECT_NAME} PRIVATE use_comm_MPI1 __prism __cpl_feom __oa3mct USE_MPI HAVE_YAXT HAVE_CDIPIO INPUT_IN_ECHAM)
if(${PDAF_ACTIVE})
   target_compile_definitions(${PROJECT_NAME} PRIVATE use_PDAF)
   target_link_libraries(${PROJECT_NAME} $ENV{PDAF_DIR}/lib/libpdaf-d.a)
   target_include_directories(${PROJECT_NAME} PRIVATE $ENV{PDAF_DIR}/src)
endif()
# be sure we use the preprocessor
if(NOT DEFINED ENV{FFLAGS})
	# CMAKE_Fortran_COMPILER_ID will also work if a wrapper is being used (e.g. mpif90 wraps ifort -> compiler id is Intel)
	if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
   	target_compile_options(${PROJECT_NAME} PRIVATE -fpp -traceback)
	elseif(${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
   	target_compile_options(${PROJECT_NAME} PRIVATE -cpp -ffree-line-length-none)
   endif()
   target_compile_options(${PROJECT_NAME} PRIVATE ${ECHAM_FFLAGS})
endif()

set(ECHAM_INSTALL_FILEPATH "~/fesom_echam6_oasis3-mct/echam-6.3.02p4/bin/echam6" CACHE FILEPATH "file path where the ECHAM binary should be put")
get_filename_component(ECHAM_INSTALL_PATH ${ECHAM_INSTALL_FILEPATH} DIRECTORY)
get_filename_component(ECHAM_INSTALL_NAME ${ECHAM_INSTALL_FILEPATH} NAME)
install(PROGRAMS ${PROJECT_BINARY_DIR}/${PROJECT_NAME} DESTINATION ${ECHAM_INSTALL_PATH} RENAME ${ECHAM_INSTALL_NAME})
