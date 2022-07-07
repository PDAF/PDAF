cmake_minimum_required(VERSION 3.4)

project(psmile Fortran)

set(PDAF_ACTIVE OFF CACHE BOOL "compile fesom with PDAF-binding active")

# get our source files
set(src_home ${CMAKE_CURRENT_LIST_DIR}) # path to src directory starting from the dir containing our CMakeLists.txt
file(GLOB all_sources ${src_home}/src/*.F90)
set(all_headers ${src_home}/include)

if(${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
	set(OASIS_FFLAGS -ansi_alias -fp-model precise -no-fast-transcendentals)
endif()
# depends on the mct library
add_subdirectory(../mct/mct ${PROJECT_BINARY_DIR}/mct)

# depends on the scrip library
add_subdirectory(../scrip/src ${PROJECT_BINARY_DIR}/scrip)

include(${CMAKE_CURRENT_LIST_DIR}/../../../cmake/FindNETCDF.cmake)

# create our library (set its name to name of this project)
add_library(${PROJECT_NAME} ${all_sources})
target_link_libraries(${PROJECT_NAME} PRIVATE ${NETCDF_Fortran_LIBRARIES} mct scrip)
target_compile_definitions(${PROJECT_NAME} PRIVATE use_libMPI use_netCDF use_comm_MPI1)
if(${PDAF_ACTIVE})
  target_compile_definitions(${PROJECT_NAME} PRIVATE use_PDAF)
endif()
target_include_directories(${PROJECT_NAME}
PRIVATE ${NETCDF_Fortran_INCLUDE_DIRECTORIES} ${all_headers}
INTERFACE ${PROJECT_BINARY_DIR}
)
target_compile_options(${PROJECT_NAME} PRIVATE ${OASIS_FFLAGS})
