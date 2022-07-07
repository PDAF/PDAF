cmake_minimum_required(VERSION 3.4)

# set default build type cache entry (do so before project(...) is called, which would create this cache entry on its own)
if(NOT CMAKE_BUILD_TYPE)
	message(STATUS "setting default build type: Release")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel.")
endif()

#set(PDAF_ACTIVE ON CACHE BOOL "compile fesom with PDAF-binding active")

project(AWI-CM)

add_subdirectory(oasis3-mct/lib/psmile)
add_subdirectory(fesom_cpl)
add_subdirectory(echam-6.3.02p4/src)
#add_subdirectory(echam-6.3.04p1/src/echam)
