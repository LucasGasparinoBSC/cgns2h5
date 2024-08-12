cmake_minimum_required(VERSION 3.20)

# Options
option(USE_MN5 "Compiling in MareNostrum5" OFF)
if(USE_MN5)
    message("Linking HDF5 as HDF5::HDF5...")
endif(USE_MN5)

set(ENV{CC} "mpicc")
set(ENV{CXX} "mpicxx")

# CMake extra file are on this path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Include the extra file
include(gpu)

project(tool_cgns2h5 LANGUAGES C CXX)

# Find HDF5
find_package(HDF5 REQUIRED COMPONENTS C HL)

# Find CGNS
find_package(CGNS REQUIRED)

# Set CMake standard
set(CMAKE_C_STANDARD 17)

# Set libraries to shared
set(LIBRARY_TYPE SHARED)

# If no build type is specified, default to Release
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Set C/CXX flags
# NVIDIA compiler
if(CMAKE_C_COMPILER_ID STREQUAL "NVHPC" OR CMAKE_C_COMPILER_ID STREQUAL "PGI")
    # Use nvidia-smi to find the compute capability of your GPU
    set_gpu()
    set(CMAKE_C_FLAGS "-acc -cuda -Minfo=all -gpu=cc${GPU_CC},cuda${CUDA_VERSION},lineinfo,managed -lnvToolsExt")
    set(CMAKE_CXX_FLAGS "-acc -cuda -Minfo=all -gpu=cc${GPU_CC},cuda${CUDA_VERSION},lineinfo,managed -lnvToolsExt")
    # Debug flags
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -traceback")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -traceback")
    # Release flags
    set(CMAKE_C_FLAGS_RELEASE "-fast")
    set(CMAKE_CXX_FLAGS_RELEASE "-fast")
    add_compile_definitions("-DUSE_GPU -lnvToolsExt")
# GNU compiler
elseif(CMAKE_C_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_C_FLAGS "-fopenmp")
    set(CMAKE_CXX_FLAGS "-fopenmp")
    # Debug flags
    set(CMAKE_C_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Werror")
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -Wall -Wextra -Wpedantic -Werror")
    # Release flags
    set(CMAKE_C_FLAGS_RELEASE "-O3")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3")
endif()

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# add the automatically determined parts of the RPATH
# which point to directories outside the build tree to the install RPATH
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# the RPATH to be used when installing, but only if it's not a system directory
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${isSystemDir}" STREQUAL "-1")

add_subdirectory(src)
