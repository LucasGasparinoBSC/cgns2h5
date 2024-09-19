message("-- Configuring HDF5...")
if(DEFINED ENV{CGNS_ROOT})
  message("-- CGNS_ROOT properly set in environmental var:" $ENV{CGNS_ROOT})
elseif(DEFINED ENV{CGNS_DIR})
  message("-- CGNS_DIR properly set in environmental var:" $ENV{CGNS_DIR})
elseif(CGNS_ROOT)
  message("-- CGNS_ROOT properly set in cmake:" ${CGNS_ROOT})
else()
  message("-- POSSIBLE ERROR! Environmental variables CGNS_ROOT/CGNS_DIR not properly set!")
  message("-- If Cmake does not find CGNS package you have two options:")
  message("--   1. Define environmental variables CGNS_ROOT and/or CGNS_DIR")
endif()

# find the parallel CGNS library
find_package(CGNS REQUIRED)

# Function for linking against the parallel CGNS lib
function(set_cgns)
    include_directories(${CGNS_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} cgns)
endfunction()