project(cgns2h5)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.cxx)
set(HEADER_DIR ${CMAKE_CURRENT_LIST_DIR}/inc)
file(GLOB HEADER_FILES ${HEADER_DIR}/*.h)

#include_directories(${HEADER_DIR} ${lib1_INCLUDE_DIRS})
include_directories(${HEADER_DIR})

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})
set_property(TARGET ${PROJECT_NAME} PROPERTY LINKER_LANGUAGE CXX)
set_mpi()
set_hdf5()
set_cgns()

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "CGNS2H5")

install(TARGETS ${PROJECT_NAME} DESTINATION bin)
install(FILES ${HEADER_FILES} DESTINATION include)