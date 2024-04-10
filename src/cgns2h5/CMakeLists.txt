project(cgns2h5)

file(GLOB_RECURSE SRC_FILES ${CMAKE_CURRENT_LIST_DIR}/*.cxx)
set(HEADER_DIR ${CMAKE_CURRENT_LIST_DIR}/inc)
file(GLOB HEADER_FILES ${HEADER_DIR}/*.h)

#include_directories(${HEADER_DIR} ${lib1_INCLUDE_DIRS})
include_directories(${HEADER_DIR} ${mpiWrap_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SRC_FILES} ${HEADER_FILES})
target_link_libraries(${PROJECT_NAME} mpiWrap)

set_target_properties(${PROJECT_NAME} PROPERTIES OUTPUT_NAME "CGNS2H5")

install(TARGETS ${PROJECT_NAME} DESTINATION bin)