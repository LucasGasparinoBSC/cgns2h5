message("-- Configuring MPI ...")
find_package(MPI REQUIRED)
if (NOT MPI_FOUND)
    #do nothing
    message(FATAL_ERROR "MPI not found!")
endif()

function(set_mpi)
    target_link_libraries(${PROJECT_NAME} MPI::MPI_CXX MPI::MPI_C)
endfunction()