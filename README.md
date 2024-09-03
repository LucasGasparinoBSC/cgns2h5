# CGNS2H5

Small application to convert CGNS files to HDF5 files (.h5) in the format of SOD2D.

## Dependencies

- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk) *Optional*
- [MPI](https://www.open-mpi.org)
- [CGNS](https://cgns.github.io)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5)

The HDF5 must be compiled in parallel using the MPI compile wrapper. CGNS must be compiled with the HDF5 support, also with parallel support enabled.

## Usage

TL;DR:
```bash
$ mpirun -np <NumProcs> ./cgns2h5 <input.cgns> <output.h5> <x.y scale factor>
```

Once the binary is built, the command line displayed above can be used to convert a CGNS file to a HDF5 file. The `<x.y scale factor>` is a float number that scales the coordinates of the mesh. All parameters are mandatory! both the input and output files must be provided in their full names, including the extensions.

Given that the code is fully parallel, it is possible to specify `<NumProc> > 1`. Alternatively, in machines that support SLURM, a parallel job script can also be used.

The code also supports multi-GPU runs; this requires that:

- NVIDIA GPUs are present in the system;
- Both HDF5 and CGNS are compiled with the NVIDIA HPC SDK;
- the code is compiled with the NVIDIA HPC SDK;

At the moment, GPU I/O is not enabled, so only the conversion methods are ported to GPU.

## Building

For now, the code can be built either with the GNU or NVHPC compilers. In either case, proper MPI support is necessary (OpenMPI for GNU, HPCX for NVHPC are suggested combinations).

Given that the HDF5 and CGNS dependencies can present some compilation issues, particularly with NVHPC, the following sections present the author's tried and tested compilation method.

Usage of a modulefile system is recommended, as it simplifies the environment setup, and ensures no conflicts with other software.

### HDF5

- Version: 1.14.4
- GNU + OpenMPI:
```bash
$ CPP=cpp CC=mpicc FC=mpif90 ./configure --with-zlib --enable-fortran --enable-parallel --enable-shared --prefix=<install_path>
$ make -j 8
$ make install
```

- NVHPC + HPCX:
```bash
$ CPP=cpp CFLAGS="-fPIC -m64" FCFLAGS="-fPIC -m64" CC=mpicc FC=mpif90 ./configure --with-zlib --enable-fortran --enable-parallel --enable-shared --prefix=<install_path>
$ make -j 8
$ make install
```

In either case, avoid using a root path as the install path. Adjust the number of processes in the `make -j` command according to the number of cores available in the system.

After installing set the following environment variables:

```bash
export HDF5_HOME=<install_path>
export HDF5_ROOT=$HDF5_HOME
export HDF5_DIR=$HDF5_HOME
export PATH=$HDF5_HOME/bin:$PATH
export LD_LIBRARY_PATH=$HDF5_HOME/lib:$LD_LIBRARY_PATH
```

### CGNS

- Version: Github `develop` branch

For either compiler, the CMake build process is recommended. Ensure that the system environment is set, then:

```bash
$ CC=mpicc CXX=mpicxx cmake -DBUILD_TESTING=ON -DCGNS_BUILD_CGNSTOOLS=ON -DCGNS_BUILD_SHARED=ON -DCGNS_BUILD_TESTING=ON -DCGNS_ENABLE_64BIT=ON -DCGNS_ENABLE_FORTRAN=OFF -DCGNS_ENABLE_BASE_SCOPE=ON -DCGNS_ENABLE_HDF5=ON -CGNS_ENABLE_MEM_DEBUG=ON -DCGNS_ENABLE_PARALLEL=ON -DCGNS_ENABLE_SCOPING=ON -DCGNS_ENABLE_TESTS=ON -DCGNS_USE_SHARED=ON -DCMAKE_INSTALL_PREFIX=<install_path> -DHDF5_NEED_MPI=ON -DHDF5_NEED_ZLIB=ON
$ make -j 8
$ ctest
$ make install
```

Given that CGNS depends on HDF5, ensure that the HDF5 environment variables are set before building CGNS. Before running the install command, ensure that all tests pass.

After installing set the following environment variables:

```bash
export CGNS_HOME=<install_path>
export CGNS_ROOT=$CGNS_HOME
export CGNS_DIR=$CGNS_HOME
export PATH=$CGNS_HOME/bin:$PATH
export LD_LIBRARY_PATH=$CGNS_HOME/lib:$LD_LIBRARY_PATH
export LIBRARY_PATH=$CGNS_HOME/lib:$LIBRARY_PATH
export C_INCLUDE_PATH=$CGNS_HOME/include:$C_INCLUDE_PATH
export CPLUS_INCLUDE_PATH=$CGNS_HOME/include:$CPLUS_INCLUDE_PATH
```

### CGNS2H5

With both dependencies installed, the code can be built. The code is built using CMake with the following commands:

```bash
$ cmake ..
$ make -j 8
```

If compiling on MN5 ACC, add `-DUSE_MN5=ON` to the CMake command. Also ensure that the HDF5 and CGNS environment variables are set before building.