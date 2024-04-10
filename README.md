# CGNS2H5

Small application to convert CGNS files to HDF5 files in the format of SOD2D.

## Dependencies

- [NVIDIA HPC SDK](https://developer.nvidia.com/hpc-sdk)
- [CGNS](https://cgns.github.io)
- [HDF5](https://www.hdfgroup.org/solutions/hdf5)

The HDF5 must be compiled in parallel using the MPI compile wrapper. CGGNS must be compiled with the HDF5 support, also in parallel.