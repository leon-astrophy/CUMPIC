CUMPIC - MPI Parallelized version of CUMC3D 

Pre-requisite - HDF5 1.12, Nvidia HPC SDK 23.9

Nvidia HPC SDK - https://developer.nvidia.com/nvidia-hpc-sdk-239-downloads

HDF5 - [https://portal.hdfgroup.org/downloads/hdf5/hdf5_1_14_4.html](https://portal.hdfgroup.org/downloads/hdf5/hdf5_1_12_3.html)

HDF5 Installation:

FCFLAGS="-fPIC -m64 -tp=px" FC="mpifort installed in your Nvidia HPC SDK directory" ./configure --prefix="Your HDF5 directory" --enable-fortran --enable-parallel

make

make install
