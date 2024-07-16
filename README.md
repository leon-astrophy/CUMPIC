!***************************************************

! CUMPIC - MPI Parallelized version of CUMC3D !

!***************************************************

Pre-requisite - OPENMPI 5.0, HDF5 1.14, and GCC-11

!******************************************************************!

OPENMPI - https://www.open-mpi.org/software/ompi/v5.0/

HDF5 - https://portal.hdfgroup.org/downloads/hdf5/hdf5_1_14_4.html

!******************************************************************!

OPENMPI Installation:

FC="GCC-11's gfortran" ./configure --prefix="Your OPENMPI directory" 

make -j4 all

make install

!******************************************************************!

HDF5 Installation:

FC="mpifort installed in your OPENMPI directory" ./configure --prefix="Your HDF5 directory" --enable-fortran --enable-parallel

make

make install

!******************************************************************!

