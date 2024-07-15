######################################################
#
# make.inc file incluing compilation flags
#
######################################################
# section for Fortran compiler #

# MPI fortran #
F90 = /home/leon/HDF5-MPIFORT/bin/h5pfc
PARALLEL = -DMPI

# The usual gfortran
#F90 = gfortran

######################################################
# Section for GPU #

# Choose between CPU or GPU threading #
#MULTICORE = -DGPU -fopenacc

######################################################
# section for math packages #

# need extra math packages? #
#MATHPACK = -llapack -lblas

######################################################
# Section for debug #

# Want to look for compilation message? #
#COMPILE = -Minfo=accel

######################################################
# these are for standard gcc

# Define F90 flags #
F90FLAGS= -cpp -Ofast -O3 -march=native -mtune=native -flto -funroll-loops \
          -ffp-contract=fast -ftree-vectorize -fprefetch-loop-arrays \
	        -ffree-line-length-none $(PARALLEL) $(HDF5) $(MULTICORE) $(COMPILE) $(MATHPACK)

LDFLAGS= -cpp -Ofast -O3 -march=native -mtune=native -flto -funroll-loops \
         -ffp-contract=fast -ftree-vectorize -fprefetch-loop-arrays \
         -ffree-line-length-none $(PARALLEL) $(HDF5) $(MULTICORE) $(COMPILE) $(MATHPACK)     
          
######################################################
