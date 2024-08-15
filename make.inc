######################################################
#
# make.inc file incluing compilation flags
#
######################################################
# section for Fortran compiler #

# MPI fortran #
F90 = /home/leon/HDF5-GPU/bin/h5pfc 

# MPI fortran #
#F90 = mpifort
#PARALLEL = -DMPI

######################################################
# Section for GPU #

# Choose between CPU or GPU threading #
#MULTICORE = -DGPU -acc

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
F90FLAGS= -cpp -fast -O4 -Mconcur -Munroll -Mdalign -Mframe -Mlre -Mstack_arrays -Mvect=short -Mvect=prefetch -Mipa=fast,inline -Mvect=levels:4 -Mvect=fuse -Mvect=gather -Mvect=idiom -Mvect=altcode -Mvect=simd $(PARALLEL) $(MULTICORE) $(MATHPACK) $(COMPILE) 

LDFLAGS= -cpp -fast -O4 -Mconcur -Munroll -Mdalign -Mframe -Mlre -Mstack_arrays -Mvect=short -Mvect=prefetch -Mipa=fast,inline -Mvect=levels:4 -Mvect=fuse -Mvect=gather -Mvect=idiom -Mvect=altcode -Mvect=simd $(PARALLEL) $(MULTICORE) $(MATHPACK) $(COMPILE)   
          
######################################################
