######################################################
#
# make.inc file incluing compilation flags
#
######################################################
# section for Fortran compiler #

# MPI fortran #
#F90=mpifort
#PARALLEL = -DMPI

# The usual gfortran
F90=h5fc #gfortran

######################################################
# section for math packages #

# need extra math packages? #
#MATHPACK = -llapack -lblas

######################################################
# Section for GPU #

# Choose between CPU or GPU threading #
#MULTICORE = -DGPU -acc

######################################################
# Section for debug #

# Want to look for compilation message? #
#COMPILE = -Minfo=accel

######################################################

# Define F90 flags #
F90FLAGS= -cpp -Ofast -O3 -march=native -mtune=native -flto -funroll-loops \
          -ffp-contract=fast -ftree-vectorize -fprefetch-loop-arrays \
	        -ffree-line-length-none $(MULTICORE) $(COMPILE) $(PARALLEL) $(MATHPACK)

LDFLAGS= -cpp -Ofast -O3 -march=native -mtune=native -flto -funroll-loops \
         -ffp-contract=fast -ftree-vectorize -fprefetch-loop-arrays \
         -ffree-line-length-none $(MULTICORE) $(COMPILE) $(PARALLEL)$(MATHPACK)

######################################################
