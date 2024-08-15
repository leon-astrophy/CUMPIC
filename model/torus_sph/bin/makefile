#include header#
include ../../../make.inc

#specify HDF5#
EXTRAINCS += $(HDF5INCS)
EXTRAOBJECTS += $(HDF5LIBS)

#specify directories#
SOLVER := ../../../solver
SOURCE := ../src

#source list#
SOURCE1 = $(shell cat $(SOLVER)/source.list)
SRC1 = $(notdir $(SOURCE1))
OBJ1 = $(patsubst %.f90,%.o,$(SRC1))

#local source files#
SOURCE2 = $(shell cat $(SOURCE)/source.list) #$(wildcard *.f90 $(SOURCE) *.f90)
SRC2 = $(notdir $(SOURCE2))
OBJ2 = $(patsubst %.f90,%.o,$(SRC2))

ALL: $(OBJ1) $(OBJ2)
	echo $(SRC1)
	echo $(SRC2)
	$(F90) $(LDFLAGS) -o ../CUMC3D $(OBJ1) $(OBJ2) $(EXTRAOBJECTS)

$(OBJ1): %.o: $(SOLVER)/%.f90 
	$(F90) $(F90FLAGS) $(EXTRAINCS) -c $< -o $@

$(OBJ2): %.o: $(SOURCE)/%.f90 
	$(F90) $(F90FLAGS) $(EXTRAINCS) -c $< -o $@
