# Copyright (C) 1997-2000 Gidon Moont

# Biomolecular Modelling Laboratory
# Imperial Cancer Research Fund
# 44 Lincoln's Inn Fields
# London WC2A 3PX

# +44 (0)20 7269 3565
# http://www.bmm.icnet.uk/

#############

# This line you will definitely have to edit

CUDA_PATH := /usr/local/cuda

#############

# You may need/want to edit some of these
#
# Hint: For the NVCC_FLAGS have a look at what the fftw build used

SHELL           = /bin/sh

NVCC              = nvcc

NVCC_FLAGS        = -arch=sm_75 -rdc=true -m64

NVCC_LINKERS      = -lm

STRIP           = strip

SECURITY	= chmod +x

LIBS = -lcufft -lcuda -lcudart -lcurand

LIBPATH := $(CUDA_PATH)/lib64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda/lib64
SOLIBS := $(LIBPATH)/libcufft.so $(LIBPATH)/libcudart.so $(LIBPATH)/libcurand.so 

####################################################

# You should not be editing anything below here

NVCC_FLAGS_FULL	= $(NVCC_FLAGS)
FFTW_LINKERS    = -L$(FFTW_DIR)

#############

.SUFFIXES:	.cu .o

.cu.o:
		$(NVCC) $(NVCC_FLAGS_FULL) -c $<

#############

LIBRARY_OBJECTS =manipulate_structures.o angles.o coordinates.o electrostatics.o grid.o qsort_scores.o

PROGRAMS = ftdock build randomspin

all:		$(PROGRAMS)

#############

ftdock:	ftdock.o $(LIBRARY_OBJECTS) structures.cuh
		$(NVCC) $(NVCC_FLAGS_FULL) ftdock.o $(LIBRARY_OBJECTS) $(LIBS) -o $@ 
		$(STRIP) $@
		$(SECURITY) $@
		
#############

build:		build.o $(LIBRARY_OBJECTS) structures.cuh
		$(NVCC) $(NVCC_FLAGS) $(LIBRARY_OBJECTS) $(LIBS) -o $@ build.o 

#############

randomspin:	randomspin.o $(LIBRARY_OBJECTS) structures.cuh
		$(NVCC) $(NVCC_FLAGS) $(LIBRARY_OBJECTS) $(LIBS) -o $@ randomspin.o

#############

clean:
		rm -f *.o core $(PROGRAMS)

#############

# dependencies

ftdock.o:			structures.cuh
build.o:			structures.cuh
randomspin.o:			structures.cuh

angles.o:			structures.cuh
coordinates.o:			structures.cuh
electrostatics.o:		structures.cuh
grid.o:				structures.cuh
manipulate_structures.o:	structures.cuh
qsort_scores.o:			structures.cuh
