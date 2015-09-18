default: Kriging

#Files
objects := FastInverse.o ByteSwapping.o IO.o dgesvd.o neighbor3.o	\
fwtd3.o vande3.o matrixio.o hb.o linearSolver.o second.o

matlabobjects := FastInverse.o ByteSwapping.o IO.o dgesvd.o		\
neighbor3.o fwtd3.o vande3.o matrixio.o hb.o matlablinearSolver.o	\
second.o

CC := mpicc
CXX := mpicxx

  VPATH := ./src/3d   # 3D version
# VPATH := ./src/2d   # 2D version

#DEBUG_MODE := yes
ifdef DEBUG_MODE
CPPFLAGS :=  -g -O0
else
CPPFLAGS := -O3 -DNDEBUG -march=native -m64 -flto
endif

#PIC := yes
ifdef PIC
CPPFLAGS := -fPIC $(CPPFLAGS)
endif


# Put PETSC directory here
PETSC_DIR := /home/julio/petsc-3.3-p5
# PETSC_ARCH := arch-linux2-cxx-opt
PETSC_ARCH := arch-linux2-c-opt

include ${PETSC_DIR}/conf/variables
CPPFLAGS := $(CPPFLAGS) -Wall -I${PETSC_DIR}/include		\
-I${PETSC_DIR}/${PETSC_ARCH}/include -I../Modifiednew3dMKL	\
-I../CXSparse/Include -I../SuiteSparse_config/ #-DCS_LONG

LDFLAGS := $(LDFLAGS) -L../CXSparse/Lib/ # -L../amdlibm/lib/dynamic 
LDFLAGS := $(LDFLAGS) -L${PETSC_DIR}/${PETSC_ARCH}/lib
LDFLAGS := $(LDFLAGS) -L$(MKLROOT)/../compiler/lib/intel64

LDLIBS :=  -lm -liomp5 -lpthread -ldl  ${PETSC_LIB} -lgsl
#STATLIBS := ../rio/exafmm-dev/wrappers/libmatern.a


STATLIBS := $(STATLIBS) ../Modifiednew3dMKL/libnew3d.a	\
../CXSparse/Lib/libcxsparse.a -Wl,--start-group		\
$(MKLROOT)/lib/intel64/libmkl_gf_lp64.a			\
$(MKLROOT)/lib/intel64/libmkl_intel_thread.a		\
$(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group
#/usr/local/lib/libfftw3.a 

# kriginglibrary: $(objects) ar -rv kriginglibrary.a $(objects)
#			kriginglibrary.a
Kriging : mainkriging.o $(objects)
	$(CXX) $(CPPFLAGS) $(LDFLAGS) mainkriging.o $(objects) $(STATLIBS) $(LDLIBS)  -o Kriging

.PHONY:libkriging
libkriging: libkriging.so
libkriging.so: libkriging.o $(matlabobjects)
	$(CXX) -shared $(CPPFLAGS) $(LDFLAGS) libkriging.o $(matlabobjects) $(STATLIBS) $(LDLIBS)  -o $@

FastInverse.o: FastInverse.h
ByteSwapping.o: ByteSwapping.h
IO.o: IO.h
matrixio.o: matrixio.h

dgesvd.o : dgesvd.h
neigbor3.o : neighbor3.h
fwtd3.o : fwtd3.h
vande3.o : vande3.h
hb.o : hb.h
linearSolver.o: linearSolver.h
matlablinearSolver.o: matlablinearSolver.h

clean:
	rm *.o 	Kriging *.dat

