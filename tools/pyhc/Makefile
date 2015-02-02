include ../makeflags.inc

LIBS = -lm -ldl 
# DBG = -g -gstabs+
CXXFLAGS = -fpic

SYMLINKS_SRC = BC.C BC.H byteconv.C byteconv.H cache.C cache.H \
cartgrid.C cartgrid.H cellcache.C cellcache.H cshmat.C cshmat.H \
fileheader.C fileheader.H gengrid.C gengrid.H grid.C grid.H HCgrid.C \
HCgrid.H maps.C maps.H mempool.C mempool.H metagrid.C metagrid.H \
realptr.C realptr.H sharr.C sharr.H shmemintf.C shmemintf.H state.H \
vec.C vec.H

SYMLINKS_HCVIS = gridcache.C gridcache.H variables.C variables.H \
constants.H

HC_OBJ =  realptr.o cache.o cellcache.o shmemintf.o vec.o sharr.o \
cshmat.o byteconv.o fileheader.o BC.o maps.o mempool.o gengrid.o \
grid.o cartgrid.o HCgrid.o metagrid.o gridcache.o

default:
	$(MAKE) links
	$(MAKE) _pyhc.so

variables.o: variables.C variables.H 
	$(CXX) $(DBG) -c $(CXXFLAGS) variables.C

pyhc.o: pyhc.cpp pyhc.h
	$(CXX) $(DBG) -c $(CXXFLAGS) pyhc.cpp

pyhc_wrap.cxx: pyhc.cpp pyhc.h 
	swig -c++ -python pyhc.i

pyhc_wrap.o: pyhc_wrap.cxx pyhc.cpp pyhc.h 
	$(CXX) $(DBG) -c $(CXXFLAGS) pyhc_wrap.cxx -I/usr/include/python2.7

_pyhc.so: pyhc_wrap.o pyhc.o  variables.o $(HC_OBJ) 
	$(CXX) -shared pyhc_wrap.o pyhc.o variables.o $(HC_OBJ) -o _pyhc.so 

links:
	for i in $(SYMLINKS_SRC); do ln -sf ../src/$$i; done;
	for i in $(SYMLINKS_HCVIS); do ln -sf ../hcvis/$$i; done;

clean:
	-rm -f *.o pyhc.py pyhc.pyc pyhc_wrap.cxx _pyhc.so *~ $(SYMLINKS_SRC) $(SYMLINKS_HCVIS)
