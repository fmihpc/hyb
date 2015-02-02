SHELL = /bin/sh

default: libhc.a

include ../makeflags.inc

include makedep.dat

HC_OBJS = realptr.o cache.o cellcache.o shmemintf.o vec.o sharr.o \
cshmat.o byteconv.o fileheader.o BC.o maps.o mempool.o gengrid.o \
grid.o cartgrid.o HCgrid.o metagrid.o 

# It is best to remove the library before adding things to it with ar cr.
# In some cases it is not properly updated otherwise.

libhc.a : $(HC_OBJS)
	-rm -f libhc.a
	$(AR$) cr libhc.a $(HC_OBJS)

../lib/libhc.a : libhc.a
	mkdir -p ../lib/; cp libhc.a ../lib/; chmod og+r ../lib/libhc.a

install: ../lib/libhc.a
	mkdir -p ../include/
	cp *.H ../include
	chmod og+r ../include/*.H

dep:
	$(CXX) -E -MM `echo $(HC_OBJS) | sed -e 's/\.o/.C/g' | tr -d '\012'` >makedep.dat

clean:
	-rm -f libhc.a *.o *.so *~
