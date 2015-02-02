include ../makeflags.inc

LIBS = -lm -ldl -lstdc++
# DBG = -g -gstabs+

SYMLINKS_HCVIS = gridcache.C gridcache.H variables.C variables.H \
constants.H

default:
	$(MAKE) links
	$(MAKE) iontracer

variables.o: variables.C variables.H $(hclibs)
	$(CXX) $(DBG) -c $(CXXFLAGS)  variables.C

gridcache.o : gridcache.C gridcache.H $(hclibs)
	$(CXX) $(DBG) -c $(CXXFLAGS) gridcache.C 

Token.o: Token.cpp Token.h
	$(CXX) $(DBG) -c $(CXXFLAGS) Token.cpp

TrackWriter.o: TrackWriter.cpp TrackWriter.h
	$(CXX) $(DBG) -c $(CXXFLAGS) TrackWriter.cpp

PointReader.o: PointReader.cpp PointReader.h Token.h $(hclibs)
	$(CXX) $(DBG) -c $(CXXFLAGS) PointReader.cpp 

Config.o: Config.cpp Config.h Token.h
	$(CXX) $(DBG) -c $(CXXFLAGS) Config.cpp 

Track.o: Track.cpp Track.h
	$(CXX) $(DBG) -c $(CXXFLAGS) Track.cpp 

Ptracer.o: Ptracer.cpp Ptracer.h Track.h Config.h variables.H $(hclibs) 
	$(CXX) $(DBG) -c $(CXXFLAGS) Ptracer.cpp 

iontracer.o: iontracer.cpp Track.h Ptracer.h
	$(CXX) $(DBG) -c $(CXXFLAGS) iontracer.cpp 

iontracer: Track.o Ptracer.o iontracer.o gridcache.o variables.o Config.o PointReader.o TrackWriter.o Token.o
	$(CXX) $(DBG) $(LDFLAGS) -o$@ iontracer.o Track.o Ptracer.o \
		Config.o variables.o gridcache.o PointReader.o TrackWriter.o Token.o \
		$(hclibs) $(LIBS)

links:
	for i in $(SYMLINKS_HCVIS); do ln -sf ../hcvis/$$i; done;

install:
	ln -sf ../tracer/iontracer ../bin/

clean:
	-rm -f *.o iontracer  *~ iontracer $(SYMLINKS_HCVIS)

