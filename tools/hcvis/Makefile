include ../makeflags.inc

# If you have the BLT library (libBLT.a) installed, uncomment the
# following two lines. If not, leave them commented. Everything works
# fine without BLT, except that the "1D view window" menu item under
# File menu does nothing. 
HAVE_BLT = false

ifeq ($(HAVE_BLT),true)
BLTLIB = -lBLT
BLTDEFS = -DHAVE_BLT=1
endif

LIBS_HCVIS = -L/usr/X11R6/lib $(BLTLIB) -ltk8.5 -ltcl8.5 $(GLLIBS) $(XLIBS) -lm -ldl

default: hcvis hcintpol

all: hcvis hcintpol hwa-hcintpol hc2vtk hc2gridxyz

togl.o: togl.c togl.h
	$(CC) -c $(CXXFLAGS) $(XINCLUDE) togl.c

hcvis.o: hcvis.C togl.h toglwin.H palette.H contour.H variables.H 3Dobj.H $(hclibs)
	$(CXX) -c $(CXXFLAGS) $(XINCLUDE) $(BLTDEFS) hcvis.C

toglwin.o: toglwin.C toglwin.H contour.H palette.H gridcache.H intpolcache.H GLaxis.H variables.H 3Dobj.H $(hclibs)
	$(CXX) -c $(CXXFLAGS) $(XINCLUDE) toglwin.C

GLaxis.o : GLaxis.C GLaxis.H
	$(CXX) -c $(CXXFLAGS) $(XINCLUDE) $(GLINCLUDE) GLaxis.C

zoomstack.o : zoomstack.C
	$(CXX) -c $(CXXFLAGS) $(XINCLUDE) $(GLINCLUDE) zoomstack.C

hcintpol.o : hcintpol.C variables.H gridcache.H $(hclibs)
hwa-hcintpol.o : hwa-hcintpol.C variables.H gridcache.H toglwin.H $(hclibs)
palette.o : palette.C palette.H
gridcache.o : gridcache.C gridcache.H $(hclibs)
intpolcache.o : intpolcache.C intpolcache.H
contour.o : contour.C contour.H
variables.o : variables.C variables.H $(hclibs)
3Dobj.o : 3Dobj.C 3Dobj.H

hcvis: togl.o hcvis.o gridcache.o intpolcache.o toglwin.o palette.o contour.o zoomstack.o GLaxis.o variables.o 3Dobj.o
	$(CXX) $(LDFLAGS) -o $@ hcvis.o togl.o toglwin.o contour.o \
		palette.o zoomstack.o GLaxis.o gridcache.o intpolcache.o variables.o \
		3Dobj.o \
		$(hclibs) $(LIBS_HCVIS)

hcintpol: hcintpol.o gridcache.o variables.o
	$(CXX) $(LDFLAGS) -o $@ hcintpol.o gridcache.o variables.o $(hclibs)

hwa-hcintpol: hwa-hcintpol.o gridcache.o variables.o
	$(CXX)  $(LDFLAGS) -o $@ hwa-hcintpol.o gridcache.o variables.o $(hclibs)

hc2vtk: hc2vtk.o gridcache.o variables.o
	$(CXX)  $(LDFLAGS) -o $@ hc2vtk.o gridcache.o variables.o $(hclibs)

hc2gridxyz: hc2gridxyz.o gridcache.o variables.o
	$(CXX)  $(LDFLAGS) -o $@ hc2gridxyz.o gridcache.o variables.o $(hclibs)

install: hcvis hcintpol
	mkdir -p ../bin/
	ln -sf ../hcvis/hcvis ../bin/
	ln -sf ../hcvis/hcintpol ../bin
	ln -sf ../hcvis/hcv ../bin
	chmod og+r *.tcl
	mkdir -p ../lib/
	cp hcvis.tcl hcvisfilesel.tcl hcvis.help hcvisobj.tcl ../lib/
	chmod og+rx hcvis hcintpol .
	chmod og+r ../lib/hcvis.tcl ../lib/hcvisfilesel.tcl ../lib/hcvisobj.tcl ../lib/hcvis.help
	if [ ! -d ../lib/palettes ]; then mkdir ../lib/palettes; cp palettes/* ../lib/palettes/; chmod og+r ../lib/palettes/*; fi
	chmod og+rx ../lib/palettes
	chmod og+rx palettes
	chmod -R og+r palettes
	chmod u+x hcv

clean:
	-rm -f *.o *~ hcvis hcintpol hwa-hcintpol hc2vtk hc2gridxyz
