OPTFLAGS = -O -fpermissive -Wall
LDFLAGS =
SHARED = -shared
GLINCLUDE =
GLUTLIB = -lglut
GLLIBS = -lGLU -lGL
#GLLIBS = -lMesaGLU -lMesaGL
#GLLIBS = -lMesaGLU /usr/X11R6/lib/libGL.so
XINCLUDE = -I/usr/X11R6/include -I/usr/include/tcl8.5 -I/usr/include/tcl8.5/tk-private/generic
XLIBS = -L/usr/X11R6/lib -lXmu -lXt -lSM -lICE -lXext -lXi -lX11
CXXFLAGS = $(OPTFLAGS) $(DEFS) -I. -I../include

.C.o:
	$(CXX) -c $(CXXFLAGS) $<

.C.s:
	$(CXX) -S $(CXXFLAGS) $<

hclibs = ../lib/libhc.a
