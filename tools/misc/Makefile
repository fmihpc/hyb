PROGS = hc2ppm hctrans hc2tecplot

SHELL = /bin/sh

default: $(PROGS)

include ../makeflags.inc

hc2ppm.o : hc2ppm.C $(hclibs)

hctrans.o : hctrans.C $(hclibs)

hc2ppm : hc2ppm.o $(hclibs)
	$(CXX) $(LDFLAGS) $(DEFS) -o $@ hc2ppm.o $(hclibs) -lm $(LIBFPE)

hc2tecplot : hc2tecplot.o $(hclibs)
	$(CXX) $(LDFLAGS) $(DEFS) -o $@ hc2tecplot.o $(hclibs) -lm $(LIBFPE)

hctrans : hctrans.o $(hclibs)
	$(CXX) $(LDFLAGS) $(DEFS) -o $@ hctrans.o $(hclibs) -lm $(LIBFPE)

install: $(PROGS)
	mkdir -p ../bin/
	cp $(PROGS) *.sh ../bin
	chmod u+rx ../bin/*

clean:
	-rm -f *.o *.so $(PROGS) *~
