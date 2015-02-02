SHELL = /bin/sh
DIRS = src hcvis misc pyhc tracer

default: libs

libs:
	mkdir -p bin/ include/ lib/
	cd src; $(MAKE) install

all:
	$(MAKE) libs
	cd hcvis; $(MAKE) all; $(MAKE) install; cd ..;
	cd misc; $(MAKE); $(MAKE) install; cd ..;
	cd pyhc; $(MAKE); cd ..;
	cd tracer; $(MAKE); $(MAKE) install; cd ..;

hcintpol:
	$(MAKE) libs
	cd hcvis; $(MAKE) hcintpol; $(MAKE) hwa-hcintpol; cd ..;

clean:
	for dir in $(DIRS); do cd $${dir}; $(MAKE) clean; cd ..; done
	-rm -rf lib/ bin/ include/
	-rm -f *~ */*~ */*/*~
