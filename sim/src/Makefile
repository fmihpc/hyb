# Compile time options
USE_SPHERICAL_COORDINATE_SYSTEM := false
USE_PARTICLE_SUBCYCLING := false
IGNORE_ELECTRIC_FIELD_HALL_TERM := false
PERIODIC_FIELDS_Y := false
RECONNECTION_GEOMETRY := false
VTK_SHOW_GHOST_CELLS := false
NO_DIAGNOSTICS := false
SAVE_POPULATION_AVERAGES := false
SAVE_PARTICLES_ALONG_ORBIT := false
SAVE_PARTICLE_CELL_SPECTRA := false

SHELL = /bin/bash

default : HYB

# Executable name
PROGRAM_NAME := hyb

# Compile information (data_user_host)
COMPILE_INFO := $(shell date --iso-8601=seconds)_$(shell whoami)_$(shell hostname)

# Remove spaces
COMPILE_INFO := \"$(shell echo $(COMPILE_INFO) | sed 's/ //g')\"

# General settings
CXX_GEN_OPTS := -Wall -ansi -pedantic
LINKINGOPTIONS := -lm -static

ifeq ($(USE_SPHERICAL_COORDINATE_SYSTEM),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DUSE_SPHERICAL_COORDINATE_SYSTEM -Wno-unused-variable
endif

ifeq ($(USE_PARTICLE_SUBCYCLING),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DUSE_PARTICLE_SUBCYCLING
endif

ifeq ($(IGNORE_ELECTRIC_FIELD_HALL_TERM),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DIGNORE_ELECTRIC_FIELD_HALL_TERM
endif

ifeq ($(PERIODIC_FIELDS_Y),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DPERIODIC_FIELDS_Y
endif

ifeq ($(RECONNECTION_GEOMETRY),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DRECONNECTION_GEOMETRY
endif

ifeq ($(VTK_SHOW_GHOST_CELLS),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DVTK_SHOW_GHOST_CELLS
endif

ifeq ($(NO_DIAGNOSTICS),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DNO_DIAGNOSTICS
endif

ifeq ($(SAVE_POPULATION_AVERAGES),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DSAVE_POPULATION_AVERAGES
endif

ifeq ($(SAVE_PARTICLES_ALONG_ORBIT),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DSAVE_PARTICLES_ALONG_ORBIT
endif

ifeq ($(SAVE_PARTICLE_CELL_SPECTRA),true)
CXX_GEN_OPTS := $(CXX_GEN_OPTS) -DSAVE_PARTICLE_CELL_SPECTRA
endif

# Compiler settings - default
HYB : CXX = g++
HYB : CXXFLAGS = -O2 -fomit-frame-pointer -ffast-math -pipe -fno-aggressive-loop-optimizations

# Compiler settings - debug
debug : CXX = g++
debug : CXXFLAGS = -g

# All program object files
OBJECTS = \
atmosphere.o backgroundcharge.o boundaries.o chemistry.o definitions.o \
detector.o diagnostics.o forbidsplitjoin.o grid.o logger.o \
magneticfield.o main.o params.o particle.o population_exospheric.o \
population_imf.o population_ionospheric.o population.o \
population_solarwind.o population_uniform.o random.o refinement.o \
resistivity.o simulation.o splitjoin.o timepool.o vectors.o \
vis_data_source_simulation.o vis_db_vtk.o

# Create and include Makefile dependencies
Makefile.deps :
	$(CXX) $(CXXFLAGS) $(CXX_GEN_OPTS) -MM *.cpp vis/*.cpp >Makefile.deps

-include Makefile.deps

# Object specific compilation instructions and dependencies

atmosphere.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) atmosphere.cpp
backgroundcharge.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) backgroundcharge.cpp 
boundaries.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) boundaries.cpp
chemistry.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) chemistry.cpp
definitions.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) definitions.cpp
detector.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) detector.cpp
diagnostics.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) diagnostics.cpp
forbidsplitjoin.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) forbidsplitjoin.cpp 
grid.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) -funroll-loops grid.cpp
logger.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) logger.cpp
magneticfield.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) magneticfield.cpp
main.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) main.cpp
params.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) params.cpp -DCOMPILE_INFO=$(COMPILE_INFO)
particle.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) particle.cpp
population_exospheric.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population_exospheric.cpp
population_imf.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population_imf.cpp
population_ionospheric.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population_ionospheric.cpp
population.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population.cpp
population_solarwind.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population_solarwind.cpp
population_uniform.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) population_uniform.cpp
random.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) random.cpp
refinement.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) refinement.cpp 
resistivity.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) resistivity.cpp
simulation.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) simulation.cpp
splitjoin.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) splitjoin.cpp
timepool.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) timepool.cpp
vectors.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) vectors.cpp
vis_data_source_simulation.o :
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) vis/vis_data_source_simulation.cpp
vis_db_vtk.o:
	$(CXX) -c $(CXXFLAGS) $(CXX_GEN_OPTS) vis/vis_db_vtk.cpp

# Main targets
HYB : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXX_GEN_OPTS) -o $(PROGRAM_NAME) $^ $(LINKINGOPTIONS)

debug : $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(CXX_GEN_OPTS) -o $(PROGRAM_NAME) $^ $(LINKINGOPTIONS)

doc :
	rm -fr doc/;
	doxygen Doxyfile;
	cd doc/; ln -s html/index.html ./;
	cd doc/latex/; $(MAKE); mv refman.pdf ../; cd ..; rm -fr latex;

clean:
	rm -f hyb Makefile.deps *.o *.hc *.vtk *.dat *.log *.err *~ */*~ vis/*.o
	rm -fr doc/

