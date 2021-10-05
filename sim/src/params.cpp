/** This file is part of the HYB simulation platform.
 *
 *  Copyright 2014- Finnish Meteorological Institute
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <sstream>
#include "params.h"
#include "simulation.h"
#include "chemistry.h"

using namespace std;

const string Params::codeVersion = "HYB simulation platform ("
#ifdef COMPILE_INFO
                                   COMPILE_INFO
#else
                                   __DATE__ " "  __TIME__
#endif
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
                                   " USE_SPHERICAL_COORDINATE_SYSTEM"
#endif
#ifdef USE_PARTICLE_SUBCYCLING
                                   " USE_PARTICLE_SUBCYCLING"
#endif
#ifdef IGNORE_ELECTRIC_FIELD_HALL_TERM
                                   " IGNORE_ELECTRIC_FIELD_HALL_TERM"
#endif
#ifdef PERIODIC_FIELDS_Y
                                   " PERIODIC_FIELDS_Y"
#endif
#ifdef RECONNECTION_GEOMETRY
                                   " RECONNECTION_GEOMETRY"
#endif
#ifdef VTK_SHOW_GHOST_CELLS
                                   " VTK_SHOW_GHOST_CELLS"
#endif
#ifdef NO_DIAGNOSTICS
                                   " NO_DIAGNOSTICS"
#endif
#ifdef SAVE_POPULATION_AVERAGES
                                   " SAVE_POPULATION_AVERAGES"
#endif
#ifdef SAVE_PARTICLES_ALONG_ORBIT
                                   " SAVE_PARTICLES_ALONG_ORBIT"
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
                                   " SAVE_PARTICLE_CELL_SPECTRA"
#endif
                                   ")";

//! Indicates whether the program execution is in initial phase
bool Params::initPhase = true;

//! Indicated whether the class is reading population config information
bool Params::readingPopulation = false;

//! Indicated whether the class is reading population config information
bool Params::readingDetector = false;

//! Indicated whether the class is reading process config information
bool Params::readingProcess = false;

//! Indicates whether the program execution is to be terminated
bool Params::stoppingPhase = false;

static const real Pi = M_PI;

//! Config file name [-]
const char* Params::configFileName = "hybrid.cfg";

//! Indicates whether whole state file is given
bool Params::wsFileGiven = false;

//! Name of the whole state file
string Params::wsFileName;

// =================================================================================
// =============================== PHYSICAL CONSTANTS ==============================
// =================================================================================

//! Speed of light in vacuum [m/s]
const real Params::c = 299792458.0;

//! Electron charge magnitude [C]
const real Params::e = 1.60217653e-19;

//! Boltzmann constant [J/K]
const real Params::k_B = 1.3806505e-23;

//! Permeability of free space [N A^-2]
const real Params::mu_0 = 4.0*pi*1e-7;

//! Permittivity of free space [F/m]
const real Params::eps_0 = 1.0/(Params::mu_0*Params::c*Params::c);

//! Gravitational constant [m^3 kg^-1 s^-2]
const real Params::G = 6.6742e-11;

//! Astronomical unit [m]
const real Params::AU = 149597870660.0;

//! Adiabatic constant []
const real Params::gamma = 5.0/3.0;

//! Atomic mass unit [kg]
const real Params::amu = 1.66053886e-27;

//! Electron mass [kg]
const real Params::m_e = 548.57990945e-6*Params::amu;

//! Proton mass [kg]
const real Params::m_p = 1.00727646688*Params::amu;

//! Atomic hydrogen mass [kg]
const real Params::m_H = 1.00794*Params::amu;

//! Atomic hydrogen molecule mass [kg]
const real Params::m_H2 = 2.0158*Params::amu;

//! Atomic deuterium (hydrogen-2) mass [kg]
const real Params::m_2H = 2.01355321270*Params::amu;

//! Atomic tritium (hydrogen-3) mass [kg]
const real Params::m_3H = 3.0160492*Params::amu;

//! Atomic helium-3 mass [kg]
const real Params::m_3He = 3.0160293*Params::amu;

//! Atomic helium-4 mass [kg]
const real Params::m_He = 4.002602*Params::amu;

//! Atomic nitrogen mass [kg]
const real Params::m_N = 14.00674*Params::amu;

//! Molecular nitrogen mass [kg]
const real Params::m_N2 = 28.0134*Params::amu;

//! Atomic oxygen mass [kg]
const real Params::m_O = 15.9994*Params::amu;

//! Molecular oxygen mass [kg]
const real Params::m_O2 = 31.9988*Params::amu;

//! Methyl radical mass [kg]
const real Params::m_CH3 = 15.0347*Params::amu;

//! Methane mass [kg]
const real Params::m_CH4 = 16.0426*Params::amu;

//! Atomic sodium (natrium) mass [kg]
const real Params::m_Na = 22.98976928*Params::amu;

//! Ethyl radical mass [kg]
const real Params::m_C2H5 = 29.0615*Params::amu;

//! Ethane mass [kg]
const real Params::m_C2H6 = 30.0694*Params::amu;

//! Propane mass [kg]
const real Params::m_C3H8 = 44.0962*Params::amu;

//! Mercury mass [kg]
const real Params::M_Me = 0.3302e24;

//! Mercury radius [m]
const real Params::R_Me = 2439.7e3;

//! Venus mass [kg]
const real Params::M_V = 4.8685e24;

//! Venus radius [m]
const real Params::R_V = 6051.8e3;

//! Moon mass [kg]
const real Params::M_Mo = 0.07349e24;

//! Moon radius [m]
const real Params::R_Mo = 1737.1e3;

//! Mars mass [kg]
const real Params::M_Ma = 0.64185e24;

//! Mars radius [m]
const real Params::R_Ma = 3390e3;

//! Titan mass [kg]
const real Params::M_T = 0.13455e24;

//! Titan radius [m]
const real Params::R_T = 2575e3;

//! Pluto mass [kg]
const real Params::M_Pl = 1.305e22;

//! Pluto radius [m]
const real Params::R_Pl = 1150e3;

// =================================================================================
// ============================= SIMULATION PARAMETERS =============================
// =================================================================================

real Params::tempRealA = 0;
real Params::tempRealB = 0;
real Params::tempRealC = 0;
int Params::tempIntA = 0;
int Params::tempIntB = 0;
int Params::tempIntC = 0;

#ifdef USE_PARTICLE_SUBCYCLING
int Params::subcycleMaxLevel = 256;
int Params::subcycleType = 1;
real Params::dt_psub[256];
real Params::accum_psubfactor[256];
#endif

// Grid functions
GridRefinementProfile Params::gridRefinementFunction;
ResistivityProfile Params::resistivityFunction;
ForbidSplitAndJoinProfile Params::forbidSplitAndJoinFunction;
BackgroundChargeDensityProfile Params::bgChargeDensityFunction;

// Magnetic field functions
vector<MagneticFieldProfile> Params::initialMagneticFieldProfile;
vector<MagneticFieldProfile> Params::constantMagneticFieldProfile;

//! Maximum number of particle populations
const int Params::MAX_POPULATIONS = 50;

//! Number of particle populations
int Params::POPULATIONS = 0;

//! Particle populations - particle population handling is wrapped behind this variable
vector<Population*> Params::pops;

//! Creates populations
PopulationFactory Params::popFactory;

// TEMPLATE POPULATION VARIABLES
string Params::population = "";
string Params::idStr = "";
string Params::hcFilePrefix = "";
string Params::boundaryFUNC = "";
real Params::m = 0;
real Params::q = 0;
real Params::T = 0;
real Params::vth = 0;
real Params::macroParticlesPerDt = 0;
bool Params::propagateV = 1;
bool Params::accumulate = 1;
bool Params::split = 1;
bool Params::join = 1;
bool Params::logParams;
real Params::n = 0;
real Params::V = 0;
real Params::backWallWeight = 0;
real Params::R = 0;
real Params::totalRate = 0;
string Params::distFunc = "";
real subcycleSteps;


//! Used in checkIdStrAlreadyFound()
vector<string> Params::idStrTbl;

//! Check if the given idStr has already been read. If not, add it to the vector.
bool Params::checkIdStrAlreadyFound(string str)
{
    bool idStrFound = false;
    for(unsigned int i = 0; i < idStrTbl.size(); ++i) {
        if(idStrTbl[i].compare(str) == 0) {
            idStrFound = true;
            break;
        }
    }

    if(idStrFound == false) {
        idStrTbl.push_back(str);
    }
    return idStrFound;
}

//! Used in checkProcIdStrAlreadyFound()
vector<string> Params::procIdStrTbl;

//! Check if the given idStr has already been read. If not, add it to the vector.
bool Params::checkProcIdStrAlreadyFound(string str)
{
    bool procIdStrFound = false;
    for(unsigned int i = 0; i < procIdStrTbl.size(); ++i) {
        if(procIdStrTbl[i].compare(str) == 0) {
            procIdStrFound = true;
            break;
        }
    }

    if(procIdStrFound == false) {
        procIdStrTbl.push_back(str);
    }
    return procIdStrFound;
}

//! Reset population variables
void Params::clearPopulationVars()
{
    // TEMPLATE POPULATION VARIABLES
    idStr = "";
    hcFilePrefix = "";
    boundaryFUNC = "";
    m = 0;
    q = 0;
    T = 0;
    vth = 0;
    macroParticlesPerDt = 0;
    propagateV = 1;
    accumulate = 1;
    split = 1;
    join = 1;
    logParams = 0;
    n = 0;
    V = 0;
    backWallWeight = 0;
    R = 0;
    totalRate = 0;
    distFunc = "";
    subcycleSteps = 10;
}

#define GETPOPVAR(popVar) var = lookupVar(#popVar); args.popVar.value = popVar; args.popVar.given = var->updatedFromFile;

//! Create argument structure for the class Population
PopulationArgs Params::getPopulationArgsStruct()
{
    PopulationArgs args;
    struct dynamicVar *var;

    // TEMPLATE POPULATION VARIABLES
    GETPOPVAR(idStr);
    GETPOPVAR(hcFilePrefix);

    // Set population distribution boundary condition functions
    var = lookupVar("boundaryFUNC");
    if(var->updatedFromFile == true) {
        vector<string> funcName;
        vector< vector<real> > funcArgs;
        getFunctionNamesAndArgs("boundaryFUNC",funcName,funcArgs);
        if( funcName.size() == funcArgs.size() ) {
            args.boundaryFUNC.name = funcName;
            args.boundaryFUNC.funcArgs = funcArgs;
        } else {
            ERRORMSG2("error when reading population boundary conditions",idStr);
            doabort();
        }
    }
    args.boundaryFUNC.given = var->updatedFromFile;

    GETPOPVAR(m);
    GETPOPVAR(q);
    GETPOPVAR(T);
    GETPOPVAR(vth);
    GETPOPVAR(macroParticlesPerDt);
    GETPOPVAR(propagateV);
    GETPOPVAR(accumulate);
    GETPOPVAR(split);
    GETPOPVAR(join);
    GETPOPVAR(logParams);
    GETPOPVAR(n);
    GETPOPVAR(V);
    GETPOPVAR(backWallWeight);
    GETPOPVAR(R);
    GETPOPVAR(totalRate);

#ifdef USE_PARTICLE_SUBCYCLING
    GETPOPVAR(subcycleSteps);
#endif

    // Set population distribution function + arguments
    var = lookupVar("distFunc");
    if(var->updatedFromFile == true) {
        vector<string> funcName;
        vector< vector<real> > funcArgs;
        getFunctionNamesAndArgs("distFunc",funcName,funcArgs);
        if(funcName.size() != funcArgs.size() || funcName.size() <= 0 || funcArgs.size() <= 0) {
            ERRORMSG2("error when reading population distribution function",idStr);
            doabort();
        } else {
            for(unsigned int i = 0; i < funcName.size(); ++i) {
                args.distFunc.name.push_back(funcName[i]);
                args.distFunc.funcArgs.push_back(funcArgs[i]);
            }
        }
    }
    args.distFunc.given = var->updatedFromFile;

    return args;
}

//! Simulation diagnostics
#ifndef NO_DIAGNOSTICS
Diagnostics Params::diag;
#endif

//! Maximum number of detectors
const int Params::MAX_DETECTORS = 200;

//! Number of detectors
int Params::DETECTORS = 0;

//! Detectors
vector<Detector*> Params::detectors;

//! Creates detectors
DetectorFactory Params::detectorFactory;

// TEMPLATE DETECTOR VARIABLES
string Params::detector = "";
string Params::detectionFile = "";
string Params::popIdStr = "";
real Params::detectionTime[2] = {0, 0};
real Params::maxCounts = 0;
string Params::coordinateFile = "";
string Params::testParticleFile = "";
string Params::detectorFUNC = "";

//! Reset detector variables
void Params::clearDetectorVars()
{
    // TEMPLATE DETECTOR VARIABLES
    detectionFile = "";
    popIdStr = "";
    detectionTime[0] = 0;
    detectionTime[1] = 0;
    maxCounts = 0;
    coordinateFile = "";
    testParticleFile = "";
    m = 0;
    q = 0;
    detectorFUNC = "";
}

#define GETDETVAR(detVar) var = lookupVar(#detVar);  args.detVar.value = detVar; args.detVar.given = var->updatedFromFile;
#define GETDETVAR2(detVar) var = lookupVar(#detVar);  args.detVar.value[0] = detVar[0]; args.detVar.value[1] = detVar[1]; args.detVar.given = var->updatedFromFile;

//! Create argument structure for the class Detector
DetectorArgs Params::getDetectorArgsStruct()
{
    DetectorArgs args;
    struct dynamicVar *var;

    // TEMPLATE DETECTOR VARIABLES
    GETDETVAR(popIdStr);
    GETDETVAR(detectionFile);
    GETDETVAR2(detectionTime);
    GETDETVAR(maxCounts);
    GETDETVAR(coordinateFile);
    GETDETVAR(testParticleFile);
    GETDETVAR(m);
    GETDETVAR(q);

    // Set detector functions
    var = lookupVar("detectorFUNC");
    if(var->updatedFromFile == true) {
        vector<string> funcName;
        vector< vector<real> > funcArgs;
        getFunctionNamesAndArgs("detectorFUNC",funcName,funcArgs);
        if( funcName.size() == funcArgs.size() ) {
            args.detectorFUNC.name = funcName;
            args.detectorFUNC.funcArgs = funcArgs;
        } else {
            ERRORMSG("error when reading detector functions");
            doabort();
        }
    }
    args.detectorFUNC.given = var->updatedFromFile;

    return args;
}

// TEMPLATE PROCESS VARIABLES
string Params::process = "";
string Params::procIdStr = "";
string Params::incidentIonIdStr = "";
string Params::exoNeutralCoronaIdStr = "";
string Params::ENAIdStr = "";
string Params::slowIonIdStr = "";
fastreal Params::crossSection = 0.0;
fastreal Params::weightFactor = 0.0;
fastreal Params::slowIonVth = 0.0;
real Params::N_limitHeavyReactions = 0.0;
real Params::probLimitHeavyReactions = 0.0;

//! Reset process variables
void Params::clearProcessVars()
{
    // TEMPLATE PROCESS VARIABLES
    procIdStr = "";
    incidentIonIdStr = "";
    exoNeutralCoronaIdStr = "";
    ENAIdStr = "";
    slowIonIdStr = "";
    crossSection = 0.0;
    weightFactor = 0.0;
    slowIonVth = 0.0;
    N_limitHeavyReactions = 0.0;
    probLimitHeavyReactions = 0.0;
}

#define GETPROCVAR(procVar) var = lookupVar(#procVar); args.procVar.value = procVar; args.procVar.given = var->updatedFromFile;

//! Create argument structure for the class Process
ProcessArgs Params::getProcessArgsStruct()
{
    ProcessArgs args;
    struct dynamicVar *var;

    // TEMPLATE PROCESS VARIABLES
    GETPROCVAR(procIdStr);
    GETPROCVAR(incidentIonIdStr);
    GETPROCVAR(exoNeutralCoronaIdStr);
    GETPROCVAR(ENAIdStr);
    GETPROCVAR(slowIonIdStr);
    GETPROCVAR(crossSection);
    GETPROCVAR(weightFactor);
    GETPROCVAR(slowIonVth);
    GETPROCVAR(N_limitHeavyReactions);
    GETPROCVAR(probLimitHeavyReactions);

    return args;
}

//! Object ID for the Hybrid Web Archive
int Params::objectIdHWA = 0;

//! Planet radius [m] (initial value constant)
real Params::R_P = 0;

//! Planet mass [kg] (initial value constant)
real Params::M_P = 0;

//! Fields U_e and U are put explicitly to zero inside this radius [m]
real Params::R_zeroFields = 0;

//! Square of the above
real Params::R_zeroFields2 = 0;

//! Polarization electric field is neglected inside this radius [m] (real)
real Params::R_zeroPolarizationField = 0;

//! Field propagation timestep - should be mostly same as dt
real Params::dtField = 0;

//! Field propagation on/off
bool Params::propagateField = true;

//! Maximum speed for an ion [m/s]
real Params::vi_max = 0;

//! Square of the above
real Params::vi_max2 = 0;

//! Maximum speed for electron fluid [m/s]
real Params::Ue_max = 0;

//! Square of the above
real Params::Ue_max2 = 0;

//! Minimum charge density in a cell, rho_q = max(rho_q, rho_q_min) [C/m^3]
real Params::rho_q_min = 0;

//! Maximum whistler wave speed [m/s]
real Params::maxVw = 0;

//! Use predictor corrector scheme in Faraday's law [-]
bool Params::fieldPredCor = 0;

//! Include electron pressure term in the electric field [-]
bool Params::electronPressure = 0;
//! Electron temperature [K]
real Params::Te = 0.;

//! Gravitational acceleration for ions [-]
bool Params::useGravitationalAcceleration = false;

//! Calculation of J using jstag scheme [-]
bool Params::useJstag = false;

//! Calculation of Ue separately in nodes [-]
bool Params::useNodeUe = false;

//! Parameter for gravitational acceleration
real Params::GMdt = 0;

//! Simulation timestep [s]
real Params::dt = 0;

//! Simulation time [s]
real Params::t = 0;

//! Simulation run duration [s]
real Params::t_max = 0;

//! Number of timesteps taken in simulation. [-]
int Params::cnt_dt;

//! Save interval for output files [s]
real Params::saveInterval = 0;

//! Whether to save HC files (0 = no, 1 = binary, 2 = ascii) [-]
int Params::saveHC = 0;

//! Whether to save VTK files (0 = no, 1 = binary, 2 = ascii) [-]
int Params::saveVTK = 0;

//! "Whether to save (1) or not (0) temporally averaged parameters [-]"
bool Params::averaging = 0;

//! Whether to save (1) or not (0) plasma hc-file [-]
bool Params::plasma_hcfile = 0;

//! Whether to save (1) or not (0) debug hc-file [-]
bool Params::dbug_hcfile = 0;

//! Include bg charge density (1) or not (0) in average hc-file [-]
bool Params::bg_in_avehcfile = 0;

// Save extra hc-files [-]
bool Params::saveExtraHcFiles = 0;

//! Breakpointing intervals (first = cyclic, second = unique file names) - PRODUCES LARGE FILES! [s]
real Params::wsDumpInterval[2] = {0,0};

//! Input parameter update interval [s]
real Params::inputInterval = 0;

//! Logging interval [s]
real Params::logInterval = 0;

#ifdef SAVE_PARTICLES_ALONG_ORBIT

//! Whether to save particles along a given spacecraft orbit (creates two files: particles_along_orbit_cellindices.dat and particles_along_orbit.dat) [-]
bool Params::saveParticlesAlongOrbit = false;

//! Orbit file name if saving particles along a spacecraft orbit (three columns: x, y, z in meters) [-]
string Params::saveParticlesAlongOrbitFile = "";

#endif

#ifdef SAVE_PARTICLE_CELL_SPECTRA
real Params::spectraEmin_eV = 0;
real Params::spectraEmax_eV = 0;
int Params::spectraNbins = 0;
bool Params::spectraLogBins = true;
bool Params::spectraEminAll = true;
bool Params::spectraEmaxAll = true;
int Params::spectraMethod = 1;
string Params::spectraUnit = "1/(m^2*s*eV*str)";
vector<real> Params::spectraEnergyBins_eV;
vector<real> Params::spectra_dE_eV;
vector< vector<real> > Params::spectraV2BinsPerPop;
#endif

//! Simulation box: minimum x [m] (initial value constant)
fastreal Params::box_xmin = 0;

//! Simulation box: maximum x [m] (initial value constant)
fastreal Params::box_xmax = 0;

//! Simulation box: minimum y [m] (initial value constant)
fastreal Params::box_ymin = 0;

//! Simulation box: maximum y [m] (initial value constant)
fastreal Params::box_ymax = 0;

//! Simulation box: minimum z [m] (initial value constant)
fastreal Params::box_zmin = 0;

//! Simulation box: maximum z [m] (initial value constant)
fastreal Params::box_zmax = 0;

//! Simulation box: tight limit factor [-] (initial value constant)
fastreal Params::box_eps = 1.0e-3;

//! Simulation box: size in x-direction [m] (initial value constant)
fastreal Params::box_X;

//! Simulation box: size in y-direction [m] (initial value constant)
fastreal Params::box_Y;

//! Simulation box: size in z-direction [m] (initial value constant)
fastreal Params::box_Z;

//! Simulation box: volume [m^3] (initial value constant)
fastreal Params::box_V;

//! Simulation box: tight minimum x [m] (initial value constant)
fastreal Params::box_xmin_tight;

//! Simulation box: tight maximum x [m] (initial value constant)
fastreal Params::box_xmax_tight;

//! Simulation box: tight minimum y [m] (initial value constant)
fastreal Params::box_ymin_tight;

//! Simulation box: tight maximum y [m] (initial value constant)
fastreal Params::box_ymax_tight;

//! Simulation box: tight minimum z [m] (initial value constant)
fastreal Params::box_zmin_tight;

//! Simulation box: tight maximum z [m] (initial value constant)
fastreal Params::box_zmax_tight;

//! Simulation box: tight size in x-direction [m] (initial value constant)
fastreal Params::box_X_tight;

//! Simulation box: tight size in x-direction [m] (initial value constant)
fastreal Params::box_Y_tight;

//! Simulation box: tight size in x-direction [m] (initial value constant)
fastreal Params::box_Z_tight;

//! Base grid side length [m] (initial value constant)
real Params::dx = 0;

//! Number of density variable smoothing.
int Params::densitySmoothingNumber = 0;

//! Number of electric field smoothing.
int Params::electricFieldSmoothingNumber = 0;

//! Maximum allowed grid refinement level [-] (initial value constant)
int Params::maxGridRefinementLevel = 0;

//! Current grid refinement level (calculated by the Tgrid::Refine) [-]
int Params::currentGridRefinementLevel = 0;

//! Grid Refinement function in Refine_resis.cpp/h
string Params::gridRefinementFUNC = "{ }";

//! Forbid split and join (spatial) function
string Params::forbidSplitAndJoinFUNC = "{ }";

//! Background charge density [-]
string Params::bgChargeDensityFUNC = "{ }";

//! Deviation allowed for split and join and the probability method (0=old,1=new)
/** Old method: The real deviation around macrosPerCell where probability is 0 is
 *	about macrosPerCell*(splitJoinDeviation + 1/(1-dev)))
 *     splitJoinDeviation should be <= 0.5.
 *     Also Negative values allowed.
 *  New method: no splitting or joining is done when the number of particles in
 *     a cell is in range [macrosPerCell*(1-dev),macrosPerCell*(1+dev)]
 *     otherwise 1 to 3 splits or joins are made (attempted)
 */
real Params::splitJoinDeviation[2] = {0.2, 0};

//! Split and join parameter a for old probability method
real Params::splitjoin_a;

//! Number of cells in x-direction including ghosts [#] (initial value constant)
int Params::nx;

//! Number of cells in y-direction including ghosts [#] (initial value constant)
int Params::ny;

//! Number of cells in z-direction including ghosts [#] (initial value constant)
int Params::nz;

//! Average amount of macroparticles per cell [#]
int Params::macroParticlesPerCell = 0;

//! Macro particle splitting [-]
bool Params::useMacroParticleSplitting = true;

//! Macro particle joining [-]
bool Params::useMacroParticleJoining = true;

//! Split
Split Params::splittingFunction;

//! Join
Join Params::joiningFunction;

//! Macroparticle splitting function
string Params::splitFUNC = "{ }";

//! Macroparticle joining function
string Params::joinFUNC = "{ }";

//! Resistivity function [-]
string Params::resistivityFUNC = "{ }";

//! Interplanetary Magnetic Field x-component [T]
real Params::SW_Bx = 0;

//! Interplanetary Magnetic Field y-component [T]
real Params::SW_By = 0;

//! Interplanetary Magnetic Field z-component [T]
real Params::SW_Bz = 0;

/** \brief Magnetic field limit for program termination [T]
 *
 * Program execution is terminated when somewhere in the simulation
 * box this limit is reached.
 */
real Params::B_limit = 0;

//! Electric field cut value. dx*max(dB/dt) ~ Ecut [V/m]
real Params::Ecut = 0;

/** \brief Magnetic field boundary conditions [-]
 *
 * If set to true (1), magnetic field in the ghost cells at the wall
 * is set to boundary_Bi (i = x,y,z) values, except the normal
 * component to the wall. If set to false (0) Neumann boundary
 * conditions are applied: magnetic field in the ghost cells at the
 * wall is copied from the nearest cells in the simulation volume.
 *
 * 1st element: frontwall (x_max)
 * 2nd element: y-walls
 * 3rd element: z-walls
 * 4th element: backweall (x_min)
 */
bool Params::Bboundaries[5] = {0, 0, 0, 0, 0};

//! Magnetic field at the box boundary [T]
real Params::boundary_Bx = 0;
real Params::boundary_By = 0;
real Params::boundary_Bz = 0;

//! Initial magnetic field functions [-]
string Params::initialMagneticFieldFUNC = "{ }";

//! Constant magnetic field functions [-]
string Params::constantMagneticFieldFUNC = "{ }";

// Titan specific
real Params::SaturnLocalTime;
real Params::SubSolarLatitude;

/** \brief Set initial values for dependant variables
 *
 * After reading the initial input parameters from the config file
 * this function calculates several dependant quantities from them.
 */
void Params::setInitialValues()
{
    box_X = box_xmax - box_xmin;
    box_Y = box_ymax - box_ymin;
    box_Z = box_zmax - box_zmin;
    // Calculate amount of grid cells and round it upwards to the
    // nearest integer (ceiling).
    if(dx > 0) {
        nx = (int)( ceil(box_X / dx) );
        ny = (int)( ceil(box_Y / dx) );
        nz = (int)( ceil(box_Z / dx) );
    }
    // Update box max limits, since ceiling rounds nx, ny and nz updwards
    box_xmax = box_xmin + nx*dx;
    box_ymax = box_ymin + ny*dx;
    box_zmax = box_zmin + nz*dx;
    // Update also box measures
    box_X = box_xmax - box_xmin;
    box_Y = box_ymax - box_ymin;
    box_Z = box_zmax - box_zmin;
    // Calculate box volume
    box_V = box_X * box_Y * box_Z;
    // Calculate tight box measures
    box_xmin_tight = box_xmin + box_eps * dx;
    box_xmax_tight = box_xmax - box_eps * dx;
    box_ymin_tight = box_ymin + box_eps * dx;
    box_ymax_tight = box_ymax - box_eps * dx;
    box_zmin_tight = box_zmin + box_eps * dx;
    box_zmax_tight = box_zmax - box_eps * dx;
    box_X_tight = box_xmax_tight - box_xmin_tight;
    box_Y_tight = box_ymax_tight - box_ymin_tight;
    box_Z_tight = box_zmax_tight - box_zmin_tight;
    // Prepare counters
    cnt_dt = 0;
    updateDependantParameters();
#ifdef USE_PARTICLE_SUBCYCLING
    //! Calculate tables for subcycle dt and accumulation factors
    dt_psub[0] = dt;
    accum_psubfactor[0] = 1;
    for(int i = 1; i < 256; i++) {
        if(subcycleType == 1) {
            dt_psub[i] = dt/(i+1.0);
            accum_psubfactor[i] = 1.0/(i+1.0);
        } else if(subcycleType == 2) {
            dt_psub[i] = dt/pow(2.0,i);
            accum_psubfactor[i] = 1.0/pow(2.0,i);
        } else {
            dt_psub[i] = dt;
            accum_psubfactor[i] = 1;
        }
    }
#endif
}

//! Check if coordinates within the simulation box
bool Params::insideBox(const gridreal coords[3])
{
    if (coords[0]<Params::box_xmin || coords[0]>Params::box_xmax ||
        coords[1]<Params::box_ymin || coords[1]>Params::box_ymax ||
        coords[2]<Params::box_zmin || coords[2]>Params::box_zmax) {
        return false;
    }
    return true;
}

//! Check if particle within the tight box
bool Params::insideBoxTight(const TLinkedParticle *part)
{
    if (part->x<Params::box_xmin_tight || part->x>Params::box_xmax_tight ||
        part->y<Params::box_ymin_tight || part->y>Params::box_ymax_tight ||
        part->z<Params::box_zmin_tight || part->z>Params::box_zmax_tight) {
        return false;
    }
    return true;
}
bool Params::insideBoxTightFrontWall(const TLinkedParticle *part)
{
    if(part->x>Params::box_xmax_tight) {
        return false;
    }
    return true;
}
bool Params::insideBoxTightBackWall(const TLinkedParticle *part)
{
    if(part->x<Params::box_xmin_tight) {
        return false;
    }
    return true;
}
bool Params::insideBoxTightSideWall(const TLinkedParticle *part)
{
    if(part->y<Params::box_ymin_tight || part->y>Params::box_ymax_tight ||
       part->z<Params::box_zmin_tight || part->z>Params::box_zmax_tight) {
        return false;
    }
    return true;
}

bool Params::insideBoxTight(const gridreal coords[3])
{
    if (coords[0]<Params::box_xmin_tight || coords[0]>Params::box_xmax_tight ||
        coords[1]<Params::box_ymin_tight || coords[1]>Params::box_ymax_tight ||
        coords[2]<Params::box_zmin_tight || coords[2]>Params::box_zmax_tight) {
        return false;
    }
    return true;
}

//! Update parameters that depend on the config file inputs
void Params::updateDependantParameters()
{
    if(R_zeroFields >= 0) {
        R_zeroFields2 = sqr(R_zeroFields);
    } else {
        R_zeroFields2 = -1.0*sqr(R_zeroFields);
    }
    vi_max2 = sqr(vi_max);
    Ue_max2 = sqr(Ue_max);
    Params::GMdt = G*M_P*dt;
    // warnings apply to both probability methods
    if( (macroParticlesPerCell*(1.0-splitJoinDeviation[0])-1.0) < 0.5 ) {
        WARNINGMSG2("splitJoinDeviation too large, setting split&join off",splitJoinDeviation[0]);
        splitJoinDeviation[0] = 1.0 - 1.5/macroParticlesPerCell;
        useMacroParticleSplitting = false;
        useMacroParticleJoining = false;
    } else if( splitJoinDeviation[0] < 0 ) {
        WARNINGMSG2("splitJoinDeviation[0] is negative",splitJoinDeviation[0]);
    }
    Params::splitjoin_a = (macroParticlesPerCell-1)/(macroParticlesPerCell*(1.0-splitJoinDeviation[0])-1.0);
    if (splitJoinDeviation[1]!=0 && splitJoinDeviation[1]!=1) {
        WARNINGMSG2("splitJoinDeviation[1] must be 0 (old method) or 1, setting s/j off",splitJoinDeviation[1]);
        useMacroParticleSplitting = false;
        useMacroParticleJoining = false;
    }
    // Field propagation on/off
    if(dtField <= 0) {
        propagateField = false;
    } else {
        propagateField = true;
    }
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    // set spectra energy bins
    if(spectraNbins > 0 && spectraEmax_eV > 0 && spectraEmin_eV >=0  && spectraEmax_eV > spectraEmin_eV) {
        spectraEnergyBins_eV.clear();
        spectraV2BinsPerPop.clear();
        if(spectraLogBins == true) { // log binning
            if(spectraEmin_eV < 1.0e-10) {
                WARNINGMSG("Emin cannot be zero or too small for logarithmic energy binning, using spectraEmin_eV = 1e-10*spectraEmax_eV");
                spectraEmin_eV = 1.0e-10*spectraEmax_eV;
            }
            real dE_eV = (log10(spectraEmax_eV) - log10(spectraEmin_eV))/spectraNbins;
            for(int i=0; i<spectraNbins+1; ++i) {
                real E_edge = log10(spectraEmin_eV) + i*dE_eV;
                spectraEnergyBins_eV.push_back( pow(10.0,E_edge) );
            }
        } else { // linear binning
            real dE_eV = (spectraEmax_eV - spectraEmin_eV)/spectraNbins;
            for(int i=0; i<spectraNbins+1; ++i) {
                real E_edge = spectraEmin_eV + i*dE_eV;
                spectraEnergyBins_eV.push_back(E_edge);
            }
        }
        spectra_dE_eV.clear();
        for(int i=0; i<spectraNbins; ++i) {
            spectra_dE_eV.push_back(spectraEnergyBins_eV[i+1] - spectraEnergyBins_eV[i]);
        }
        for(int i=0; i<POPULATIONS; ++i) {
            vector<real> temp;
            spectraV2BinsPerPop.push_back(temp);
            for(unsigned int j=0; j<spectraEnergyBins_eV.size(); ++j) {
                // E_eV = 0.5*m*v^2/q_e => v^2 = 2*E_kev*q_e/m
                spectraV2BinsPerPop[i].push_back(2*spectraEnergyBins_eV[j]*Params::e/pops[i]->m);
            }
        }
    } else {
        if(Params::t > 0) {
            WARNINGMSG("cannot determine cell spectra energy binning");
        } else {
            ERRORMSG("cannot determine cell spectra energy binning");
            doabort();
        }
    }
#endif
}


/** \brief Introduces input variables and constants to the configuration system
 *
 * Input variables need to be introduced here, since otherwise they are not
 * recognised by the configuration system. Constant values should be marked
 * flagged using makeConstant("varname") and initial value constant using
 * makeInitConstant("varname"). Action functions (dependant variable updates)
 * are assigned using addAction("varname", actionFunction).
 */
void Params::initVariables()
{
    ADD_REAL(Pi,"Pi []");
    makeConstant("Pi");

    ADD_REAL(c, "Speed of light in vacuum [m/s]");
    makeConstant("c");

    ADD_REAL(e,"Electron charge magnitude [C]");
    makeConstant("e");

    ADD_REAL(k_B, "Boltzmann constant [J/K]");
    makeConstant("k_B");

    ADD_REAL(mu_0, "Permeability of free space [N A^-2]");
    makeConstant("mu_0");

    ADD_REAL(eps_0, "Permittivity of free space [F/m]");
    makeConstant("eps_0");

    ADD_REAL(G, "Gravitational constant [m^3 kg^-1 s^-2]");
    makeConstant("G");

    ADD_REAL(AU, "Astronomical unit [m]");
    makeConstant("AU");

    ADD_REAL(gamma, "Adiabatic constant []");
    makeConstant("gamma");

    ADD_REAL(amu, "Atomic mass unit [kg]");
    makeConstant("amu");

    ADD_REAL(m_e, "Electron mass [kg]");
    makeConstant("m_e");

    ADD_REAL(m_p, "Proton mass [kg]");
    makeConstant("m_p");

    ADD_REAL(m_H, "Atomic hydrogen mass [kg]");
    makeConstant("m_H");

    ADD_REAL(m_H2, "Hydrogen molecule mass [kg]");
    makeConstant("m_H2");

    ADD_REAL(m_2H, "Atomic deuterium (hydrogen-2) mass [kg]");
    makeConstant("m_2H");

    ADD_REAL(m_3H, "Atomic tritium (hydrogen-3) mass [kg]");
    makeConstant("m_3H");

    ADD_REAL(m_3He, "Atomic helium-3 mass [kg]");
    makeConstant("m_3He");

    ADD_REAL(m_He, "Atomic helium mass [kg]");
    makeConstant("m_He");

    ADD_REAL(m_N, "Atomic nitrogen mass [kg]");
    makeConstant("m_N");

    ADD_REAL(m_N2, "Molecular nitrogen mass [kg]");
    makeConstant("m_N2");

    ADD_REAL(m_O, "Atomic oxygen mass [kg]");
    makeConstant("m_O");

    ADD_REAL(m_O2, "Molecular oxygen mass [kg]");
    makeConstant("m_O2");

    ADD_REAL(m_CH4, "Molecular methane mass [kg]");
    makeConstant("m_CH4");

    ADD_REAL(m_Na, "Atomic sodium (natrium) mass [kg]");
    makeConstant("m_Na");

    ADD_REAL(M_Me, "Mercury mass [kg]");
    makeConstant("M_Me");

    ADD_REAL(R_Me, "Mercury radius [m]");
    makeConstant("R_Me");

    ADD_REAL(M_V, "Venus mass [kg]");
    makeConstant("M_V");

    ADD_REAL(R_V, "Venus radius [m]");
    makeConstant("R_V");

    ADD_REAL(M_Mo, "Moon mass [kg]");
    makeConstant("M_Mo");

    ADD_REAL(R_Mo, "Moon radius [m]");
    makeConstant("R_Mo");

    ADD_REAL(M_Ma, "Mars mass [kg]");
    makeConstant("M_Ma");

    ADD_REAL(R_Ma, "Mars radius [m]");
    makeConstant("R_Ma");

    ADD_REAL(M_T, "Titan mass [kg]");
    makeConstant("M_T");

    ADD_REAL(R_T, "Titan radius [m]");
    makeConstant("R_T");

    ADD_REAL(M_Pl, "Pluto mass [kg]");
    makeConstant("M_Pl");

    ADD_REAL(R_Pl, "Pluto radius [m]");
    makeConstant("R_Pl");

    ADD_REAL(t, "Simulation time [s]");
    makeConstant("t");

    ADD_REAL(tempRealA,"Template real variable [-]");
    ADD_REAL(tempRealB,"Template real variable [-]");
    ADD_REAL(tempRealC,"Template real variable [-]");
    ADD_INT(tempIntA,"Template int variable [-]");
    ADD_INT(tempIntB,"Template int variable [-]");
    ADD_INT(tempIntC,"Template int variable [-]");

    ADD_INT(objectIdHWA,"Object ID for the Hybrid Web Archive: 0 = Unspecified, 1 = Mercury, 2 = Venus, 3 = Earth, 4 = Moon, 5 = Mars, 6 = Asteroid, 7 = Ganymede, 8 = Titan, 9 = Pluto, 10 = Comet, 11 = Heliosphere/solar wind, 12 = Exoplanet [-]");
    makeInitConstant("objectIdHWA");

    ADD_REAL(R_P, "Planet radius [m]");
    makeInitConstant("R_P");

    ADD_REAL(M_P, "Planet mass [kg]");
    makeInitConstant("M_P");

    // TEMPLATE POPULATION VARIABLES
    ADD_POPULATION(population, "");
    setVarDumppingOff("population");

    ADD_STRING(idStr,"-");
    setVarDumppingOff("idStr");

    ADD_STRING(hcFilePrefix,"-");
    setVarDumppingOff("hcFilePrefix");

    ADD_FUNCTION(boundaryFUNC,"-");
    setVarDumppingOff("boundaryFUNC");

    ADD_REAL(m,"-");
    setVarDumppingOff("m");

    ADD_REAL(q,"-");
    setVarDumppingOff("q");

    ADD_REAL(T,"-");
    setVarDumppingOff("T");

    ADD_REAL(vth,"-");
    setVarDumppingOff("vth");

    ADD_REAL(macroParticlesPerDt,"-");
    setVarDumppingOff("macroParticlesPerDt");

    ADD_BOOL(propagateV,"-");
    setVarDumppingOff("propagateV");

    ADD_BOOL(accumulate,"-");
    setVarDumppingOff("accumulate");

    ADD_BOOL(split,"-");
    setVarDumppingOff("split");

    ADD_BOOL(join,"-");
    setVarDumppingOff("join");

    ADD_BOOL(logParams,"-");
    setVarDumppingOff("logParams");

    ADD_REAL(n,"-");
    setVarDumppingOff("n");

    ADD_REAL(V,"-");
    setVarDumppingOff("V");

    ADD_REAL(backWallWeight,"-");
    setVarDumppingOff("backWallWeight");

    ADD_REAL(R,"-");
    setVarDumppingOff("R");

    ADD_REAL(totalRate,"-");
    setVarDumppingOff("totalRate");

    ADD_FUNCTION(distFunc,"-");
    setVarDumppingOff("distFunc");

    ADD_REAL(subcycleSteps,"-");
    setVarDumppingOff("subcycleSteps");

    // TEMPLATE DETECTOR VARIABLES
    ADD_DETECTOR(detector, "");
    setVarDumppingOff("detector");

    ADD_STRING(popIdStr,"-");
    setVarDumppingOff("popIdStr");

    ADD_STRING(detectionFile,"-");
    setVarDumppingOff("detectionFile");

    ADD_REAL_TBL(detectionTime,"-", 2);
    setVarDumppingOff("detectionTime");

    ADD_REAL(maxCounts,"-");
    setVarDumppingOff("maxCounts");

    ADD_STRING(coordinateFile,"-");
    setVarDumppingOff("coordinateFile");

    ADD_STRING(testParticleFile,"-");
    setVarDumppingOff("testParticleFile");

    ADD_FUNCTION(detectorFUNC,"-");
    setVarDumppingOff("detectorFUNC");

    // TEMPLATE PROCESS VARIABLES
    ADD_PROCESS(process, "");
    setVarDumppingOff("process");

    ADD_STRING(procIdStr,"-");
    setVarDumppingOff("procIdStr");

    ADD_STRING(incidentIonIdStr,"-");
    setVarDumppingOff("incidentIonIdStr");

    ADD_STRING(exoNeutralCoronaIdStr,"-");
    setVarDumppingOff("exoNeutralCoronaIdStr");

    ADD_STRING(ENAIdStr,"-");
    setVarDumppingOff("ENAIdStr");

    ADD_STRING(slowIonIdStr,"-");
    setVarDumppingOff("slowIonIdStr");

    ADD_FASTREAL(crossSection,"-");
    setVarDumppingOff("crossSection");

    ADD_FASTREAL(weightFactor,"-");
    setVarDumppingOff("weightFactor");

    ADD_FASTREAL(slowIonVth,"-");
    setVarDumppingOff("slowIonVth");

    ADD_REAL(N_limitHeavyReactions,"-");
    setVarDumppingOff("N_limitHeavyReactions");

    ADD_REAL(probLimitHeavyReactions,"-");
    setVarDumppingOff("probLimitHeavyReactions");

    ADD_REAL(R_zeroFields, "Fields U_e and U are put explicitly to zero inside this radius [m]");
    ADD_REAL(R_zeroPolarizationField, "Polarization electric field is neglected inside this radius [m]");
    ADD_BOOL(fieldPredCor, "Field propagation using predictor corrector scheme []");
    ADD_BOOL(electronPressure, "Include electron pressure term in the electric field []");
    ADD_REAL(Te, "Electron temperature [K]");
    ADD_BOOL(useGravitationalAcceleration, "Gravitational acceleration for ions [-]");
    ADD_BOOL(useJstag, "Calculation of J using jstag scheme [-]");
    ADD_BOOL(useNodeUe, "Calculation of Ue separately in nodes [-]");
    ADD_REAL(dt, "Simulation timestep [s]");
    ADD_REAL(dx, "Base grid cell size [m]");
    makeInitConstant("dx");
    ADD_REAL(dtField, "Field propagation timestep - should be mostly same as dt [s]");
    ADD_REAL(vi_max, "Constraint: maximum ion velocity [m/s]");
    ADD_REAL(Ue_max, "Constraint: maximum electron velocity [m/s]");
    ADD_REAL(rho_q_min, "Constraint: minimum charge density in a cell, rho_q = max(rho_q, rho_q_min) [C/m^3]");
    ADD_REAL(maxVw, "Constraint: maximum whistler wave speed [m/s]");
    ADD_REAL(t_max, "Duration of simulation run [s]");
    ADD_REAL(saveInterval, "Save interval for output files [s]");
    ADD_INT(saveHC, "Whether to save HC files (0 = no, 1 = binary, 2 = ascii) [-]");
    ADD_INT(saveVTK, "Whether to save VTK files (0 = no, 1 = binary, 2 = ascii) [-]");
    ADD_BOOL(averaging, "Whether to save (1) or not (0) temporally averaged parameters [-]");
    ADD_BOOL(plasma_hcfile, "Whether to save (1) or not (0) plasma hc-file [-]");
    ADD_BOOL(dbug_hcfile, "Whether to save (1) or not (0) dbug hc-file [-]");
    ADD_BOOL(bg_in_avehcfile, "Include bg charge density (1) or not (0) in average hc-file [-]");
    ADD_BOOL(saveExtraHcFiles, "Save extra hc-files [-]");
    makeInitConstant("saveExtraHcFiles");
    ADD_REAL_TBL(wsDumpInterval, "Breakpointing intervals (first = cyclic, second = unique file names) - PRODUCES LARGE FILES! [s]",2);
    ADD_REAL(inputInterval, "Input parameter dynamics interval [s]");
    ADD_REAL(logInterval, "Logging interval [s]");
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    ADD_BOOL(saveParticlesAlongOrbit,"Particles along orbit: Whether to save particles along a given spacecraft orbit (creates two files: particles_along_orbit_cellindices.dat and particles_along_orbit.dat)  [-]");
    ADD_STRING(saveParticlesAlongOrbitFile,"Particles along orbit: Orbit file name if saving particles along a spacecraft orbit (three columns: x, y, z in meters) [-]");
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    ADD_REAL(spectraEmin_eV,"Particle cell spectra: Emin [eV]");
    ADD_REAL(spectraEmax_eV,"Particle cell spectra: Emax [eV]");
    ADD_INT(spectraNbins,"Particle cell spectra: number of bins [-]");
    ADD_BOOL(spectraLogBins,"Particle cell spectra: logarithmic binning [-]");
    ADD_BOOL(spectraEminAll,"Particle cell spectra: include all E<Emin energies in the Emin bin [-]");
    ADD_BOOL(spectraEmaxAll,"Particle cell spectra: include all E>Emax energies in the Emax bin [-]");
    ADD_INT(spectraMethod,"Particle cell spectra: method to create spectra [-]");
#endif
    ADD_FUNCTION(gridRefinementFUNC, "Grid Refinement function [-]");
    makeInitConstant("gridRefinementFUNC");
    ADD_INT(maxGridRefinementLevel, "Maximum allowed grid refinement level []");
    makeInitConstant("maxGridRefinementLevel");
    ADD_INT(densitySmoothingNumber, "Number of density variable smoothing []");
    ADD_INT(electricFieldSmoothingNumber, "Number of electric field smoothing []");
    ADD_FUNCTION(forbidSplitAndJoinFUNC, "Forbid split and join (spatial) function [-]");
    ADD_FUNCTION(bgChargeDensityFUNC, "Background charge density [-]");
    makeInitConstant("bgChargeDensityFUNC");
    ADD_INT(macroParticlesPerCell, "Average amount of macroparticles per cell [#]");
    ADD_BOOL(useMacroParticleSplitting, "Macro particle splitting [-]");
    ADD_BOOL(useMacroParticleJoining, " Macro particle joining [-]");
    ADD_REAL_TBL(splitJoinDeviation, "Deviation allowed in splitting and joining, and probability method (0=old,1=new)",2);
    ADD_FUNCTION(splitFUNC, "Macroparticle splitting function");
    ADD_FUNCTION(joinFUNC, "Macroparticle joining function");
    ADD_FUNCTION(resistivityFUNC, "Resistivity function [-]");

    ADD_FASTREAL(box_xmin, "Simulation Box: x_min [m]");
    makeInitConstant("box_xmin");

    ADD_FASTREAL(box_xmax, "Simulation Box: x_max [m]");
    makeInitConstant("box_xmax");

    ADD_FASTREAL(box_ymin, "Simulation Box: y_min [m]");
    makeInitConstant("box_ymin");

    ADD_FASTREAL(box_ymax, "Simulation Box: y_max [m]");
    makeInitConstant("box_ymax");

    ADD_FASTREAL(box_zmin, "Simulation Box: z_min [m]");
    makeInitConstant("box_zmin");

    ADD_FASTREAL(box_zmax, "Simulation Box: z_max [m]");
    makeInitConstant("box_zmax");

    ADD_FASTREAL(box_eps, "Constraint: simulation box tight limit factor []");
    makeInitConstant("box_eps");

    ADD_REAL(B_limit, "Limit for the maximum magnetic field strength [T]");
    ADD_REAL(Ecut, "Electric field cut value (used if Ecut > 0). dx*max(dB/dt) ~ Ecut [V/m]");
    ADD_REAL(SW_Bx, "Interplanetary Magnetic Field: x-component [T]");
    makeInitConstant("SW_Bx");
    ADD_REAL(SW_By, "Interplanetary Magnetic Field: y-component [T]");
    ADD_REAL(SW_Bz, "Interplanetary Magnetic Field: z-component [T]");
    ADD_BOOL_TBL(Bboundaries, "Magnetic field boundary conditions (+X, +/-Y, +/-Z, -X, perpendicular) []", 5);
    ADD_REAL(boundary_Bx, "Magnetic field x-components at the box boundary [T]");
    ADD_REAL(boundary_By, "Magnetic field y-components at the box boundary [T]");
    ADD_REAL(boundary_Bz, "Magnetic field z-components at the box boundary [T]");
    ADD_FUNCTION(initialMagneticFieldFUNC, "Initial magnetic field functions [-]");
    makeInitConstant("initialMagneticFieldFUNC");
    ADD_FUNCTION(constantMagneticFieldFUNC, "Constant magnetic field functions [-]");
    makeInitConstant("constantMagneticFieldFUNC");
    ADD_REAL(SaturnLocalTime, "SLT in hours [h]");
    ADD_REAL(SubSolarLatitude, "Declination of the Sun rel. to Titan's orb. plane [deg]");
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    ADD_REAL(sph_dy, "Base grid cell size [m]");
    makeInitConstant("sph_dy");
    ADD_REAL(sph_dz, "Base grid cell size [m]");
    makeInitConstant("sph_dz");
    // Truncation angles along theta and phi
    ADD_FASTREAL(sph_alpha_theta, "Simulation Box");
    makeInitConstant("sph_alpha_theta");

    ADD_FASTREAL(sph_alpha_phi, "Simulation Box");
    makeInitConstant("sph_alpha_phi");

    ADD_FASTREAL(sph_theta_min, "Simulation Box");
    makeInitConstant("sph_theta_min");

    ADD_FASTREAL(sph_theta_max, "Simulation Box");
    makeInitConstant("sph_theta_max");

    ADD_FASTREAL(sph_phi_min, "Simulation Box");
    makeInitConstant("sph_phi_min");

    ADD_FASTREAL(sph_phi_max, "Simulation Box");
    makeInitConstant("sph_phi_max");

    ADD_BOOL(sph_BC_type, "Type of spherical boundary conditions");
    ADD_BOOL(sph_BC_use_ghost_cell, "Type of spherical boundary conditions. Use or not ghost cells");
    ADD_BOOL(sph_propagation_type, "Type of particles propagation: radial or flat");
    ADD_BOOL(sph_propagation_dir, "direction of particles propagation: along x or z axis");
    ADD_BOOL(sph_coordinate_grid_visual, "Which coordinate grid is used for visualisation. Hybrid (0) or spherical (1)");
#endif
#ifdef USE_PARTICLE_SUBCYCLING
    ADD_INT(subcycleMaxLevel, "Particle subcycling: Subcycling maximum level");
    makeInitConstant("subcycleMaxLevel");
    ADD_INT(subcycleType, "Particle subcycling: Subcycling type");
    makeInitConstant("subcycleType");
#endif
}

// =================================================================================
// ============================= CLASS PRIVATE CONTENT =============================
// =================================================================================

//! Used to prevent the creation of multiple instances of this class
bool Params::onlyOneObject = false;

//! Use Params(false) for normal construction of the object and (true) for variable dumping purposes
Params::Params()
{
    MSGFUNCTIONCALL("Params::Params");
    if(onlyOneObject == false) {
        onlyOneObject = true;
    } else {
        errorlog << "ERROR [Params::Params]: only one Params object can be defined\n";
        doabort();
    }
    MSGFUNCTIONEND("Params::Params");
}

Params::~Params()
{
    MSGFUNCTIONCALL("Params::~Params");
    struct dynamicVar *next=varList;
    struct dynamicVar *last;
    while (next) {
        last=next;
        next=next->next;
        delete last;
    }
    MSGFUNCTIONEND("Params::~Params");
}

void Params::init(const bool dumpVariables)
{
    initPhase = true;
    stoppingPhase = false;
    varList = (struct dynamicVar *)0;
    initVariables();
    popFactory.init();
}

//! Returns simulation time in milliseconds with leading zeros
string Params::getSimuTimeStr()
{
    char fn[80];
    long int t_ms = static_cast<long int>(1000*Params::t + 0.5);
    sprintf(fn, "%.8ld",t_ms);
    string res(fn);
    if(Params::stoppingPhase == true) {
        res += "_terminated";
    }
    return res;
}

// ***** Variable handling
void Params::addVariable(void *ptr, DATATYPE type, const char *name, const char *comment, int tableSize)
{
    struct dynamicVar *newVar = new struct dynamicVar;
    struct dynamicVar *last;
    newVar->varPtr.ptr_void = ptr;
    newVar->type   = type;
    newVar->tableSize = getSize(type);
    if (getSize(type) == -1) {
        errorlog << "PARAMS: Cannot add variable " << *name << ": invalid type\n";
        delete newVar;
        return;
    }
    newVar->tableSize = tableSize;
    // + 1 needed for \0 which is exluded by strlen
    newVar->name = new char[ strlen(name) +1];
    strcpy(newVar->name, name);
    // + 1 needed for \0 which is exluded by strlen
    newVar->comment = new char[ strlen(comment)+1];
    strcpy(newVar->comment, comment);
    newVar->next = (struct dynamicVar *)0;
    newVar->constant = 0;
    newVar->initconstant = 0;
    newVar->dumpping = true;
    newVar->updatedFromFile = false;
    newVar->action = (void(*)(void*))0;
    newVar->expressions = new char*[tableSize];
    last = getLast();
    if (last) {
        last->next = newVar;
    } else {
        varList = newVar;
    }
}

/** \brief Assigns action function to be called when variable updated.
 *
 * Example declaration of valid function: void action(void *value){ };
 */
void Params::addAction(const char varName[], void (*action)(void *value))
{
    struct dynamicVar *var;
    var = lookupVar(varName);
    if (var) {
        var->action=action;
    }
}

//! Makes variable constant
void Params::makeConstant(const char varName[])
{
    struct dynamicVar *var;
    var = lookupVar(varName);
    if (var) {
        var->constant=1;
    }
}

//! Makes variable initial value constant
void Params::makeInitConstant(const char varName[])
{
    struct dynamicVar *var;
    var = lookupVar(varName);
    if (var) {
        var->initconstant=1;
    }
}

//! Set variable dumpping flag off
void Params::setVarDumppingOff(const char varName[])
{
    struct dynamicVar *var;
    var = lookupVar(varName);
    if (var) {
        var->dumpping = false;
    }
}

void Params::dumpVars(const char *dumpFile)
{
    ofstream varDump;
    struct dynamicVar *next = varList;
    varDump.open(dumpFile, ios::out);
    varDump.precision(10);
    varDump << scientific << setiosflags(ios::left);
    while (next) {
        if(next->dumpping == false) {
            next = next->next;
            continue;
        }
        varDump << "# " << next->comment << " (" << getType(next->type) << ")\n";
        if(next->type != TYPE_POPULATION && next->type != TYPE_DETECTOR && next->type != TYPE_PROCESS) {
            varDump << "#";
        }
        if(next->constant) {
            varDump << "const " << next->name << " ";
        } else if(next->initconstant) {
            varDump << "iniconst " << next->name << " ";
        } else if(next->type != TYPE_POPULATION && next->type != TYPE_DETECTOR && next->type != TYPE_PROCESS) {
            varDump << next->name << " ";
        }
        for (int ii = 0; ii < next->tableSize; ii++) {
            switch (next->type) {
            case TYPE_BOOL:
                varDump << (*( next->varPtr.ptr_bool + ii));
                break;
            case TYPE_INT:
                varDump << (*(next->varPtr.ptr_int + ii));
                break;
            case TYPE_GRIDREAL:
                varDump << (*( next->varPtr.ptr_gridreal + ii));
                break;
            case TYPE_DATAREAL:
                varDump << (*( next->varPtr.ptr_datareal + ii));
                break;
            case TYPE_TPDF_ID:
                varDump << (*( next->varPtr.ptr_TPDF_ID + ii));
                break;
            case TYPE_SHORTREAL:
                varDump << (*( next->varPtr.ptr_shortreal + ii));
                break;
            case TYPE_FASTREAL:
                varDump << (*( next->varPtr.ptr_fastreal + ii));
                break;
            case TYPE_REAL:
                varDump << (*( next->varPtr.ptr_real + ii));
                break;
            case TYPE_STRING:
                varDump << (*( next->varPtr.ptr_string + ii));
                break;
            case TYPE_CHAR:
                varDump << (*( next->varPtr.ptr_char + ii));
                break;
            case TYPE_FUNCTION:
                varDump << (*( next->varPtr.ptr_string + ii));
                break;
            case TYPE_POPULATION:
                if(ii == 0) {
                    varDump << "\n";
                } else if (ii < next->tableSize - 1) {
                    varDump << "\n\n";
                    varDump << "#";
                }
                varDump << (*( next->varPtr.ptr_string + ii));
                break;
            case TYPE_DETECTOR:
                if(ii == 0) {
                    varDump << "\n";
                } else if (ii < next->tableSize - 1) {
                    varDump << "\n\n";
                    varDump << "#";
                }
                varDump << (*( next->varPtr.ptr_string + ii));
                break;
            case TYPE_PROCESS:
                if(ii == 0) {
                    varDump << "\n";
                } else if (ii < next->tableSize - 1) {
                    varDump << "\n\n";
                    varDump << "#";
                }
                varDump << (*( next->varPtr.ptr_string + ii));
                break;
            default:
                cerr << "ERROR [Params::dumpVars]: invalid variable type (" << next->name << ")\n";
                doabort();
                break;
            }
            if (ii < next->tableSize - 1 &&
                next->type != TYPE_POPULATION && next->type != TYPE_DETECTOR && next->type != TYPE_PROCESS) {
                varDump << " ";
            }
        }
        varDump << "\n\n";
        next = next->next;
    }
    // Dump populations
    for(unsigned int i=0; i < Params::pops.size(); ++i) {
        varDump << Params::pops[i]->configDump();
        varDump << "\n\n";
    }
    // Dump detectors
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        varDump << Params::detectors[i]->configDump();
        varDump << "\n\n";
    }
}

void Params::outParams()
{
    // Write header lines
    static bool writeHeader = true;
    if (writeHeader == true) {
        // General information
        // File structure
        int nn = 1;
        string separator(1,'\t');
        paramslog << "% ";
        paramslog << int2string(nn++,2) << ". t [s]         " << separator;
        paramslog << int2string(nn++,2) << ". vi_max [m/s]    " << separator;
        paramslog << int2string(nn++,2) << ". Ue_max [m/s]    " << separator;
        paramslog << int2string(nn++,2) << ". rhoqmin [C/m^3] " << separator;
        paramslog << int2string(nn++,2) << ". mples/Cell [#]  " << separator;
        paramslog << int2string(nn++,2) << ". useSplit [-]    " << separator;
        paramslog << int2string(nn++,2) << ". useJoin [-]     " << separator;
        paramslog << int2string(nn++,2) << ". Ecut [Vm]       " << separator;
        paramslog << int2string(nn++,2) << ". B_lim [T]       " << separator;
        paramslog << int2string(nn++,2) << ". SW_Bx [T]       " << separator;
        paramslog << int2string(nn++,2) << ". SW_By [T]       " << separator;
        paramslog << int2string(nn++,2) << ". SW_Bz [T]       " << separator;
        paramslog << "\n";
        paramslog << scientific << showpos;
        paramslog.precision(10);
        writeHeader = false;
    }
    string separator(3,' ');
    separator.append(1,'\t');
    paramslog << t << separator;
    paramslog << vi_max << separator;
    paramslog << Ue_max << separator;
    paramslog << rho_q_min << separator;
    paramslog << static_cast<real>(macroParticlesPerCell) << separator;
    paramslog << static_cast<real>(useMacroParticleSplitting) << separator;
    paramslog << static_cast<real>(useMacroParticleJoining) << separator;
    paramslog << Ecut << separator;
    paramslog << B_limit << separator;
    paramslog << SW_Bx << separator;
    paramslog << SW_By << separator;
    paramslog << SW_Bz << separator;
    paramslog << "\n";
}

//! Set new variable value
void Params::setVar(struct dynamicVar *var,  void * value, int index)
{
    if(var->constant) {
        return;
    }
    if(var->initconstant == true && initPhase == false) {
        return;
    }
    switch (var->type) {
    case TYPE_BOOL:
        *(var->varPtr.ptr_bool  + index) = *((bool *) value);
        break;
    case TYPE_INT:
        *(var->varPtr.ptr_int  + index) = *((int *) value);
        break;
    case TYPE_GRIDREAL:
        *(var->varPtr.ptr_gridreal  + index) = *((gridreal *) value);
        break;
    case TYPE_DATAREAL:
        *(var->varPtr.ptr_datareal  + index) = *((datareal*) value);
        break;
    case TYPE_TPDF_ID:
        *(var->varPtr.ptr_TPDF_ID  + index) = *((TPDF_ID*) value);
        break;
    case TYPE_SHORTREAL:
        *(var->varPtr.ptr_shortreal  + index ) = *((shortreal*) value);
        break;
    case TYPE_FASTREAL:
        *(var->varPtr.ptr_fastreal  + index) = *((fastreal*) value);
        break;
    case TYPE_REAL:
        *(var->varPtr.ptr_real + index)= *((real*) value);
        break;
    case TYPE_STRING:
        *(var->varPtr.ptr_string + index)= *((string*) value);
        break;
    case TYPE_CHAR:
        *(var->varPtr.ptr_char + index)= *((char*) value);
        break;
    case TYPE_FUNCTION:
        *(var->varPtr.ptr_string + index)= *((string*) value);
        break;
    case TYPE_POPULATION:
        *(var->varPtr.ptr_string + index)= *((string*) value);
        break;
    case TYPE_DETECTOR:
        *(var->varPtr.ptr_string + index)= *((string*) value);
        break;
    case TYPE_PROCESS:
        *(var->varPtr.ptr_string + index)= *((string*) value);
        break;
    default:
        errorlog << "ERROR [Params::setVar]: invalid variable type (" << var->name << ")\n";
        doabort();
        break;
    }
}

int Params::getSize(DATATYPE type)
{
    int size;
    switch (type) {
    case TYPE_BOOL:
        size = sizeof(bool);
        break;
    case TYPE_INT:
        size = sizeof(int);
        break;
    case TYPE_GRIDREAL:
        size = sizeof(gridreal);
        break;
    case TYPE_DATAREAL:
        size = sizeof(datareal);
        break;
    case TYPE_TPDF_ID:
        size = sizeof(TPDF_ID);
        break;
    case TYPE_SHORTREAL:
        size = sizeof(shortreal);
        break;
    case TYPE_FASTREAL:
        size = sizeof(fastreal);
        break;
    case TYPE_REAL:
        size = sizeof(real);
        break;
    case TYPE_STRING:
        size = sizeof(string);
        break;
    case TYPE_CHAR:
        size = sizeof(char);
        break;
    case TYPE_FUNCTION:
        size = sizeof(string);
        break;
    case TYPE_POPULATION:
        size = sizeof(string);
        break;
    case TYPE_DETECTOR:
        size = sizeof(string);
        break;
    case TYPE_PROCESS:
        size = sizeof(string);
        break;
    default:
        return -1;
    }
    return size;
}

string Params::getType(DATATYPE type)
{
    string typeStr;
    switch (type) {
    case TYPE_BOOL:
        typeStr = "boolean";
        break;
    case TYPE_INT:
        typeStr = "integer";
        break;
    case TYPE_GRIDREAL:
        typeStr = "gridreal";
        break;
    case TYPE_DATAREAL:
        typeStr = "datareal";
        break;
    case TYPE_TPDF_ID:
        typeStr = "TPDF_ID";
        break;
    case TYPE_SHORTREAL:
        typeStr = "shortreal";
        break;
    case TYPE_FASTREAL:
        typeStr = "fastreal";
        break;
    case TYPE_REAL:
        typeStr = "real";
        break;
    case TYPE_STRING:
        typeStr = "string";
        break;
    case TYPE_CHAR:
        typeStr = "char";
        break;
    case TYPE_FUNCTION:
        typeStr = "function";
        break;
    case TYPE_POPULATION:
        typeStr = "population";
        break;
    case TYPE_DETECTOR:
        typeStr = "detector";
        break;
    case TYPE_PROCESS:
        typeStr = "process";
        break;
    default:
        return "type";
    }
    return typeStr;
}

real Params::getRealValue(const char*varName, int index)
{
    struct dynamicVar *var=lookupVar(varName);
    real value;
    switch (var->type) {
    case TYPE_BOOL:
        value=real(*(var->varPtr.ptr_bool + index));
        break;
    case TYPE_INT:
        value=real(*(var->varPtr.ptr_int + index));
        break;
    case TYPE_GRIDREAL:
        value=real(*(var->varPtr.ptr_gridreal + index));
        break;
    case TYPE_DATAREAL:
        value=real(*(var->varPtr.ptr_datareal + index));
        break;
    case TYPE_TPDF_ID:
        value=real(*(var->varPtr.ptr_TPDF_ID + index));
        break;
    case TYPE_SHORTREAL:
        value=real(*(var->varPtr.ptr_shortreal + index));
        break;
    case TYPE_FASTREAL:
        value=real(*(var->varPtr.ptr_fastreal + index));
        break;
    case TYPE_REAL:
        value=*(var->varPtr.ptr_real + index);
        break;
    case TYPE_STRING:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a string variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    case TYPE_CHAR:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a char variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    case TYPE_FUNCTION:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a function variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    case TYPE_POPULATION:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a population variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    case TYPE_DETECTOR:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a detector variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    case TYPE_PROCESS:
        errorlog << "ERROR [Params::getRealValue]: trying to get numerical value of a process variable (" << var->name << ")\n";
        doabort();
        value = 0;
        break;
    default:
        errorlog << "ERROR [Params::getRealValue]: cannot recognize variable type (" << var->name << ")\n";
        doabort();
        // Dummy statement, just to get rid of error messages...
        value=0;
    }
    return value;
}

struct dynamicVar *Params::getLast()
{
    struct dynamicVar *next;
    next = varList;
    while (next && next->next) {
        next = next->next;
    }
    return next;
}

struct dynamicVar *Params::lookupVar(const char *name)
{
    struct dynamicVar *var;
    var = varList;
    while (var) {
        if (!strcmp(name, var->name)) {
            break;
        }
        var = var->next;
    }
    return var;
}

// This method returns 1 if variable has been changed really. Slightly complex
// code in it, but who cares...
int Params::readAndUpdateVariable(ifstream &varDump, struct dynamicVar *var, int index)
{
    int changed=0;
    char testStr[128];
    char expression[512];
    int pos;
    // Do not update constant
    if(var->constant) {
        return changed;
    }
    // Do not update initial value constant outside initial phase
    if(var->initconstant && initPhase == false) {
        return changed;
    }
    pos = varDump.tellg();
    varDump >> testStr;
    varDump.seekg(pos);
    // Select variables type
    switch (var->type) {
    case TYPE_BOOL: {
        bool temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = bool( readAndEvaluateRPNExpression(expression) );
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_bool + index)) {
            *(var->varPtr.ptr_bool + index) =temp;
            changed=1;
        }
    }
    break;
    case TYPE_INT: {
        int temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = int( readAndEvaluateRPNExpression(expression) );
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_int + index)) {
            *(var->varPtr.ptr_int + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_GRIDREAL: {
        gridreal temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = readAndEvaluateRPNExpression(expression);
        } else {
            varDump >> temp;
        }

        if (temp != *(var->varPtr.ptr_gridreal + index)) {
            *(var->varPtr.ptr_gridreal + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_DATAREAL: {
        datareal temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = readAndEvaluateRPNExpression(expression);
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_datareal + index)) {
            *(var->varPtr.ptr_datareal + index)  = temp;
            changed=1;
        }
    }
    break;
    case  TYPE_TPDF_ID: {
        TPDF_ID temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_TPDF_ID + index)) {
            *(var->varPtr.ptr_TPDF_ID + index)  = temp;
            changed=1;
        }
    }
    break;
    case TYPE_SHORTREAL: {
        shortreal temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = readAndEvaluateRPNExpression(expression);
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_shortreal + index)) {
            *(var->varPtr.ptr_shortreal + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_FASTREAL: {
        fastreal temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = readAndEvaluateRPNExpression(expression);
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_fastreal + index)) {
            *(var->varPtr.ptr_fastreal + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_REAL: {
        real temp;
        if(testStr[0]=='=') {
            varDump.getline(expression, 512, ';');
            temp = readAndEvaluateRPNExpression(expression);
        } else {
            varDump >> temp;
        }
        if (temp != *(var->varPtr.ptr_real + index)) {
            *(var->varPtr.ptr_real + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_STRING: {
        string temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_string + index)) {
            *(var->varPtr.ptr_string + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_CHAR: {
        char temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_char + index)) {
            *(var->varPtr.ptr_char + index) = temp;
            changed=1;
        }
    }
    break;
    case TYPE_FUNCTION: {
        string str = readFunctionTypeVar(varDump,var);
        if (str != *(var->varPtr.ptr_string + index)) {
            *(var->varPtr.ptr_string + index) = str;
            changed=1;
        }
    }
    break;
    case TYPE_POPULATION: {
        // Read population type string
        string temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_string + index)) {
            *(var->varPtr.ptr_string + index) = temp;
            changed=1;
        }
        // Get characters before the opening brace and check there are no illegal ones
        const unsigned int maxCh = 1024;
        char tempB[maxCh];
        varDump.getline(tempB,maxCh,'{');
        string str(tempB);
        if (str.find_first_not_of(" \n\t") != string::npos) {
            errorlog << "ERROR [Params::readAndUpdateVariable(" << var->name << ")]: bad character before the population opening brace (" << temp << ")\n";
            doabort();
        }
        // Set population reading flag on
        readingPopulation = true;
        resetUpdatedFromFileFlags();
        clearPopulationVars();
    }
    break;
    case TYPE_DETECTOR: {
        // Read detector type string
        string temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_string + index)) {
            *(var->varPtr.ptr_string + index) = temp;
            changed=1;
        }
        // Get characters before the opening brace and check there are no illegal ones
        const unsigned int maxCh = 1024;
        char tempB[maxCh];
        varDump.getline(tempB,maxCh,'{');
        string str(tempB);
        if (str.find_first_not_of(" \n\t") != string::npos) {
            errorlog << "ERROR [Params::readAndUpdateVariable(" << var->name << ")]: bad character before the detector opening brace (" << temp << ")\n";
            doabort();
        }
        // Set population reading flag on
        readingDetector = true;
        resetUpdatedFromFileFlags();
        clearDetectorVars();
    }
    break;
    case TYPE_PROCESS: {
        // Read population type string
        string temp;
        varDump >> temp;
        if (temp != *(var->varPtr.ptr_string + index)) {
            *(var->varPtr.ptr_string + index) = temp;
            changed=1;
        }
        // Get characters before the opening brace and check there are no illegal ones
        const unsigned int maxCh = 1024;
        char tempB[maxCh];
        varDump.getline(tempB,maxCh,'{');
        string str(tempB);
        if (str.find_first_not_of(" \n\t") != string::npos) {
            errorlog << "ERROR [Params::readAndUpdateVariable(" << var->name << ")]: bad character before the process opening brace (" << temp << ")\n";
            doabort();
        }
        // Set process reading flag on
        readingProcess = true;
        resetUpdatedFromFileFlags();
        clearProcessVars();
    }
    break;
    default:
        errorlog << "ERROR [Params::readAndUpdateVariable(" << var->name << ")]: invalid variable type\n";
        doabort();
        break;
    }
    return changed;
}

//! Reads function type variable from the stream and returns it in a string
string Params::readFunctionTypeVar(ifstream &fileStream, struct dynamicVar *var)
{
    const unsigned int maxCh = 1024;
    char temp[maxCh];
    // Get characters before the opening brace
    fileStream.getline(temp,maxCh,'{');
    // Check there are no illegal characters
    string str(temp);
    if (str.find_first_not_of(" \n\t") != string::npos) {
        errorlog << "ERROR [Params::readFunctionTypeVar(" << var->name << ")]: bad character before the opening brace (" << temp << ")\n";
        doabort();
    }
    // Begin constructing the string
    str = "{";
    // Read characters between the braces
    fileStream.getline(temp, maxCh, '}');
    // Remove unnecessary space etc. from the string
    string tempB(temp);
    tempB = cropPrecedingAndTrailingSpaces(tempB);
    // Construct the string
    str.append(tempB);
    str.append("}");
    return str;
}

//! Go thru the config file and update variables
void Params::readAndUpdateVariables(const char *dumpFile)
{
    ifstream varDump;
    char varName[100];
    int ii;
    int changed;
    struct dynamicVar *var;
    varDump.open(dumpFile);
    if (varDump.is_open()) {
        resetUpdatedFromFileFlags();
        // Reset population and process id string lists (used when checking the population or process has not been read twice)
        idStrTbl.clear();
        procIdStrTbl.clear();
        // Go thru the file and read&update variables
        do {
            varDump >> varName;
            // Skip const and iniconst prefixes
            if (!strcmp(varName ,"const") || !strcmp(varName ,"iniconst")) {
                varDump >> varName;
            }
            // Getting rid of comments and empty lines
            if (strstr(varName, "#") || !strcmp(varName, "")) {
                varDump.ignore(2048, '\n');
                continue;
            }
            // Getting rid of braced content
            if ( strstr(varName, "{") ) {
                varDump.unget();
                skipBracedContent(varDump);
                continue;
            }
            if ( strstr(varName, "}") ) {
                // Class is in the population reading mode
                if(readingPopulation == true) {
                    if(checkIdStrAlreadyFound(idStr) == true) {
                        ERRORMSG2("population already read (=same id string used twice)",idStr);
                        doabort();
                    }
                    // Check whether the population already exists
                    bool popFound = false;
                    for(unsigned int i=0; i < pops.size(); ++i) {
                        if(idStr.compare( pops[i]->getIdStr() ) == 0) {
                            // Save new population arguments
                            pops[i]->saveArgs(getPopulationArgsStruct());
                            popFound = true;
                            break;
                        }
                    }
                    // If the population does not exist yet, add it to the vector
                    if(popFound == false) {
                        PopulationArgs tempArgs = getPopulationArgsStruct();
                        Population* tempPop = popFactory.createPopulation(population,tempArgs);
                        pops.push_back(tempPop);
                    }
                    // Set population reading flag off
                    readingPopulation = false;
                }
                // Class is in the detector reading mode
                if(readingDetector == true) {
                    bool detIdStrFound = false;
                    for(unsigned int i = 0; i < idStrTbl.size(); ++i) {
                        if(idStrTbl[i].compare(popIdStr) == 0) {
                            detIdStrFound = true;
                            break;
                        }
                    }
                    if(detIdStrFound == false) {
                        ERRORMSG2("population not found for particle detectors",popIdStr);
                        doabort();
                    }
                    DetectorArgs tempArgs = getDetectorArgsStruct();
                    //Detectors created only at initialization!
                    // no updating!
                    if(initPhase == true) {
                        Detector* tempDetector = detectorFactory.createDetector(detector,tempArgs);
                        detectors.push_back(tempDetector);
                    }
                    // Set population reading flag off
                    readingDetector = false;
                }
                // Class is in the process reading mode
                if(readingProcess == true) {
                    if(checkProcIdStrAlreadyFound(procIdStr) == true) {
                        ERRORMSG2("process already read (=same proc id string used twice)",procIdStr);
                        doabort();
                    }
                    // Check whether the process already exists
                    if(ParticleProcesses::checkProcIdStrExists(procIdStr) == true) {
                        // update process args
                        ParticleProcesses::updateReaction(process,getProcessArgsStruct());
                    } else {
                        // create new process
                        ParticleProcesses::createReaction(process,getProcessArgsStruct());
                    }
                    // Set process reading flag off
                    readingProcess = false;
                }
                varDump.ignore(2048, '\n');
                continue;
            }
            var=lookupVar(varName);
            // Check if the variable exist
            if (var) {
                var->updatedFromFile = true;
                // Do not try to update constant and initial value constant outside initial phase
                if( (var->constant) || (var->initconstant == true && initPhase == false) ) {
                    // Skip braced function content
                    if (var->type == TYPE_FUNCTION) {
                        skipBracedContent(varDump);
                        continue;
                    }
                } else {
                    changed=0;
                    for( ii=0; ii<(var->tableSize); ii++) {
                        changed |= readAndUpdateVariable(varDump, var, ii);
                    }
                    if (changed) {
                        // Execute action only if not in initial phase
                        if (var->action && initPhase == false) {
                            var->action(var->varPtr.ptr_void);
                        }
                    }
                }
            } else {
                errorlog << "ERROR [Params::readAndUpdateVariables]: variable not found (" << varName << ")\n";
                doabort();
            }
            varDump.ignore(2048, '\n');
            varName[0]='\0';
            var = (struct dynamicVar *)0;
        } while( !varDump.eof());
        // Config file reading end
        // Update variable dependencies
        if(initPhase == true) {
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
            setInitialValues();
#else
            sph_setInitialValues();
#endif
        } else {
            updateDependantParameters();
        }
        // Update population argumets and cross dependencies
        for(unsigned int i=0; i < pops.size(); ++i) {
            pops[i]->updateArgs();
        }
        // Write parameter log
        outParams();
    } else {
        errorlog << "ERROR [Params::readAndUpdateVariables]: no config file available (" << dumpFile << ")\n";
        doabort();
    }
}

void Params::resetUpdatedFromFileFlags()
{
    struct dynamicVar *next = varList;
    while (next) {
        next->updatedFromFile = false;
        next = next->next;
    }
}

//! Reads variables from config file and updates them. Initial constants can be updated only during this function.
void Params::readAndInitVariables(const char *dumpFile)
{
    MSGFUNCTIONCALL("Params::readAndInitVariables");
    readAndUpdateVariables(dumpFile);
    // Set initial phase flag to false
    initPhase = false;
    MSGFUNCTIONEND("Params::readAndInitVariables");
}

//! Reverse Polish Notation expression handling
double Params::readAndEvaluateRPNExpression(const char *expression)
{
    char *rest;
    const char *line;
    double stack[64];
    char tempStr[64];
    int stackI=-1;
    int ii;
    double value, temp, temp2;
    double retval;
    int equ_seen=0;
    line=expression;
    rest=(char *)line;
    retval=0;
    while(line[0]) {
        value = strtod(line, &rest);
        if (strcmp(line, rest)) {
            //lines differ, value is valid
            stack[++stackI]=value;
        } else {
            ii=0;
            switch(line[ii]) {
            case ' ':
            case '\t':
            case '(':	//Parenthesis just to support 'commenting'
            case ')':       // in expressions -> possible to mark
            case '}':       // stack operations by two sets of
            case '{':       // parenthesis
                rest++;
                break;
            case '=':
                if (equ_seen) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack=)\n";
                    doabort();
                }
                equ_seen=1;
                rest++;
                break;
            case '+':
                if (stackI<1) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack+)\n";
                    doabort();
                }
                temp = stack[stackI--];
                temp += stack[stackI--];
                stack[++stackI]=temp;
                rest++;
                break;
            case '-':
                if (stackI<1) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack-)\n";
                    doabort();
                }
                temp = -stack[stackI--];
                temp += stack[stackI--];
                stack[++stackI]=temp;
                rest++;
                break;
            case '*':
                if (stackI<1) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack*)\n";
                    doabort();
                }
                temp = stack[stackI--];
                temp *= stack[stackI--];
                stack[++stackI]=temp;
                rest++;
                break;
            case '/':
                if (stackI<1) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack/)\n";
                    doabort();
                }
                temp2 = stack[stackI--];
                temp = stack[stackI--];
                stack[++stackI]=temp/temp2;
                rest++;
                break;
            case '^':
                if (stackI==-1) {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack^)\n";
                    doabort();
                }
                temp = stack[stackI--];
                stack[stackI] = pow(stack[stackI], temp);
                rest++;
                break;
            default:
                sscanf(line+ii,"%s", tempStr);
                rest += strlen(tempStr);
                if (!strcmp(tempStr, "sin")) {
                    stack[stackI] = sin(stack[stackI]);
                } else if (!strcmp(tempStr, "cos")) {
                    stack[stackI] = cos(stack[stackI]);
                } else if (!strcmp(tempStr, "exp")) {
                    stack[stackI] = exp(stack[stackI]);
                } else if (!strcmp(tempStr, "stepup")) {
                    if (stack[stackI] >= 0) {
                        stack[stackI]=1;
                    } else {
                        stack[stackI]=0;
                    }
                } else if (!strcmp(tempStr, "stepdown")) {
                    if (stack[stackI] >=0) {
                        stack[stackI]=0;
                    } else {
                        stack[stackI]=1;
                    }
                } else if (!strcmp(tempStr, "fetch")) {
                    // =t fetch; reads a value
                    // from <varName>.dat
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (fetch N/A yet)\n";
                    doabort();
                    //stack[stackI]
                    //    = exp(stack[stackI]);
                } else if(lookupVar(tempStr)) {
                    stack[++stackI] = getRealValue(tempStr);
                } else {
                    errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error (stack other)\n";
                    doabort();
                }
                break;
            }
        }
        line=rest;
    }
    if(stackI) {
        errorlog << "ERROR [Params::readAndEvaluateRPNExpression(\"" << expression  << ")]: parse error stack (" << stackI <<" items)\n";
        doabort();
    }
    retval =stack[0];
    return retval;
}

//! Get the function names and arguments if any. Return false if no names and arguments.
bool Params::getFunctionNamesAndArgs(const char varName[], vector<string> &funcNames, vector< vector<real> > &args)
{
    struct dynamicVar *var;
    var = lookupVar(varName);
    // Check whether the variable exists and is of function type
    if (var) {
        switch (var->type) {
        case TYPE_FUNCTION:
            // This is ok
            break;
        default:
            errorlog << "ERROR [ Params::getFunctionNamesAndArgs(" << varName <<  ") ]: cannot parse arguments for non-function type variable\n";
            doabort();
            break;
        }
    }
    // Empty result vectors
    funcNames.clear();
    args.clear();
    string str = getStringValue(varName);
    // Check the opening brace is ok
    string::size_type nBraceOpen = str.find_first_of('{');
    string temp = str.substr(0,nBraceOpen+1);
    if(nBraceOpen == string::npos || temp.find_first_not_of("{ \n\t") != string::npos) {
        errorlog << "ERROR [ Params::getFunctionNamesAndArgs(" << varName <<  ") ]: cannot find the opening brace\n";
        doabort();
    }
    // Check the closing brace is ok
    string::size_type nBraceClose = str.find_last_of('}');
    temp = str.substr(nBraceClose,str.length() - nBraceClose);
    if(nBraceClose == string::npos || nBraceOpen >= nBraceClose || temp.find_first_not_of("} \n\t") != string::npos) {
        errorlog << "ERROR [ Params::getFunctionNamesAndArgs(" << varName <<  ") ]: cannot find the closing brace\n";
        doabort();
    }
    // Check the size
    if(nBraceOpen >= nBraceClose-1) {
        funcNames.clear();
        args.clear();
        // Zero size function variable
        return false;
    }
    str = str.substr(nBraceOpen+1,nBraceClose-nBraceOpen-1);
    // Check the string is non-empty
    if(str.find_first_not_of("{} \n\t") == string::npos) {
        funcNames.clear();
        args.clear();
        // Empty function variable
        return false;
    }
    unsigned int nn = 0;
    string::size_type mm = 0;
    unsigned int ll = 0;
    // Go thru the string lines
    do {
        mm = str.find_first_of("\n",nn);
        if(mm == string::npos) {
            ll = str.length() - 1;
        } else {
            ll = mm - nn;
        }
        // Check the line is non-empty and store the name and arguments
        if (ll > 0) {
            temp = str.substr(nn,ll+1);
            if(temp.find_first_not_of("{} \n\t") != string::npos) {
                string tempFuncName;
                vector<real> tempArgs;
                tempArgs = parseFunctionNameAndArguments(temp,tempFuncName);
                // Store the function name and arguments
                funcNames.push_back(tempFuncName);
                args.push_back(tempArgs);
            }
        }
        nn = mm+1;
    } while(mm != string::npos);
    return true;
}

//! Returns function name string. In case of name multiple functions, the names are catenated using a separator.
string Params::getFunctionName(const char varName[])
{
    vector<string> funcNames;
    vector< vector<real> > args;
    bool nonZeroFuncs = getFunctionNamesAndArgs(varName, funcNames, args);
    if(nonZeroFuncs == false) {
        return "";
    }
    string result = "";
    vector<string>::iterator funcNamesIter = funcNames.begin();
    while( funcNamesIter != funcNames.end() ) {
        string str = *funcNamesIter++;
        result.append(str);
        if( funcNamesIter < funcNames.end() ) {
            result.append("_");
        }
    }
    return result;
}

//! Parse one line of function argumets
vector<real> Params::parseFunctionNameAndArguments(string str, string &funcName)
{
    vector<real> args;
    str = cropPrecedingAndTrailingSpaces(str);
    // First string entry is the function name
    string::size_type startPos = str.find_first_not_of(' ',0);
    string::size_type endPos = str.find_first_of(' ',startPos);
    if (endPos == string::npos) {
        endPos = str.length();
    }
    if (endPos > startPos && startPos != string::npos) {
        funcName = str.substr(startPos,endPos-startPos);
    } else {
        funcName = "";
    }
    // Go thru the string and store the values as reals in a vector array
    startPos = str.find_first_not_of(' ',endPos);
    real tempValue;
    while (startPos != string::npos) {
        string strTest = str.substr(startPos,1);
        // If the argument is a RPN expression (for example, "=R_P 1000e3 +;")
        if (strTest.compare("=") == 0) {
            endPos = str.find_first_of(";",startPos) + 1;
            string val = str.substr(startPos+1,endPos-startPos-2);
            tempValue = readAndEvaluateRPNExpression( val.c_str() );
        }
        // If the value is purely numerical
        else {
            endPos = str.find_first_of(' ',startPos);
            string val = str.substr(startPos,endPos-startPos);
            tempValue = string2double(val);
        }
        args.push_back(tempValue);
        startPos = str.find_first_not_of(' ',endPos);
    }
    return args;
}

//! Returns a string value of variable
string Params::getStringValue(const char*varName, int index)
{
    struct dynamicVar *var=lookupVar(varName);
    string value;
    switch (var->type) {
    case TYPE_STRING:
        value = *(var->varPtr.ptr_string + index);
        break;
    case TYPE_CHAR:
        value = *(var->varPtr.ptr_char + index);
        break;
    case TYPE_FUNCTION:
        value = *(var->varPtr.ptr_string + index);
        break;
    case TYPE_POPULATION:
        value = *(var->varPtr.ptr_string + index);
        break;
    case TYPE_DETECTOR:
        value = *(var->varPtr.ptr_string + index);
        break;
    case TYPE_PROCESS:
        value = *(var->varPtr.ptr_string + index);
        break;
    default:
        errorlog << "ERROR [Params::getStringValue]: cannot get variable string value (" << var->name << ")\n";
        doabort();
        // Dummy statement, just to get rid of error messages...
        value="";
    }
    return value;
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Base grid side length [m] (initial value constant)
real Params::sph_dy = 0;

//! (SPHERICAL) Base grid side length [m] (initial value constant)
real Params::sph_dz = 0;
real Params::sph_dtheta = 0;
real Params::sph_dphi   = 0;

//! (SPHERICAL) Truncation theta angle
fastreal Params::sph_alpha_theta = 0;

//! (SPHERICAL) Truncation phi angle
fastreal Params::sph_alpha_phi = 0;

//! (SPHERICAL) Simulation box: angle of theta axis (0; pi/2), 0 - parallel to solar wind, pi/2 - perpendicular [rad] (initial value constant)
fastreal Params::sph_theta_min = 0;
fastreal Params::sph_theta_max = 0;
fastreal Params::sph_phi_min = 0;
fastreal Params::sph_phi_max = 0;
fastreal Params::sph_box_theta = 0;
fastreal Params::sph_box_phi = 0;

//! (SPHERICAL) Spherical simulation box: volume [m^3] (initial value constant)
fastreal Params::sph_box_V;

//! (SPHERICAL) Type of spherical boundary conditions
int Params::sph_BC_type;

//! (SPHERICAL) Type of spherical boundary conditions. Use or not ghost cell
int Params::sph_BC_use_ghost_cell;

//! (SPHERICAL) Type of propagation: radial (0), flat (Cartesian) (1)
int Params::sph_propagation_type;

//! (SPHERICAL) direction of particles propagation: (0) along x or (1) z axis.
int Params::sph_propagation_dir;

//! (SPHERICAL) Which coordinate grid is used for visualisation. Hybrid (0) or spherical (1)
int Params::sph_coordinate_grid_visual;

//! (SPHERICAL) Spherical version of "setInitialValues"
void Params::sph_setInitialValues()
{
    box_X = box_xmax - box_xmin;
    box_Y = box_ymax - box_ymin;
    box_Z = box_zmax - box_zmin;
    sph_box_theta = sph_theta_max - sph_theta_min;
    sph_box_phi   = sph_phi_max   - sph_phi_min;
    // Calculate amount of grid cells and round it upwards to the
    // nearest integer (ceiling).
    if(dx > 0) {
        // Using this block there is a problem with large scale grid in spherical coordinates (gap in spherical grid near the borders)
        //nx = (int)( ceil(box_X / dx) );
        //ny = (int)( ceil(box_Y / sph_dy) );
        //nz = (int)( ceil(box_Z / sph_dz) );
        // Using this block you must be very carefull with fitting box size and number of sells (in config file)
        nx = (int)( floor(box_X / dx) );
        ny = (int)( floor(box_Y / sph_dy) );
        nz = (int)( floor(box_Z / sph_dz) );
    }
    // Update box max limits, since ceiling rounds nx, ny and nz updwards
    box_xmax = box_xmin + nx*dx;
    box_ymax = box_ymin + ny*sph_dy;
    box_zmax = box_zmin + nz*sph_dz;
    // Update also box measures
    box_X = box_xmax - box_xmin;
    box_Y = box_ymax - box_ymin;
    box_Z = box_zmax - box_zmin;
    // This very important block for spherical coordinate boundary
    // conditions. Using it we fit the sizes dx, sph_dy and sph_dz
    // to initial box sizes, that is very important to spherical
    // and hybrid transformatins.
    if(Params::sph_BC_type == 1) {
        dx = box_X/nx;
        sph_dy = box_Y/ny;
        sph_dz = box_Z/nz;
        sph_dtheta = sph_box_theta/ny;
        sph_dphi = sph_box_phi/nz;
    }
    // Calculate box volume
    box_V = box_X * box_Y * box_Z;
    sph_box_V = 4/3*pi*(box_xmax*box_xmax*box_xmax - box_xmin*box_xmin*box_xmin);
    // Calculate tight box measures
    box_xmin_tight = box_xmin + box_eps * dx;
    box_xmax_tight = box_xmax - box_eps * dx;
    box_ymin_tight = box_ymin + box_eps * sph_dy;
    box_ymax_tight = box_ymax - box_eps * sph_dy;
    box_zmin_tight = box_zmin + box_eps * sph_dz;
    box_zmax_tight = box_zmax - box_eps * sph_dz;
    box_X_tight = box_xmax_tight - box_xmin_tight;
    box_Y_tight = box_ymax_tight - box_ymin_tight;
    box_Z_tight = box_zmax_tight - box_zmin_tight;
    // Prepare counters
    cnt_dt = 0;
    updateDependantParameters();
}

#endif

