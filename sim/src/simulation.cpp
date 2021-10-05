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
#include <iomanip>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <csignal>
#include <sstream>
#include "simulation.h"
#include "random.h"
#include "magneticfield.h"
#include "atmosphere.h"
#include "refinement.h"
#include "resistivity.h"
#include "forbidsplitjoin.h"
#include "backgroundcharge.h"
#include "detector.h"
#include "boundaries.h"
#include "vis/vis_data_source_simulation.h"
#include "vis/vis_db_factory.h"
#include "templates.h"
#include "chemistry.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

//! Simulation main log
Logger mainlog("simu.log", 100000);

//! Simulation error log
Logger errorlog("simu.err", 100000);

//! Simulation parameter log
Logger paramslog("params.log", 1000000, false);

//! Simulation grid. The most essential object in the program.
Tgrid g;

//! Simulation parameters
Params simuConfig;

//! Output file pointer
FILE *outputfp;

//! Constructor
Simulation::Simulation()
{
    timepool("Init");
    mainlog.init();
    errorlog.init();
    paramslog.init();
    // Start program execution time counter
    getExecutionSecs();
    MSGFUNCTIONCALL("Simulation::Simulation");
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    mainlog << "\nCARTESIAN COORDINATE SYSTEM\n";
#else
    mainlog << "\nSPHERICAL COORDINATE SYSTEM\n";
#ifdef WRITE_POPULATION_AVERAGES
    WARNINGMSG("WRITE_POPULATION_AVERAGES defined but it may not work in spherical coordinates.");
#endif
#ifdef PERIODIC_FIELDS_Y
#error PERIODIC_FIELDS_Y does not work with USE_SPHERICAL_COORDINATE_SYSTEM.
#endif
#ifdef RECONNECTION_GEOMETRY
#error RECONNECTION_GEOMETRY does not work with USE_SPHERICAL_COORDINATE_SYSTEM.
#endif
#endif
    mainlog << "\n";
    static bool onlyOneSimulationObject = false;
    if (onlyOneSimulationObject == true) {
        ERRORMSG("only one Simulation object can be defined");
        doabort();
    } else {
        onlyOneSimulationObject = true;
    }
    simuConfig.init();
    simuConfig.readAndInitVariables(Params::configFileName);
    mainlog << "HWA object ID = " << Params::objectIdHWA << " (";
    if(Params::objectIdHWA == 0) {
        mainlog << "Unspecified";
    } else if(Params::objectIdHWA == 1) {
        mainlog << "Mercury";
    } else if(Params::objectIdHWA == 2) {
        mainlog << "Venus";
    } else if(Params::objectIdHWA == 3) {
        mainlog << "Earth";
    } else if(Params::objectIdHWA == 4) {
        mainlog << "Moon";
    } else if(Params::objectIdHWA == 5) {
        mainlog << "Mars";
    } else if(Params::objectIdHWA == 6) {
        mainlog << "Asteroid";
    } else if(Params::objectIdHWA == 7) {
        mainlog << "Ganymede";
    } else if(Params::objectIdHWA == 8) {
        mainlog << "Titan";
    } else if(Params::objectIdHWA == 9) {
        mainlog << "Pluto";
    } else if(Params::objectIdHWA == 10) {
        mainlog << "Comet";
    } else if(Params::objectIdHWA == 11) {
        mainlog << "Heliosphere/solar wind";
    } else if(Params::objectIdHWA == 12) {
        mainlog << "Exoplanet";
    } else {
        mainlog << "N/A)\n";
        ERRORMSG("Bad HWA object ID");
        doabort();
    }
    mainlog << ")\n";
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    mainlog << "\nCREATING PARTICLE SPECTRA IN CELLS\n\n";
    mainlog << "|---------------- CELL SPECTRA ----------------|\n"
            << "| Emin = " << Params::spectraEmin_eV << " eV\n"
            << "| Emax = " << Params::spectraEmax_eV << " eV\n"
            << "| Nbins = " << Params::spectraNbins << "\n"
            << "| logBins = " << Params::spectraLogBins << "\n"
            << "| EminAllLowEnergies = " << Params::spectraEminAll << "\n"
            << "| EmaxAllHighEnergies = " << Params::spectraEmaxAll << "\n"
            << "| spectraMethod = " << Params::spectraMethod << "\n"
            << "| spectraUnit = " << Params::spectraUnit << "\n"
            << "|----------------------------------------------|\n";

#endif
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    g.init(Params::nx, Params::ny, Params::nz, Params::box_xmin, Params::box_ymin, Params::box_zmin, Params::dx);
#else
    g.sph_init(Params::nx, Params::ny, Params::nz, Params::box_xmin, Params::box_ymin, Params::box_zmin, Params::sph_theta_min, Params::sph_phi_min, Params::dx, Params::sph_dy, Params::sph_dz);
#endif
    initializeSimulation();
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    g.set_save_particles_orbit(Params::saveParticlesAlongOrbitFile.c_str());
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    if (Params::averaging == false) {
        ERRORMSG("cannot create cell spectra HC files because averaging is turned off");
        doabort();
    }
#endif
    MSGFUNCTIONEND("Simulation::Simulation");
}

//! Destructor
Simulation::~Simulation()
{
    MSGFUNCTIONCALL("Simulation::~Simulation");
    for (std::vector<VisDB*>::iterator visDB = visWriters.begin(); visDB != visWriters.end(); ++visDB) {
        delete *visDB;
    }
    delete visDataSourceImpl;
    MSGFUNCTIONEND("Simulation::~Simulation");
}

//! Kinda dummy function but needed
bool Simulation::AlwaysTrue(TLinkedParticle&)
{
    return true;
}

//! Particle dumping
bool Simulation::output(TLinkedParticle& part)
{
    fprintf(outputfp,"%g %g %g %g %g %g %g %d\n",
            double(part.x),double(part.y),double(part.z),
            double(part.vx),double(part.vy),double(part.vz),
            double(part.w),int(part.popid));
    return true;
}

//! Initialize the simulation
void Simulation::initializeSimulation()
{
    MSGFUNCTIONCALL("Simulation::initializeSimulation");
    macroParticlePropagations = 0.0;
    // Initialize our portable random number generator with some
    // seed (always the same ==> repeatable)
    portrand.init(1024);
    initializeGridRefinement();
    initializeForbidSplitJoin();
    initializeResistivity();
    initializeBackgroundChargeDensity();
    initializeMagneticField();
    initializeSplitJoin();
    mainlog << "|---------------- GENERAL SIMULATION INFORMATION ----------------|\n"
            << "| R_P = " << Params::R_P/1e3 << " km\n"
            << "| R_zeroFields (inner boundary) = " << Params::R_zeroFields/1e3 << " km = " << Params::R_zeroFields/Params::R_P << " R_P = " << (Params::R_zeroFields-Params::R_P)/1e3 << " km + R_P\n"
            << "| dt = " << Params::dt << " s\n"
            << "| t_max = " << Params::t_max << " s\n"
            << "| save interval = " << Params::saveInterval << " s\n"
            << "| input interval = " << Params::inputInterval << " s\n"
            << "| log interval = " << Params::logInterval << " s\n"
            << "| breakpoint interval (cyclic) = " << Params::wsDumpInterval[0] << " s\n"
            << "| breakpoint interval (keep) = " << Params::wsDumpInterval[1] << " s\n"
            << "| Average output files = " << (Params::averaging ? "yes" : "no") << "\n"
            << "| MacroParticlesPerCell = "  << Params::macroParticlesPerCell << "\n"
            << "| Maximun number of grid refinement levels = " << Params::maxGridRefinementLevel << "\n"
            << "| Current number of grid refinement levels = " << Params::currentGridRefinementLevel << "\n"
            << "| min_dx = " << Params::dx/real(1 << Params::currentGridRefinementLevel)*1e-3 << " km = " << Params::dx/(Params::R_P*real(1 << Params::currentGridRefinementLevel)) << " R_P\n"
            << "| max_dx = " << Params::dx*1e-3 << " km = " << Params::dx/Params::R_P << " R_P\n"
            << "| vi_max = " << Params::vi_max/1e3 << " km/s\n"
            << "| Ue_max = " << Params::Ue_max/1e3 << " km/s\n"
            << "| rho_q_min  = " << Params::rho_q_min/(Params::e*1e6) << " e/cm^3\n"
            << "| maxVw      = " << Params::maxVw/1e3 << " km/s\n"
            << "| min(dx)/dt = " << Params::dx/real(1 << Params::currentGridRefinementLevel)/Params::dt*1e-3 << " km/s\n"
            << "| M_P = " << Params::M_P << " kg\n"
            << "| Use Jstag = " << (Params::useJstag ? "yes" : "no") << "\n"
            << "| Use nodeUe = " << (Params::useNodeUe ? "yes" : "no") << "\n"
            << "| Use predictor-corrector = " << (Params::fieldPredCor ? "yes" : "no") << "\n"
            << "| Use gravity = " << (Params::useGravitationalAcceleration ? "yes" : "no") << "\n"
            << "| Use macroparticle splitting = " << (Params::useMacroParticleSplitting ? "yes" : "no") << "\n"
            << "| Use macroparticle joining = " << (Params::useMacroParticleJoining ? "yes" : "no") << "\n"
            << "| Split&join deviation allowed = " << Params::splitJoinDeviation[0] << "\n"
            << "| Split&join probability method = " << Params::splitJoinDeviation[1] << "\n"
            << "| Ecut = " << Params::Ecut << " V/m\n"
            << "| E filter = " << Params::electricFieldSmoothingNumber << "\n"
            << "| rho filter = " << Params::densitySmoothingNumber << "\n"
            << "|----------------------------------------------------------------|\n";
    mainlog
            << "|---------------- ELECTRON PRESSURE ----------------|\n"
            << "| Use electron pressure E term = " << (Params::electronPressure ? "yes" : "no") << "\n";
    if(Params::electronPressure == true) {
        mainlog
                << "| Te = " << Params::Te << " K (isothermal electrons)\n"
                << "| R_zeroPolarizationField = " << Params::R_zeroPolarizationField/1e3 << " km = " << Params::R_zeroPolarizationField/Params::R_P << " R_P\n";
    }
    mainlog << "|---------------------------------------------------|\n";
    if(Params::inputInterval < Params::dt && Params::inputInterval > 0) {
        ERRORMSG("inputInterval < dt");
        doabort();
    }
    if(Params::saveInterval < Params::dt && Params::saveInterval > 0) {
        ERRORMSG("saveInterval < dt");
        doabort();
    }
    if(Params::logInterval < Params::dt && Params::logInterval > 0) {
        ERRORMSG("logInterval < dt");
        doabort();
    }
    if((Params::wsDumpInterval[0] < Params::dt && Params::wsDumpInterval[0] > 0) || (Params::wsDumpInterval[1] < Params::dt && Params::wsDumpInterval[1] > 0)) {
        ERRORMSG("wsDumpInterval < dt");
        doabort();
    }
    mainlog
            << "|-------------------- IMF --------------------|\n"
            << "|   B  = [" << Params::SW_Bx/1e-9 << ", " << Params::SW_By/1e-9 << ", " << Params::SW_Bz/1e-9 << "] nT\n"
            << "|  |B| = " << norm(Params::SW_Bx,Params::SW_By,Params::SW_Bz)/1e-9 << " nT\n"
            << "|---------------------------------------------|\n";
    mainlog
            << "|----------- MAGNETIC FIELD BOUNDARY CONDITIONS -----------|\n"
            << "| Front wall (+X) = " << (Params::Bboundaries[0] ? "IMF" : "Neumann zero-gradient") << "\n"
            << "| Side wall (+-Y) = " << (Params::Bboundaries[1] ? "IMF" : "Neumann zero-gradient") << "\n"
            << "| Side wall (+-Z) = " << (Params::Bboundaries[2] ? "IMF" : "Neumann zero-gradient") << "\n"
            << "| Back wall (-X)  = " << (Params::Bboundaries[3] ? "IMF" : "Neumann zero-gradient") << "\n"
            << "| IMF norm. comp. = " << (Params::Bboundaries[4] ? "yes" : "no") << "\n"
            << "|----------------------------------------------------------|\n";
    // Initialize particle populations
    for(unsigned int i=0; i < Params::pops.size(); ++i) {
        Params::pops[i]->initialize();
    }
    // Write population log entries
    mainlog << "\n";
    for(unsigned int i=0; i < Params::pops.size(); ++i) {
        mainlog << "|--------------- PARTICLE POPULATION ---------------|\n";
        mainlog << Params::pops[i]->toString();
        mainlog << "|---------------------------------------------------|\n\n";
    }
    // Write detector log entries
    mainlog << "\n";
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        mainlog << "|--------------- DETECTOR SET ---------------|\n";
        mainlog << Params::detectors[i]->toString();
        mainlog << "|---------------------------------------------------|\n\n";
    }
    ParticleProcesses::writeLog();
#ifndef NO_DIAGNOSTICS
    // Initialize diagnostics
    Params::diag.init();
#endif
    // Show hc-file configurations
    vector<string> hcFilePrefix;
    vector< vector<int> > popId;
    Population::getHcFileConfigs(hcFilePrefix,popId);
    mainlog << "|--------------- HC-FILE CONFIGURATIONS ---------------|\n";
    for(unsigned int i = 0; i < hcFilePrefix.size(); ++i) {
        mainlog << "| "  << hcFilePrefix[i] << "*.hc: ";
        for(unsigned int j = 0; j < popId[i].size(); ++j) {
            mainlog << Params::pops[ popId[i][j] ]->getIdStr() << " ";
        }
        mainlog << "\n";
    }
    mainlog << "|------------------------------------------------------|\n\n";
#ifdef USE_PARTICLE_SUBSTEPPING
    mainlog << "|---------------- PARTICLE SUBCYCLING -----------------|\n\n";
    mainlog << "Substepping scheme: " << Params::subcycleType << "\n";
    mainlog << "Substepping maximum level: " << Params::subcycleMaxLevel << "\n";
    mainlog << "|------------------------------------------------------|\n\n";
#endif
    //! If user terminates with kill or ctrl-c
    signal(SIGTERM,&TermHandler);
    signal(SIGINT,&TermHandler);
    if (Params::saveExtraHcFiles == true) {
        saveExtraHcFiles();
    }
    // Start hc-file averaging
    if (Params::averaging == true) {
        g.begin_average();
    }
    // Initialize visualization data writer(s)
    visDataSourceImpl = new SimulationVisDataSourceImpl(simuConfig, g);
    if(Params::saveVTK == 1) {
        visWriters.push_back(&(VisDBFactory::getVisDB(VTK_LEGACY_BINARY, *visDataSourceImpl)));
    } else if(Params::saveVTK == 2) {
        visWriters.push_back(&(VisDBFactory::getVisDB(VTK_LEGACY_ASCII, *visDataSourceImpl)));
    }
    MSGFUNCTIONEND("Simulation::initializeSimulation");
}

//! Run simulation
void Simulation::run()
{
    Params::t = 0.0;
    // Read previous state, if it was given
    if (Params::wsFileGiven) {
        timepool("LoadBreakpoint");
        readState( Params::wsFileName.c_str() );
        finalizeTimestep(false);
    }
    mainlog << "TIMELOOP START\n";
    bool noTermination = true;
    while(Params::t <= Params::t_max && noTermination == true) {
        // Propagate simulation forward one timestep.
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        stepForward();
#else
        sph_stepForward();
#endif
        timepool("Misc");
        noTermination = finalizeTimestep();
    }
    mainlog << "TIMELOOP END\n";
}

//! Run simulation but construct only the input parameter logfile
void Simulation::runOnlyParameterDynamics()
{
    Params::t = 0.0;
    while(Params::t <= Params::t_max) {
        Params::cnt_dt++;
        Params::t += Params::dt;
        if (Params::inputInterval > 0 && (Params::cnt_dt % int(Params::inputInterval/Params::dt + 0.5) == 0) && Params::t > 0) {
            updateParams();
        }
    }
}

//! Forward simulation one timestep
void Simulation::stepForward()
{
    timepool("Newparticle");
    for (unsigned int i = 0; i < Params::pops.size(); ++i) {
        Params::pops[i]->createParticles();
    }
    timepool("Field");
    if(Params::propagateField == true) {
        g.zero_rhoq_nc_Vq();
    }
    timepool("Xpropag");
#ifndef USE_PARTICLE_SUBCYCLING
    g.particle_pass(&PropagateX);
#else
    g.particle_pass(&PropagatePart1);
#endif
    g.particle_pass_with_relocation(&AlwaysTrue);
    timepool("Field");
    if(Params::propagateField == true) {
        g.finalize_accum();
        g.smoothing();//smooth the particle related variables before propagating the field.
        if (Params::fieldPredCor == false) {
            fieldpropagate(Tgrid::FACEDATA_B, Tgrid::FACEDATA_B, Tgrid::FACEDATA_B,Params::dtField,true);
        } else {
            fieldpropagate(Tgrid::FACEDATA_BSTAR, Tgrid::FACEDATA_B, Tgrid::FACEDATA_B,Params::dtField/2,false);
            fieldpropagate(Tgrid::FACEDATA_B, Tgrid::FACEDATA_B, Tgrid::FACEDATA_BSTAR,Params::dtField,true);
        }
        if(Params::electronPressure==true) {
            g.FC(Tgrid::FACEDATA_B,Tgrid::CELLDATA_B);//update cell-data B to calculate E for particle acceleration.
            // g.Neumann_rhoq();//set up the boundary condition of rho_q for Cell-to-node interpolation. This is redundant because it is already taken care of in finalize_accum().
            g.CN_rhoq();//Interpolate rho_q from cell to node (saved temporarily in NODEDATA_UE[0])
            g.NF_rhoq();//Interpolate rho_q from node to face (saved temporarily in FACEDATA_MINUSDB)
            g.CalcGradient_rhoq();//Calculate the gradient of rho_q (saved temporarily in CELLDATA_TEMP1)
            g.calc_cell_E();//calculate E at the center of the cell for particle acceleration. (Saved temporarily in CELLDATA_TEMP2).
        }
    }
    timepool("Vpropag");
#ifndef USE_PARTICLE_SUBCYCLING
    g.particle_pass(&PropagateV);
#else
    g.particle_pass(&PropagatePart2);
#endif
    timepool("splitjoin");
    if(Params::useMacroParticleSplitting == true || Params::useMacroParticleJoining == true) {
        int nsplit, njoined;
        g.split_and_join(nsplit,njoined);
    }
    if(ParticleProcesses::isInitialized() == true) {
        timepool("ParticleProcesses");
        g.particle_pass(ParticleProcesses::run);
    }
    timepool("Misc");
    // Run field detectors
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        Params::detectors[i]->runFieldDetects();
        Params::detectors[i]->runTestParticles();
    }
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    if(Params::saveParticlesAlongOrbit == true) {
        g.particles_write();
    }
#endif
    // Count particle propagations
    macroParticlePropagations += g.Nparticles();
    // particle dumping
    static const int Ndumps = 4;
    static real particle_dump_start[Ndumps] = {1e20,1e20,1e20,1e20};
    static bool particle_dump_done[Ndumps] = {0};
    for(int i=0; i<Ndumps; ++i) {
        if ((Params::cnt_dt % int(particle_dump_start[i]/Params::dt+0.5) == 0) && (particle_dump_done[i] == false)  && (Params::cnt_dt > 0)) {
            string fn = "pdump_" + Params::getSimuTimeStr() + ".dat";
            outputfp = fopen(fn.c_str(),"w");
            //fprintf(outputfp,"D=2r %d %d\n",g.Nparticles(),8);
            fprintf(outputfp,"# x y z vx vy vz w popid\n");
            g.particle_pass(&output);
            fclose(outputfp);
            particle_dump_done[i] = true;
            mainlog << "Dumped particles in \"" << fn << "\" at t=" << Params::t << "\n";
        }
    }
}

//! Finalize time step
bool Simulation::finalizeTimestep(bool doBreakpointing)
{
    // Save step for hc-files
    if (Params::saveInterval > 0 && (Params::cnt_dt % int(Params::saveInterval/Params::dt+0.5) == 0)) {
        timepool("SaveStep");
        saveStep();
    }
    // Breakpointing
    if (doBreakpointing == true && Params::wsDumpInterval[0] > 0 && (Params::cnt_dt % int(Params::wsDumpInterval[0]/Params::dt+0.5) == 0) && Params::cnt_dt >0) {
        timepool("SaveBreakpoint");
        // single file name, older file overwritten
        dumpState("breakpoint.dat");
    }
    if (doBreakpointing == true && Params::wsDumpInterval[1] > 0 && (Params::cnt_dt % int(Params::wsDumpInterval[1]/Params::dt+0.5) == 0) && Params::cnt_dt > 0) {
        timepool("SaveBreakpoint");
        // unique file name, older file not overwritten
        string fn = "breakpoint_" + Params::getSimuTimeStr() + ".dat";
        dumpState(fn.c_str());
    }
    timepool("Misc");
#ifndef NO_DIAGNOSTICS
    // Logging
    if (Params::logInterval > 0 && (Params::cnt_dt % int(Params::logInterval/Params::dt+0.5) == 0)) {
        Params::diag.run();
    }
#endif
    // Check program termination flag
    if (Params::stoppingPhase == true) {
        mainlog << "STOPPING: program termination flag detected\n";
        saveVisualizationFiles();
        return false;
    }
    // Increase time
    Params::cnt_dt++;
    Params::t += Params::dt;
    // Update input parameters
    if (Params::inputInterval > 0 && (Params::cnt_dt % int(Params::inputInterval/Params::dt + 0.5) == 0) && Params::t > 0) {
        updateParams();
    }
    return true;
}

//! Update simulation input parameters
void Simulation::updateParams()
{
    simuConfig.readAndUpdateVariables(Params::configFileName);
    setResistivity();
}

//! Do save step 
void Simulation::saveStep()
{
    MSGFUNCTIONCALL("Simulation::saveStep");
    // Write visualization files
    saveVisualizationFiles();
    // Output running time counters
    time_t rawtime;
    time(&rawtime);
    mainlog << "|--------------- RUNNING TIME ---------------|\n";
    mainlog << "| Time now      : " << ctime(&rawtime);
    mainlog << "| From startup  : " << secsToTimeStr(getExecutionSecs()) << "\n"
            << "| Last savestep : " << secsToTimeStr(getLastIntervalSecs()) << "\n";
    mainlog << "|--------------------------------------------|\n";
    MSGFUNCTIONEND("Simulation::saveStep");
    // Zero logging counters
    mainlog.zeroAllCounters();
    errorlog.zeroAllCounters();
}

//! Save hc-files
void Simulation::saveVisualizationFiles()
{
    bool averageOk = false;
    if (Params::averaging == true) {
        averageOk = g.end_average();
    }
    if(Params::saveHC > 0) {
        // Get hc-file configurations
        vector<string> hcFilePrefix;
        vector< vector<int> > popId;
        Population::getHcFileConfigs(hcFilePrefix,popId);
        // Get hc-file format
        string hcFileFormat;
        if(Params::saveHC == 1) {
            hcFileFormat = "binary";
        } else if(Params::saveHC == 2) {
            hcFileFormat = "ascii";
        } else {
            hcFileFormat = "binary";
        }
        // Save hc-files
        for (unsigned int i = 0; i < hcFilePrefix.size(); i++) {
            // Do not save hcfile marked with string "-"
            if (hcFilePrefix[i].compare("-") == 0) {
                continue;
            }
            string fn = hcFilePrefix[i] + "_hybstate_" + Params::getSimuTimeStr() + ".hc";
            g.hcwrite_MHD(fn.c_str(),hcFileFormat,"populations",popId[i]);
        }
        // Save temporally averaged hc-file
        if (averageOk == true) {
            string fn = "average_hybstate_" + Params::getSimuTimeStr() + ".hc";
            g.hcwrite_MHD(fn.c_str(),hcFileFormat,"average",vector<int>());
#ifdef SAVE_POPULATION_AVERAGES
            for (unsigned int i = 0; i < hcFilePrefix.size(); i++) {
                // Do not save hcfile marked with string "-"
                if (hcFilePrefix[i].compare("-") == 0) {
                    continue;
                }
                string fn = hcFilePrefix[i] + "_ave_hybstate_" + Params::getSimuTimeStr() + ".hc";
                g.hcwrite_MHD(fn.c_str(),hcFileFormat,"populations_ave",popId[i]);
            }
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
            for (unsigned int i = 0; i < hcFilePrefix.size(); i++) {
                if (hcFilePrefix[i].compare("-") == 0) {
                    continue;
                }
                string fn = hcFilePrefix[i] + "_spectra_" + Params::getSimuTimeStr() + ".hc";
                g.hcwrite_SPECTRA(fn.c_str(),hcFileFormat,popId[i]);
            }
#endif
        }
        // Save plasma hc-file
        if (Params::plasma_hcfile == true) {
            string fn = "plasma_hybstate_" + Params::getSimuTimeStr() + ".hc";
            g.hcwrite_MHD(fn.c_str(),hcFileFormat,"plasma",vector<int>());
        }
        // Save debug hc-file
        if (Params::dbug_hcfile == true) {
            string fn = "DBUG_hybstate_" + Params::getSimuTimeStr() + ".hc";
            g.hcwrite_DBUG(fn.c_str());
        }
    } // if(saveHC)
    if(Params::saveVTK > 0) {
        string fn = "hybstate_" + Params::getSimuTimeStr();
        for (std::vector<VisDB*>::iterator visDB = visWriters.begin(); visDB != visWriters.end(); ++visDB) {
            //(*visDB)->writeVisValues(fn.c_str(),std::vector<DB::Particles>(1, DB::Particles()));
            (*visDB)->writeVisValues(fn.c_str());
        }
    }
    if (Params::averaging == true) {
        g.begin_average();
    }
}

//! Save extra hc-files
void Simulation::saveExtraHcFiles()
{
    for(unsigned int i=0; i < Params::pops.size(); ++i) {
        Params::pops[i]->writeExtraHcFile();
    }
}

//! Save a breakpoint
void Simulation::dumpState(const char *fileName)
{
    MSGFUNCTIONCALL("Simulation::dumpState");
    mainlog << "Creating a break point: " << fileName << "\n";
    ofstream dumpFile(fileName, ios::out);
    g.dumpState(dumpFile);
    dumpFile.close();
    MSGFUNCTIONEND("Simulation::dumpState");
}

//! Load a breakpoint
void Simulation::readState(const char *fileName)
{
    MSGFUNCTIONCALL("Simulation::readState");
    mainlog << "Loading a break point: " << fileName << "\n";
    ifstream is(fileName);
    g.readState(is);
    MSGFUNCTIONEND("Simulation::readState");
}

//! Finalize simulation main class
int Simulation::finalize()
{
    MSGFUNCTIONCALL("Simulation::finalize");
    const double cpu = timepool.cputime();
    mainlog << "|-------------------------------------------\n"
            << "| " << macroParticlePropagations << " macroparticles propagated in " << cpu << " seconds\n"
            << "| " << macroParticlePropagations/cpu << " macros/second\n"
            << "|-------------------------------------------\n";
    //portrand.save("portrand.state");
    MSGFUNCTIONEND("Simulation::finalize");
    return 0;
}

/** \brief Magnetic field at the box boundaries. Implements e.g. the IMF entering the box.
 *
 * dim = 0 => Set boundary B_y and B_z (don't touch B_x, since it is the normal component to the boundary)
 * dim = 1 => Set boundary B_x and B_z (don't touch B_y, since it is the normal component to the boundary)
 * dim = 2 => Set boundary B_x and B_y (don't touch B_z, since it is the normal component to the boundary)
 *
 * B_x: celldata[Tgrid::CELLDATA_B][0]
 * B_y: celldata[Tgrid::CELLDATA_B][1]
 * B_z: celldata[Tgrid::CELLDATA_B][2]
 *
 */
void Simulation::BoundaryB(datareal celldata[Tgrid::NCELLDATA][3], int dim)
{
    switch(dim) {
        // Walls in x-direction
    case 0:
        celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz;
        if (Params::Bboundaries[4]) { //if also for the perpendicular component
            celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        }
        break;
        // Walls in y-direction
    case 1:
        celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz;
        if (Params::Bboundaries[4]) {
            celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        }
        break;
        // Walls in z-direction
    case 2:
        celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        if (Params::Bboundaries[4]) {
            celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz;
        }
        break;
#ifdef RECONNECTION_GEOMETRY
        // Back wall in x-direction
    case 3:
        celldata[Tgrid::CELLDATA_B][1] = -Params::boundary_By;
        celldata[Tgrid::CELLDATA_B][2] = -Params::boundary_Bz;
        if (Params::Bboundaries[4]) { //if also for the perpendicular component
            celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        }
        break;
#endif
    default:
        ERRORMSG("Unknown dim in BoundaryB function");
        doabort();
        break;
    }
}

#ifdef USE_PARTICLE_SUBCYCLING

//! (SUBCYCLING)
bool Simulation::PropagatePart1(TLinkedParticle& part, ParticlePassArgs a)
{
    if(Params::pops[part.popid]->getPropagateV() == false) {
        return true;
    }
    bool part_kept = false;
    // Check Particle dtlevel here
    int popstep;
    popstep = Params::pops[part.popid]->subcycleSteps;
    real pv, dv;
    pv = sqrt(part.vx*part.vx + part.vy*part.vy + part.vz*part.vz);
    dv = a.size/Params::dt;
    if(Params::subcycleType == 1) {
        part.dtlevel = floor(pv/dv*popstep);
    } else if(Params::subcycleType == 2) {
        part.dtlevel = max((int)floor( log( pv/dv*popstep ) /log(2)),0);
    } else {
        ERRORMSG2("Invalid substepping scheme", Params::subcycleType);
        doabort();
    }
    if(part.dtlevel < 0) part.dtlevel = 0;
    if(part.dtlevel > Params::subcycleMaxLevel) {
        part.dtlevel = Params::subcycleMaxLevel;
    }
    real tol = 1e-5;
    if(abs(part.accumed - 1) > tol) {
        errorlog << "Particle pop: " << part.popid << " Particle dtlevel: " << part.dtlevel << " steps/cell:" << (int)(a.size/pv/Params::dt_psub[part.dtlevel]) << " weight accumed: " << part.accumed << "\n";
    }
    part.accumed = 0;
    part_kept = PropagateX(part);
    if(!part_kept) return part_kept;
    if(part.dtlevel >= 1) {
        PropagateV(part);
        if(Params::subcycleType == 1) {
            for(int i = 1; i <= part.dtlevel/2; i++) {
                part_kept = PropagateX(part);
                if(!part_kept) return part_kept;
                PropagateV(part);
            }
        } else if(Params::subcycleType == 2) {
            for(int lev = 2; lev <= part.dtlevel; lev++)
                for(int i = 0; i < pow(2,lev-2); i++) {
                    part_kept = PropagateX(part);
                    if(!part_kept) return part_kept;
                    PropagateV(part);
                }
        }
    }
    return part_kept;
}

//! (SUBCYCLING)
bool Simulation::PropagatePart2(TLinkedParticle& part)
{
    if(Params::pops[part.popid]->getPropagateV() == false) {
        return true;
    }
    bool part_kept = true;
    if(part.dtlevel >= 1) {
        part_kept = PropagateX(part);
        if(!part_kept) return part_kept;
        if(Params::subcycleType == 1) {
            for(int i = 1; i <= (part.dtlevel-(1-part.dtlevel%2))/2 ; i++) {
                PropagateV(part);
                part_kept = PropagateX(part);
                if(!part_kept) return part_kept;
            }
        } else if(Params::subcycleType == 2) {
            for(unsigned char lev = 2; lev <= part.dtlevel; lev++)
                for(unsigned char i = 0; i < pow(2,lev-2); i++) {
                    PropagateV(part);
                    part_kept = PropagateX(part);
                    if(!part_kept) return part_kept;
                }
        }
    }
    PropagateV(part);
    return part_kept;
}

#endif

//! Move particle (r = v*dt)
bool Simulation::PropagateX(TLinkedParticle& part)
{
#ifndef USE_PARTICLE_SUBCYCLING
    static real pdt, pw;
    pdt = Params::dt;
    pw = 1;
#else
    real pdt, pw;
    pdt = Params::dt_psub[part.dtlevel];
    pw = Params::accum_psubfactor[part.dtlevel];
    part.accumed += Params::accum_psubfactor[part.dtlevel];
#endif
    fastreal r_new[3], rave[3];
    r_new[0] = part.x + part.vx*pdt;
    r_new[1] = part.y + part.vy*pdt;
    r_new[2] = part.z + part.vz*pdt;
    // Check particle detectors
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        Params::detectors[i]->runPartDetects(&part,r_new);
    }
    part.x = r_new[0];
    part.y = r_new[1];
    part.z = r_new[2];
    // Check boundary conditions
    const bool part_kept = Params::pops[part.popid]->checkBoundaries(part,rave);
    if (part_kept) {
        const fastreal v[3] = {part.vx, part.vy, part.vz};
        if(Params::propagateField == true) {
            g.accumulate_PIC(rave, v, part.w*pw, part.popid);
        }
    }
    return part_kept;
}


//! Accelerate particle (Lorentz force)
bool Simulation::PropagateV(TLinkedParticle& part)
{
    if(Params::pops[part.popid]->getPropagateV() == false) {
        return true;
    }
#ifndef USE_PARTICLE_SUBCYCLING
    static real pdt;
    pdt = Params::dt;
#else
    real pdt;
    pdt = Params::dt_psub[part.dtlevel];
#endif
    // Particle's centroid coordinates and velocity vectors
    const fastreal r[3] = {part.x, part.y, part.z};
    fastreal v[3] = {part.vx, part.vy, part.vz};
    real B[3],Ue[3];
    // Self-consistent B1 field from cell faces + constant B0 field => B(r) = B1(r) + B0(r)
    g.faceintpol(r, Tgrid::FACEDATA_B, B);
    addConstantMagneticField(r, B);
    // Velocity field of the electron fluid from cells
    // Do not pass r => uses saved_cellptr and avoids findcell call
    g.cellintpol(Tgrid::CELLDATA_UE,Ue);
    if(Params::electronPressure==true) {
        real Efield[3],tx,ty,tz,sx,sy,sz,dvx,dvy,dvz,vmx,vmy,vmz,v0x,v0y,v0z,vpx,vpy,vpz,qmideltT2,t2,b2;
        // get the NGP electric field
        // Do not pass r => uses saved_cellptr and avoids findcell call
        g.cellintpol(Tgrid::CELLDATA_TEMP2,Efield);
        Efield[0] += B[1]*Ue[2] - B[2]*Ue[1];
        Efield[1] += B[2]*Ue[0] - B[0]*Ue[2];
        Efield[2] += B[0]*Ue[1] - B[1]*Ue[0];
        qmideltT2= 0.5*Params::pops[part.popid]->q*pdt/Params::pops[part.popid]->m;
        dvx=qmideltT2*Efield[0];
        dvy=qmideltT2*Efield[1];
        dvz=qmideltT2*Efield[2];
        tx=qmideltT2*B[0];
        ty=qmideltT2*B[1];
        tz=qmideltT2*B[2];
        t2=tx*tx+ty*ty+tz*tz;
        b2=2./(1.+t2);
        sx=b2*tx;
        sy=b2*ty;
        sz=b2*tz;
        vmx=v[0]+dvx;
        vmy=v[1]+dvy;
        vmz=v[2]+dvz;
        v0x=vmx+vmy*tz-vmz*ty;
        v0y=vmy+vmz*tx-vmx*tz;
        v0z=vmz+vmx*ty-vmy*tx;
        vpx=vmx+v0y*sz-v0z*sy;
        vpy=vmy+v0z*sx-v0x*sz;
        vpz=vmz+v0x*sy-v0y*sx;
        v[0]=vpx+dvx;
        v[1]=vpy+dvy;
        v[2]=vpz+dvz;
    } else {
        // Vector: dU = v_i - U_e
        real dU[3] = { v[0]-Ue[0], v[1]-Ue[1], v[2]-Ue[2] };
        // Constant: alpha/2 = q*dt/(2*m)
        const real half_alpha = 0.5*Params::pops[part.popid]->q*pdt/Params::pops[part.popid]->m;
        // Vector: W = q*dt*B/(2*m)
        real b[3] = {half_alpha*B[0], half_alpha*B[1], half_alpha*B[2]};
        // |W|^2
        const real b2 = sqr(b[0]) + sqr(b[1]) + sqr(b[2]);
        // Constant: beta = 2/(1+|b|^2)
        const real beta = 2.0/(1.0 + b2);
        // Cross products
        real dUxb[3],dUxbxb[3];
        crossProduct(dU,b,dUxb);
        crossProduct(dUxb,b,dUxbxb);
        // Add velocity components
        v[0] += beta*( dUxb[0] + dUxbxb[0] );
        v[1] += beta*( dUxb[1] + dUxbxb[1] );
        v[2] += beta*( dUxb[2] + dUxbxb[2] );
    }
    // Gravity correction
    if(Params::useGravitationalAcceleration == true) {
        real rLength = sqrt( sqr(part.x) + sqr(part.y) + sqr(part.z) );
#ifndef USE_PARTICLE_SUBCYCLING
        real s = -Params::GMdt/cube(rLength);
#else
        real s = -Params::GMdt/cube(rLength)*Params::dt/pdt;
#endif
        v[0] += s*part.x;
        v[1] += s*part.y;
        v[2] += s*part.z;
    }
    const real v2 = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    // Check particle maximum speed (CONSTRAINT)
    if (v2 > Params::vi_max2) {
        const real norm = Params::vi_max/sqrt(v2);
        v[0] *= norm;
        v[1] *= norm;
        v[2] *= norm;
#ifndef NO_DIAGNOSTICS
        // Increase particle speed cutting rate counter
        Params::diag.pCounter[part.popid]->cutRateV += 1.0;
#endif
    }
    part.vx = v[0];
    part.vy = v[1];
    part.vz = v[2];
    return true;
}

//! Propagate magnetic field (Faraday's law)
void Simulation::fieldpropagate(Tgrid::TFaceDataSelect fsBnew,
                                Tgrid::TFaceDataSelect fsBold, Tgrid::TFaceDataSelect fsBrhs,
                                real fp_dt, bool do_upwinding)
{
    // 1. Calculate the current density based on Brhs
    // Interpolate B to cell quantity from faces (interior cells only)
    g.FC(fsBrhs,Tgrid::CELLDATA_B);
    // Extend B to ghostcells by applying homogeneous Neumann boundary conditions
    g.Neumann(Tgrid::CELLDATA_B);
    // Set user given boundary B fields
    if (Params::Bboundaries[0] == true) { // Front x-boundary
        g.boundarypass(0,1,BoundaryB); // Frontwall
    }
    if (Params::Bboundaries[1] == true) { // y-boundaries
        g.boundarypass(1,0,BoundaryB);
        g.boundarypass(1,1,BoundaryB);
    }
    if (Params::Bboundaries[2] == true) { // z-boundaries
        g.boundarypass(2,0,BoundaryB);
        g.boundarypass(2,1,BoundaryB);
    }
    if (Params::Bboundaries[3] == true) { // Back x-boundary
        g.boundarypass(0,0,BoundaryB); // Backwall
    }
    // Interpolate B from cells to nodes
    g.CN(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B);
    if(Params::useNodeUe == true) {
       g.Neumann_rhoq();  // should allready be done?
       g.CN_ne();
    }
    if(Params::useJstag == true) { // alternative j calculation
        // J TO NODES CALCULATION.
        g.boundary_faces(Tgrid::FACEDATA_B);
        g.calc_node_j();
        g.NC(Tgrid::NODEDATA_J,Tgrid::CELLDATA_J);
    } else {
        // Calculate J in cell faces by taking curl(B)/mu_0 from the nodes
        g.FaceCurl(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
        // Interpolate J from cell faces to cells (interior cells only)
        g.FC(Tgrid::FACEDATA_J,Tgrid::CELLDATA_J);
    }
    // 2. ue = u - j/(e*n)
    // Calculate U_e in cells (CELLDATA_UE)
    // (note: This must be done even if using useNodeUe, because velocity propagation needs
    // cell value and interpolation will not do.)
    g.calc_ue();
    // Extend U_e to ghost cells
    g.Neumann(Tgrid::CELLDATA_UE);
    if(Params::useNodeUe == true) { // Calculate Ue separately to nodes
        //g.Neumann_rhoq();  // should allready be done?
        //g.CN_ne();
        g.Neumann(Tgrid::CELLDATA_Ji);
        g.CN(Tgrid::CELLDATA_Ji,Tgrid::NODEDATA_Ji);
        g.calc_node_ue();
    } else { // no separate node Ue calculation
        // Interpolate U_e from cells to nodes
        g.CN(Tgrid::CELLDATA_UE,Tgrid::NODEDATA_UE);
    }
    // Upwind B. Approximate node_B in upstream of the node_UE. Displacement = mindx/2,
    // where mindx is the smallest dx of the cells touching the node. dt is not used.
    if (do_upwinding == true) {
        g.CN_donor(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B, Tgrid::NODEDATA_UE,fp_dt);
    }
    // Interpolate J from cells to nodes, for the same reason (but dont overwrite if using jstag).
    if(Params::useJstag == false) {
        g.CN(Tgrid::CELLDATA_J,Tgrid::NODEDATA_J);
    }
    // Calculate E = -U_e x (B + B_0) + eta*J at the nodes
    g.calc_node_E();
    g.smoothing_E();//smooth the electric field before propagating B field.
    // Calculate -dB/fp_dt in cell faces by taking curl(E) from the nodes
    g.FaceCurl(Tgrid::NODEDATA_E,Tgrid::FACEDATA_MINUSDB,1);
    // B_new = B_old - curl(E)*fp_dt
    g.FacePropagate(fsBold,fsBnew,fp_dt);
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "stepForward"
void Simulation::sph_stepForward()
{
    timepool("Newparticle");
    for(unsigned int i = 0; i < Params::pops.size(); ++i) {
        Params::pops[i]->createParticles();
    }
    timepool("Field");
    if(Params::propagateField == true) {
        g.zero_rhoq_nc_Vq();
    }
    timepool("Xpropag");
    g.particle_pass(&sph_PropagateX);
    g.particle_pass_with_relocation(&AlwaysTrue);
    timepool("Field");
    if(Params::propagateField == true) {
        g.sph_finalize_accum();
        //smooth the particle related variables before propagating the field.
        g.sph_smoothing();
        if(Params::fieldPredCor == true) {
            sph_fieldpropagate(Tgrid::FACEDATA_BSTAR, Tgrid::FACEDATA_B, Tgrid::FACEDATA_B,Params::dtField/2,false);
            sph_fieldpropagate(Tgrid::FACEDATA_B, Tgrid::FACEDATA_B, Tgrid::FACEDATA_BSTAR,Params::dtField,true);
        } else {
            sph_fieldpropagate(Tgrid::FACEDATA_B, Tgrid::FACEDATA_B, Tgrid::FACEDATA_B,Params::dtField,true);
        }
        if(Params::electronPressure == true) {
            if(Params::sph_BC_use_ghost_cell == 0) {
                //update cell-data B to calculate E for particle acceleration.
                g.sph_FC(Tgrid::FACEDATA_B,Tgrid::CELLDATA_B);
                //set up the boundary condition of rho_q for Cell-to-node interpolation. This is redundant because it is already taken care of in finalize_accum().
                g.sph_Neumann_rhoq_0();
                //Interpolate rho_q from cell to node (saved temporarily in NODEDATA_UE[0])
                g.sph_CNb_rhoq();
                //Interpolate rho_q from node to face (saved temporarily in FACEDATA_MINUSDB)
                g.sph_NF_rhoq();
                //Calculate the gradient of rho_q (saved temporarily in CELLDATA_TEMP1)
                g.sph_CalcGradient_rhoq();
                //calculate E at the center of the cell for particle acceleration. (Saved temporarily in CELLDATA_TEMP2).
                g.sph_calc_cell_E();
            } else if(Params::sph_BC_use_ghost_cell == 1) {
                //update cell-data B to calculate E for particle acceleration.
                g.sph_FC(Tgrid::FACEDATA_B,Tgrid::CELLDATA_B);
                //set up the boundary condition of rho_q for Cell-to-node interpolation. This is redundant because it is already taken care of in finalize_accum().
                g.sph_Neumann_rhoq();
                //Interpolate rho_q from cell to node (saved temporarily in NODEDATA_UE[0])
                g.sph_CN_rhoq();
                //Interpolate rho_q from node to face (saved temporarily in FACEDATA_MINUSDB)
                g.sph_NF_rhoq();
                //Calculate the gradient of rho_q (saved temporarily in CELLDATA_TEMP1)
                g.sph_CalcGradient_rhoq();
                //calculate E at the center of the cell for particle acceleration. (Saved temporarily in CELLDATA_TEMP2).
                g.sph_calc_cell_E();
            }
        }
    }
    timepool("Vpropag");
    g.particle_pass(&sph_PropagateV);
    timepool("splitjoin");
    if(Params::useMacroParticleSplitting == true || Params::useMacroParticleJoining == true) {
        int nsplit, njoined;
        g.split_and_join(nsplit,njoined);
    }
    timepool("Misc");
    // Run field detectors
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        Params::detectors[i]->runFieldDetects();
        Params::detectors[i]->runTestParticles();
    }
    // Count particle propagations
    macroParticlePropagations += g.Nparticles();
    /*int particle_dump_start[1] = {30000};
    if ((Params::cnt_dt % int(particle_dump_start[0]/Params::dt+0.5) == 0) && particle_dump1_done == 0 && Params::cnt_dt > 0) {
    	string fn = "pdump_" + Params::getSimuTimeStr() + ".dat";
    	outputfp = fopen(fn.c_str(),"w");
    	fprintf(outputfp,"D=2r %d %d\n",g.Nparticles(),8);
    	fprintf(outputfp,"# x y z vx vy vz w m\n");
    	g.particle_pass(&output);
    	fclose(outputfp);
    	particle_dump1_done = 1;
    	mainlog << "Dumped particles in \"" << fn << "\" at t=" << Params::t << "\n";
    }*/
}

//! (SPHERICAL) Toroidal boundary conditions
void Simulation::sph_BoundaryB(datareal celldata[Tgrid::NCELLDATA][3], gridreal sph_centroid[3], int dim)
{
    switch(dim) {
        // Walls in x-direction
    case 0:
        // Magnetic field is toroidal and depends of theta coordinate as sin(theta)
        celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz*pow(sin(sph_centroid[1]), 1);
        if(Params::Bboundaries[4]) { //if also for the perpendicular component
            celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        }
        break;
        // Walls in y-direction
    case 1:
        // Magnetic field is toroidal and depends of theta coordinate as sin(theta)
        celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz;
        if (Params::Bboundaries[4]) {
            celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        }
        break;
        // Walls in z-direction
    case 2:
        celldata[Tgrid::CELLDATA_B][0] = Params::boundary_Bx;
        celldata[Tgrid::CELLDATA_B][1] = Params::boundary_By;
        if (Params::Bboundaries[4]) {
            celldata[Tgrid::CELLDATA_B][2] = Params::boundary_Bz;
        }
        break;
    }
}

/** \brief (SPHERICAL) Spherical version of "PropagateX"
 *
 * Here we have propagateX in Hybrid coordinates. r = [r, theta_h, phi_h] and v = [r_dot, theta_h_dot, phi_h_dot].
 * In this case we don't need to change it
 */
bool Simulation::sph_PropagateX(TLinkedParticle& part)
{
    fastreal r_new[3], rave[3];
    fastreal r[3] = {part.x, part.y, part.z};
    fastreal v[3];
    real dt = Params::dt;
    // We take velocities and positions of the particle from sph_addparticle in Grid.cpp in hybrid coordinates.
    sph_transf_H2S_R(r);
    sph_transf_S2C_R(r);
    r_new[0] = r[0] + part.vx*dt;
    r_new[1] = r[1] + part.vy*dt;
    r_new[2] = r[2] + part.vz*dt;
    sph_transf_C2S_r(r_new);
    sph_transf_S2H_R(r_new);
    // Check particle detectors
    for(unsigned int i=0; i < Params::detectors.size(); ++i) {
        Params::detectors[i]->runPartDetects(&part,r_new);
    }
    bool part_kept = false;
    part.x = r_new[0];
    part.y = r_new[1];
    part.z = r_new[2];
    //Check boundary conditions
    part_kept = Params::pops[part.popid]->checkBoundaries(part,rave);
    if(part_kept) {
        v[0] = part.vx;
        v[1] = part.vy;
        v[2] = part.vz;
        if(Params::propagateField == true) {
            g.sph_accumulate_PIC(rave, v, part.w, part.popid);
        }
    }
    return part_kept;
}

/** (SPHERICAL) Spherical version of "PropagateV"
 *
 * NOTE: function copied to detect.h/cpp as void ParticleDetect::PropagateV(void)
 * This so that the test particle _mass_ and _charge_ can be directly used there
 * Copy changes made in this function to detect.cpp!
 */
bool Simulation::sph_PropagateV(TLinkedParticle& part)
{
    if(Params::pops[part.popid]->getPropagateV() == false) {
        return true;
    }
    // Particle's centroid coordinates and velocity vectors
    fastreal r[3] = {part.x, part.y, part.z};
    fastreal v[3] = {part.vx, part.vy, part.vz};
    real B[3],Ue[3];
    real dt = Params::dt;
    fastreal v_new[3];
    // Self-consistent B1 field from cell faces + constant B0 field => B(r) = B1(r) + B0(r)
    g.faceintpol(r, Tgrid::FACEDATA_B, B);
    //!g.sph_faceintpol(r, Tgrid::FACEDATA_B, B);
    addConstantMagneticField(r, B);
    // Velocity field of the electron fluid from cells
    // Do not pass r => uses saved_cellptr and avoids findcell call
    //!g.cellintpol(Tgrid::CELLDATA_UE,Ue);
    g.faceintpol(r, Tgrid::FACEDATA_UE, Ue);
    //g.sph_faceintpol(r, Tgrid::FACEDATA_UE, Ue); // Here we use first order of approximation of Ue insread of zero order
    // We take velocities and positions of the particle from sph_addparticle in Grid.cpp in hybrid coordinates.
    sph_transf_H2S_R(r);
    sph_transf_S2C_A(r,B);
    sph_transf_S2C_A(r,Ue);
    if(Params::electronPressure == true) {
        // New BB algorithm version
        real Efield[3] = {0 ,0, 0};
        real alpha, b2;
        real v1[3], v2[3], v3[3];
        g.cellintpol(Tgrid::CELLDATA_TEMP2,Efield);
        sph_transf_S2C_A(r,Efield);
        Efield[0] += B[1]*Ue[2] - B[2]*Ue[1];
        Efield[1] += B[2]*Ue[0] - B[0]*Ue[2];
        Efield[2] += B[0]*Ue[1] - B[1]*Ue[0];
        alpha  = 0.5*Params::pops[part.popid]->q*Params::dt/Params::pops[part.popid]->m;
        b2     = 2.0/(1.0 + sqr(alpha*B[0]) + sqr(alpha*B[1]) + sqr(alpha*B[2]));
        // Hockney&Eastwood p. 113 (4-99)
        v1[0] = v[0] + alpha*Efield[0];
        v1[1] = v[1] + alpha*Efield[1];
        v1[2] = v[2] + alpha*Efield[2];
        // p. 113 (4-103)
        v3[0] = v1[0] + alpha*(v1[1]*B[2] - v1[2]*B[1]);
        v3[1] = v1[1] + alpha*(v1[2]*B[0] - v1[0]*B[2]);
        v3[2] = v1[2] + alpha*(v1[0]*B[1] - v1[1]*B[0]);
        // p. 113 (4-102)
        v2[0] = v1[0] + b2*alpha*(v3[1]*B[2] - v3[2]*B[1]);
        v2[1] = v1[1] + b2*alpha*(v3[2]*B[0] - v3[0]*B[2]);
        v2[2] = v1[2] + b2*alpha*(v3[0]*B[1] - v3[1]*B[0]);
        // p. 101 (4-99)
        v[0] = v2[0] + alpha*Efield[0];
        v[1] = v2[1] + alpha*Efield[1];
        v[2] = v2[2] + alpha*Efield[2];
        // Old BB algorithm version
        /*real Efield[3],tx,ty,tz,sx,sy,sz,dvx,dvy,dvz,vmx,vmy,vmz,v0x,v0y,v0z,vpx,vpy,vpz,qmideltT2,t2,b2;
        // get the NGP electric field
        // Do not pass r => uses saved_cellptr and avoids findcell call
        g.cellintpol(Tgrid::CELLDATA_TEMP2,Efield);
        sph_transf_S2C_A(r,Efield); // Tranformation from spherical coordinates to Cartesian
        Efield[0] += B[1]*Ue[2] - B[2]*Ue[1];
        Efield[1] += B[2]*Ue[0] - B[0]*Ue[2];
        Efield[2] += B[0]*Ue[1] - B[1]*Ue[0];
        qmideltT2= 0.5*Params::pops[part.popid]->q*Params::dt/Params::pops[part.popid]->m;
        dvx=qmideltT2*Efield[0]; dvy=qmideltT2*Efield[1]; dvz=qmideltT2*Efield[2];
        tx=qmideltT2*B[0]; ty=qmideltT2*B[1]; tz=qmideltT2*B[2];
        t2=tx*tx+ty*ty+tz*tz;
        b2=2./(1.+t2);
        sx=b2*tx; sy=b2*ty; sz=b2*tz;
        vmx=v[0]+dvx; vmy=v[1]+dvy; vmz=v[2]+dvz;
        v0x=vmx+vmy*tz-vmz*ty; v0y=vmy+vmz*tx-vmx*tz; v0z=vmz+vmx*ty-vmy*tx;
        vpx=vmx+v0y*sz-v0z*sy; vpy=vmy+v0z*sx-v0x*sz; vpz=vmz+v0x*sy-v0y*sx;
        v[0]=vpx+dvx;
        v[1]=vpy+dvy;
        v[2]=vpz+dvz;*/
    } else {
        // Vector: dU = v_i - U_e
        real dU[3] = { v[0]-Ue[0], v[1]-Ue[1], v[2]-Ue[2] };
        // Constant: alpha/2 = q*dt/(2*m)
        const real half_alpha = 0.5*Params::pops[part.popid]->q*Params::dt/Params::pops[part.popid]->m;
        // Vector: W = q*dt*B/(2*m)
        real b[3] = {half_alpha*B[0], half_alpha*B[1], half_alpha*B[2]};
        // |W|^2
        const real b2 = sqr(b[0]) + sqr(b[1]) + sqr(b[2]);
        // Constant: beta = 2/(1+|b|^2)
        const real beta = 2.0/(1.0 + b2);
        // Cross products
        real dUxb[3],dUxbxb[3];
        crossProduct(dU,b,dUxb);
        crossProduct(dUxb,b,dUxbxb);
        // Add velocity components
        v[0] += beta*( dUxb[0] + dUxbxb[0] );
        v[1] += beta*( dUxb[1] + dUxbxb[1] );
        v[2] += beta*( dUxb[2] + dUxbxb[2] );
    }
    // Gravity correction
    // Velocities are calculated in Cartesian coordinates
    if (Params::useGravitationalAcceleration == true) {
        sph_transf_S2C_R(r);
        real rLength = sqrt(sqr(r[0]) + sqr(r[1]) + sqr(r[2]));
        real s = -Params::GMdt/cube(rLength);
        v[0] += s*r[0];
        v[1] += s*r[1];
        v[2] += s*r[2];
        // In spherical coordinates we need correst only r component of velocities, theta and phi component is perpendicular to gravitational field.
        //real rLength = abs(part.x);
        //real s = -Params::GMdt/cube(rLength);
        //v[0] += s*part.x;
    }
    const real v2 = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    // Check particle maximum speed (CONSTRAINT)
    if (v2 > Params::vi_max2) {
        const real norm = Params::vi_max/sqrt(v2);
        v[0] *= norm;
        v[1] *= norm;
        v[2] *= norm;
#ifndef NO_DIAGNOSTICS
        // Increase particle speed cutting rate counter
        Params::diag.pCounter[part.popid]->cutRateV += 1.0;
#endif
    }
    part.vx = v[0];
    part.vy = v[1];
    part.vz = v[2];
    return true;
}

//! (SPHERICAL) Spherical version of "fieldpropagate". Propagates magnetic field (Faraday's law).
void Simulation::sph_fieldpropagate(Tgrid::TFaceDataSelect fsBnew,
                                    Tgrid::TFaceDataSelect fsBold, Tgrid::TFaceDataSelect fsBrhs,
                                    real fp_dt, bool do_upwinding)
{
    if(Params::sph_propagation_type == 0) { // Flat front propagation
        if(Params::sph_BC_use_ghost_cell == 0) {
            g.sph_FC(fsBrhs,Tgrid::CELLDATA_B);
            //!g.sph_Cell_BC(Tgrid::CELLDATA_B); // Magnetic boundary conditions
            //!g.sph_YCave(Tgrid::CELLDATA_B);   // Theta Boundary conditions for cells (averaging)
            // Calculate B in cells (CELLDATA_UE)
            g.sph_Neumann_0(Tgrid::CELLDATA_B);
            g.sph_CNb(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B);
            g.sph_NC_copy(Tgrid::NODEDATA_B, Tgrid::CELLDATA_B_TEMP2);
            // Calculate J in cell faces by taking curl(B)/mu_0 from the nodes
            //!g.sph_FaceCurl(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
            g.sph_FaceCurl_2(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
            // Interpolate J from cell faces to cells (interior cells only)
            g.sph_FC(Tgrid::FACEDATA_J,Tgrid::CELLDATA_J);
            //g.FC_copy(Tgrid::FACEDATA_J, Tgrid::CELLDATA_J_TEMP);
            // Calculate U_e in cells (CELLDATA_UE)
            g.sph_calc_ue();
            //g.sph_YCave(Tgrid::CELLDATA_UE); // Theta Boundary conditions for cells (averaging)
            g.sph_Neumann_0(Tgrid::CELLDATA_UE);
            g.sph_CNb(Tgrid::CELLDATA_UE,Tgrid::NODEDATA_UE);
            // Ue - Spherical coordinates
            // Interpolation from nodes to faces (To calculate first order of accuracy of Lorents forse)
            g.sph_NF(Tgrid::NODEDATA_UE,Tgrid::FACEDATA_UE);
            g.sph_Neumann_0(Tgrid::CELLDATA_B);
            if(do_upwinding) {
                g.sph_CNb_donor(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B, Tgrid::NODEDATA_UE,fp_dt);
            }
            //g.sph_NC_copy(Tgrid::NODEDATA_B, Tgrid::CELLDATA_B_TEMP2);
            g.sph_Node_BC(Tgrid::NODEDATA_B); // Magnetic boundary conditions
            //!g.sph_YNave(Tgrid::NODEDATA_B); // Theta Boundary conditions for nodes (averaging)
            g.sph_NC_copy(Tgrid::NODEDATA_B, Tgrid::CELLDATA_B_TEMP2);
            // Interpolate J from cells to nodes, for the same reason
            g.sph_Neumann_0(Tgrid::CELLDATA_J);
            g.sph_CNb(Tgrid::CELLDATA_J,Tgrid::NODEDATA_J);
            g.sph_calc_node_E();
            //g.sph_YNave(Tgrid::NODEDATA_E); // Theta Boundary conditions for nodes (averaging)
            //!g.sph_NC_copy(Tgrid::NODEDATA_E, Tgrid::CELLDATA_E);
            // E - Spherical coordiantes
            //smooth the electric field before propagating B field.
            g.sph_smoothing_E();
            g.sph_NC_copy(Tgrid::NODEDATA_E, Tgrid::CELLDATA_E);
            g.sph_FaceCurl_2(Tgrid::NODEDATA_E,Tgrid::FACEDATA_MINUSDB,1);
            // B_new = B_old - curl(E)*fp_dt
            g.FacePropagate(fsBold,fsBnew,fp_dt);
        } //end of BC_use_ghost_cell == 0
        else if(Params::sph_BC_use_ghost_cell == 1) {
            // Interpolate B to cell quantity from faces (interior cells only)
            g.sph_FC(fsBrhs,Tgrid::CELLDATA_B);
            // Extend B to ghostcells by applying homogeneous Neumann boundary conditions
            g.sph_Neumann(Tgrid::CELLDATA_B);
            // Set user given boundary B fields
            if(Params::Bboundaries[0]) { // Front x-boundary
                g.sph_boundarypass(0,1,sph_BoundaryB); // Frontwall
            }
            if(Params::Bboundaries[1]) { // y-boundaries
                g.sph_boundarypass(1,0,sph_BoundaryB);
                g.sph_boundarypass(1,1,sph_BoundaryB);
            }
            if(Params::Bboundaries[2]) { // z-boundaries
                g.sph_boundarypass(2,0,sph_BoundaryB);
                g.sph_boundarypass(2,1,sph_BoundaryB);
            }
            if(Params::Bboundaries[3]) { // Back x-boundary
                g.sph_boundarypass(0,0,sph_BoundaryB); // Backwall
            }
            // Interpolate B from cells to nodes
            g.sph_CN(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B);
            g.sph_YNave(Tgrid::NODEDATA_B); // Theta Boundary conditions for nodes (averaging)
            // Calculate J in cell faces by taking curl(B)/mu_0 from the nodes
            g.sph_FaceCurl_2(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
            // Interpolate J from cell faces to cells (interior cells only)
            g.sph_FC(Tgrid::FACEDATA_J,Tgrid::CELLDATA_J);
            // 2. ue = u - j/(e*n)
            // Calculate U_e in cells (CELLDATA_UE)
            g.sph_calc_ue();
            // Extend U_e to ghost cells
            g.sph_Neumann(Tgrid::CELLDATA_UE);
            // Interpolate U_e from cells to nodes
            g.sph_CN(Tgrid::CELLDATA_UE,Tgrid::NODEDATA_UE);
            //g.sph_YNave(Tgrid::NODEDATA_UE); // Theta Boundary conditions for nodes (averaging)
            g.sph_NF(Tgrid::NODEDATA_UE,Tgrid::FACEDATA_UE); // Interpolation from nodes to faces
            // Upwind B. Approximate node_B in upstream of the node_UE. Displacement = mindx/2,
            // where mindx is the smallest dx of the cells touching the node. dt is not used.
            if(do_upwinding) {
                g.sph_CN_donor(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B, Tgrid::NODEDATA_UE,fp_dt);
            }
            g.sph_YNave(Tgrid::NODEDATA_B); // Theta Boundary conditions for nodes (averaging)
            // Interpolate J from cells to nodes, for the same reason
            g.sph_CN(Tgrid::CELLDATA_J,Tgrid::NODEDATA_J);
            //g.sph_YNave(Tgrid::NODEDATA_J); // Theta Boundary conditions for nodes (averaging)
            // Calculate E = -U_e x (B + B_0) + eta*J at the nodes
            g.sph_calc_node_E();
            g.sph_NC_copy(Tgrid::NODEDATA_E, Tgrid::CELLDATA_E);
            //g.sph_YNave(Tgrid::NODEDATA_E); // Theta Boundary conditions for nodes (averaging)
            //smooth the electric field before propagating B field.
            g.sph_smoothing_E();
            // Calculate -dB/fp_dt in cell faces by taking curl(E) from the nodes
            g.sph_FaceCurl_2(Tgrid::NODEDATA_E,Tgrid::FACEDATA_MINUSDB,1);
            // B_new = B_old - curl(E)*fp_dt
            g.FacePropagate(fsBold,fsBnew,fp_dt);
        } //end of BC_use_ghost_cell == 1
    } //end of propagation_type == 0
    else if(Params::sph_propagation_type == 1) { // Flat front propagation
        if(Params::sph_BC_use_ghost_cell == 0) {
            g.sph_FC(fsBrhs,Tgrid::CELLDATA_B);
            // Calculate B in cells (CELLDATA_UE)
            g.sph_Neumann_0(Tgrid::CELLDATA_B);
            g.sph_CNb(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B);
            g.sph_YNave(Tgrid::NODEDATA_B); // Theta Boundary conditions for nodes (averaging)
            //g.sph_NC_copy(Tgrid::NODEDATA_B, Tgrid::CELLDATA_B_TEMP2);
            if(Params::useJstag == 1) { // alternative j calculation
                // J TO NODES CALCULATION.
                g.sph_boundary_faces(Tgrid::FACEDATA_B);
                //g.sph_calc_node_j();
                g.sph_calc_node_j_2();
                g.sph_YNave(Tgrid::NODEDATA_J);
                g.sph_NC(Tgrid::NODEDATA_J,Tgrid::CELLDATA_J);
                g.sph_NC_copy(Tgrid::NODEDATA_J, Tgrid::CELLDATA_J_TEMP);
            } else {
                // Calculate J in cell faces by taking curl(B)/mu_0 from the nodes
                g.sph_FaceCurl_2(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
                // Interpolate J from cell faces to cells (interior cells only)
                g.sph_FC(Tgrid::FACEDATA_J,Tgrid::CELLDATA_J);
            }
            // Calculate U_e in cells (CELLDATA_UE)
            g.sph_calc_ue();
            g.sph_Neumann_0(Tgrid::CELLDATA_UE);
            if(Params::useNodeUe == 1) { // Calculate Ue separately to nodes
                g.sph_Neumann_rhoq_0();  // should already be done?
                g.sph_CNb_ne();
                g.sph_Neumann_0(Tgrid::CELLDATA_Ji);
                g.sph_CNb_Ji(Tgrid::CELLDATA_Ji,Tgrid::NODEDATA_Ji);
                g.sph_YNave_Ji(Tgrid::NODEDATA_Ji);
                g.sph_calc_node_ue();
                g.sph_YNave(Tgrid::NODEDATA_UE); // Theta Boundary conditions for nodes (averaging)
            } else { // no separate node Ue calculation
                // Interpolate U_e from cells to nodes
                g.sph_CNb(Tgrid::CELLDATA_UE,Tgrid::NODEDATA_UE);
                g.sph_YNave(Tgrid::NODEDATA_UE);
                // Ue - Spherical coordinates
            }
            // Interpolation from nodes to faces (To calculate first order of accuracy of Lorents forse)
            g.sph_NF(Tgrid::NODEDATA_UE,Tgrid::FACEDATA_UE);
            g.sph_Neumann_0(Tgrid::CELLDATA_B);
            if(do_upwinding) {
                g.sph_CNb_donor(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B, Tgrid::NODEDATA_UE,fp_dt);
            }
            g.sph_YNave(Tgrid::NODEDATA_B);
            //g.sph_NC_copy(Tgrid::NODEDATA_B, Tgrid::CELLDATA_B_TEMP2);
            g.sph_Node_BC(Tgrid::NODEDATA_B); // Magnetic boundary conditions
            g.sph_YNave(Tgrid::NODEDATA_B); // Theta Boundary conditions for nodes (averaging)
            // Interpolate J from cells to nodes, for the same reason (but dont overwrite if using jstag).
            if(Params::useJstag == 0) {
                g.sph_Neumann_0(Tgrid::CELLDATA_J);
                g.sph_CNb(Tgrid::CELLDATA_J,Tgrid::NODEDATA_J);
                g.sph_YNave(Tgrid::NODEDATA_J); // Theta Boundary conditions for nodes (averaging)
            }
            g.sph_calc_node_E();
            g.sph_YNave(Tgrid::NODEDATA_E); // Theta Boundary conditions for nodes (averaging)
            // E - Spherical coordiantes
            //smooth the electric field before propagating B field.
            g.sph_smoothing_E();
            if(Params::electricFieldSmoothingNumber > 0) {
                g.sph_YNave(Tgrid::NODEDATA_E);
            }
            g.sph_FaceCurl_2(Tgrid::NODEDATA_E,Tgrid::FACEDATA_MINUSDB,1);
            // B_new = B_old - curl(E)*fp_dt
            g.FacePropagate(fsBold,fsBnew,fp_dt);
        } //end of BC_use_ghost_cell == 0
        else if(Params::sph_BC_use_ghost_cell == 1) {
            g.sph_FC(fsBrhs,Tgrid::CELLDATA_B);
            //g.sph_Cell_BC(Tgrid::CELLDATA_B); // Magnetic boundary conditions
            // Calculate B in cells (CELLDATA_UE)
            g.sph_Neumann_SCS(Tgrid::CELLDATA_B);
            g.sph_CN(Tgrid::CELLDATA_B,Tgrid::NODEDATA_B);
            g.sph_YNave(Tgrid::NODEDATA_B);
            if(Params::useJstag == 1) { // alternative j calculation
                // J TO NODES CALCULATION.
                g.sph_boundary_faces(Tgrid::FACEDATA_B);
                //g.sph_calc_node_j();
                g.sph_calc_node_j_2();
                g.sph_YNave(Tgrid::NODEDATA_J);
                g.sph_NC(Tgrid::NODEDATA_J,Tgrid::CELLDATA_J);
            } else {
                // Calculate J in cell faces by taking curl(B)/mu_0 from the nodes
                g.sph_FaceCurl_2(Tgrid::NODEDATA_B,Tgrid::FACEDATA_J,1/Params::mu_0);
                // Interpolate J from cell faces to cells (interior cells only)
                g.sph_FC(Tgrid::FACEDATA_J,Tgrid::CELLDATA_J);
            }
            // Calculate U_e in cells (CELLDATA_UE)
            g.sph_calc_ue();
            g.sph_Neumann_SCS(Tgrid::CELLDATA_UE);
            if(Params::useNodeUe == 1) { // Calculate Ue separately to nodes
                g.sph_Neumann_rhoq();  // should already be done?
                //g.sph_CN_ne();
                g.CN_ne();
                g.sph_Neumann(Tgrid::CELLDATA_Ji);
                g.sph_CN_Ji(Tgrid::CELLDATA_Ji,Tgrid::NODEDATA_Ji);
                g.sph_YNave_Ji(Tgrid::NODEDATA_Ji);
                g.sph_calc_node_ue();
                g.sph_YNave(Tgrid::NODEDATA_UE); // Theta Boundary conditions for nodes (averaging)
            } else { // no separate node Ue calculation
                // Interpolate U_e from cells to nodes
                g.sph_CN(Tgrid::CELLDATA_UE,Tgrid::NODEDATA_UE);
                g.sph_YNave(Tgrid::NODEDATA_UE);
                // Ue - Spherical coordinates
            }
            // Interpolation from nodes to faces (To calculate first order of accuracy of Lorents forse)
            g.sph_NF(Tgrid::NODEDATA_UE,Tgrid::FACEDATA_UE);
            g.sph_Neumann_SCS(Tgrid::CELLDATA_B);
            if(do_upwinding) {
                g.sph_CN_donor(Tgrid::CELLDATA_B, Tgrid::NODEDATA_B, Tgrid::NODEDATA_UE,fp_dt);
            }
            g.sph_YNave(Tgrid::NODEDATA_B);
            g.sph_Node_BC(Tgrid::NODEDATA_B); // Magnetic boundary conditions
            g.sph_YNave(Tgrid::NODEDATA_B);
            // Interpolate J from cells to nodes, for the same reason (but dont overwrite if using jstag).
            if(Params::useJstag == 0) {
                g.sph_Neumann_SCS(Tgrid::CELLDATA_J);
                g.sph_CN(Tgrid::CELLDATA_J,Tgrid::NODEDATA_J);
                g.sph_YNave(Tgrid::NODEDATA_J);
            }
            g.sph_calc_node_E();
            g.sph_YNave(Tgrid::NODEDATA_E);
            g.sph_NC_copy(Tgrid::NODEDATA_E, Tgrid::CELLDATA_E);
            // E - Spherical coordiantes
            //smooth the electric field before propagating B field.
            g.sph_smoothing_E();
            if (Params::electricFieldSmoothingNumber > 0) {
                g.sph_YNave(Tgrid::NODEDATA_E);
            }
            g.sph_FaceCurl_2(Tgrid::NODEDATA_E,Tgrid::FACEDATA_MINUSDB,1);
            g.sph_FC(Tgrid::FACEDATA_MINUSDB,Tgrid::CELLDATA_B_TEMP1); // Only for diagnostic
            // B_new = B_old - curl(E)*fp_dt
            g.FacePropagate(fsBold,fsBnew,fp_dt);
        } //end of BC_use_ghost_cell == 1
    } //end of propagation_type == 1
}

#endif

