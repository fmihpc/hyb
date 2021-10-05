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

#include "diagnostics.h"
#include "simulation.h"

using namespace std;

//! Constructor
Diagnostics::Diagnostics()
{
    static bool onlyOneInstance = false;
    if(onlyOneInstance == true) {
        ERRORMSG("only one Diagnostics instance can be created");
        doabort();
    } else {
        onlyOneInstance = true;
    }
}

//! Destructor
Diagnostics::~Diagnostics()
{
    // Delete dynamically allocated memory
    for(unsigned int i=0; i < pCounter.size(); ++i) {
        delete pCounter[i];
    }
    pCounter.clear();
    // Close files
    for(unsigned int i=0; i < plog.size(); ++i) {
        plog[i]->flush();
        plog[i]->close();
        delete plog[i];
    }
    plog.clear();
    flog.flush();
    flog.close();
}

//! Initialize particle counters and create log files
void Diagnostics::init()
{
    for(unsigned int i=0; i < Params::pops.size(); ++i) {
        // particle counter
        pCounter.push_back(new ParticleCounter(Params::pops[i]->getPopId()));
        // log file
        stringstream fn;
        fn << "pop" << int2string(i+1,3) << "_" << Params::pops[i]->getIdStr() << ".log";
        plog.push_back(new ofstream(fn.str().c_str()));
        (*plog[i]) << scientific << showpos;
        (*plog[i]).precision(10);
    }
}

//! Run diagnostics
void Diagnostics::run()
{
    logParticles();
    logFields();
}

//! Calculate particle parameters (used when passing through particle list)
bool Diagnostics::particleAnalyzeFunction(TLinkedParticle& part)
{
#ifndef NO_DIAGNOSTICS
    // Number of real physical particles [#]
    Params::diag.pCounter[part.popid]->totalWeight += part.w;
    // Number of macroparticles [#]
    Params::diag.pCounter[part.popid]->macroParticles += 1;
    // Average velocity component and magnitude in particle populations [m/s]
    Params::diag.pCounter[part.popid]->avgVx += part.vx*part.w;
    Params::diag.pCounter[part.popid]->avgVy += part.vy*part.w;
    Params::diag.pCounter[part.popid]->avgVz += part.vz*part.w;
    Params::diag.pCounter[part.popid]->avgV += sqrt( sqr(part.vx) + sqr(part.vy) + sqr(part.vz))*part.w;
    // Kinetic energy in real physical particles [J]
    Params::diag.pCounter[part.popid]->kineticEnergy += 0.5*part.w*Params::pops[part.popid]->m*( sqr(part.vx) + sqr(part.vy) + sqr(part.vz));
#endif

    return true;
}

//! Write particle population log files
void Diagnostics::logParticles()
{
    static bool initDone = false;
    if(initDone == false) {
        // Create particle log file headers
        for(unsigned int i=0; i < pCounter.size(); ++i) {
            (*plog[i])
                    << "% " << Params::pops[i]->getIdStr() << "\n"
                    << "% m [kg] = " << Params::pops[i]->m << "\n"
                    << "% q [C]  = " << Params::pops[i]->q << "\n"
                    << "% columns = 43\n"
                    << "% 01. Time [s]\n"
                    << "% 02. Particles [#]\n"
                    << "% 03. Macroparticles [#]\n"
                    << "% 04. Splitting rate [#/dt]\n"
                    << "% 05. Joining rate [#/dt]\n"
                    << "% 06. Charge exchange rate [#/dt]\n"
                    << "% 07. Electron impact ionization rate [#/dt]\n"
                    << "% 08. avg(vx) [m/s]\n"
                    << "% 09. avg(vy) [m/s]\n"
                    << "% 10. avg(vz) [m/s]\n"
                    << "% 11. avg(|v|) [m/s]\n"
                    << "% 12. cutV rate [#/dt]\n"
                    << "% 13. Kinetic energy [J]\n"
                    << "% 14. Escape rate (total) [#/s]\n"
                    << "% 15. Escape rate (front wall) [#/s]\n"
                    << "% 16. Escape rate (back wall) [#/s]\n"
                    << "% 17. Escape rate (side walls) [#/s]\n"
                    << "% 18. Escape x-momentum rate (total) [kgm/s^2]\n"
                    << "% 19. Escape y-momentum rate (total) [kgm/s^2]\n"
                    << "% 20. Escape z-momentum rate (total) [kgm/s^2]\n"
                    << "% 21. Escape x-momentum rate (front wall) [kgm/s^2]\n"
                    << "% 22. Escape y-momentum rate (front wall) [kgm/s^2]\n"
                    << "% 23. Escape z-momentum rate (front wall) [kgm/s^2]\n"
                    << "% 24. Escape x-momentum rate (back wall) [kgm/s^2]\n"
                    << "% 25. Escape y-momentum rate (back wall) [kgm/s^2]\n"
                    << "% 26. Escape z-momentum rate (back wall) [kgm/s^2]\n"
                    << "% 27. Escape x-momentum rate (side walls) [kgm/s^2]\n"
                    << "% 28. Escape y-momentum rate (side walls) [kgm/s^2]\n"
                    << "% 29. Escape z-momentum rate (side walls) [kgm/s^2]\n"
                    << "% 30. Escape kinetic energy rate (total) [J/s]\n"
                    << "% 31. Escape kinetic energy rate (front wall) [J/s]\n"
                    << "% 32. Escape kinetic energy rate (back wall) [J/s]\n"
                    << "% 33. Escape kinetic energy rate (side walls) [J/s]\n"
                    << "% 34. Impact rate [#/s]\n"
                    << "% 35. Impact x-momentum rate [kgm/s^2]\n"
                    << "% 36. Impact y-momentum rate [kgm/s^2]\n"
                    << "% 37. Impact z-momentum rate [kgm/s^2]\n"
                    << "% 38. Impact kinetic energy rate [J/s]\n"
                    << "% 39. Inject rate [#/s]\n"
                    << "% 40. Inject x-momentum rate [kgm/s^2]\n"
                    << "% 41. Inject y-momentum rate [kgm/s^2]\n"
                    << "% 42. Inject z-momentum rate [kgm/s^2]\n"
                    << "% 43. Inject kinetic energy rate [J/s]\n"
                    << flush;
        }
        initDone = true;
    }
    // Go through all particles and do analysis stuff
    g.particle_pass(&particleAnalyzeFunction);
    // Finalize counters
    for (unsigned int i = 0; i < pCounter.size(); i++) {
        pCounter[i]->finalizeCounters();
    }
    // Go thru all populations
    for (unsigned int i = 0; i < pCounter.size(); i++) {
        (*plog[i])
                << Params::t << "  " << "\t"
                << pCounter[i]->totalWeight << "\t"
                << pCounter[i]->macroParticles << "\t"
                << pCounter[i]->splittingRate << "\t"
                << pCounter[i]->joiningRate << "\t"
                << pCounter[i]->chargeExchangeRate << "\t"
                << pCounter[i]->electronImpactIonizationRate << "\t"
                << pCounter[i]->avgVx << "\t"
                << pCounter[i]->avgVy << "\t"
                << pCounter[i]->avgVz << "\t"
                << pCounter[i]->avgV << "\t"
                << pCounter[i]->cutRateV << "\t"
                << pCounter[i]->kineticEnergy << "\t"
                << pCounter[i]->escapeRateParticles << "\t"
                << pCounter[i]->escapeRateParticlesFrontWall << "\t"
                << pCounter[i]->escapeRateParticlesBackWall << "\t"
                << pCounter[i]->escapeRateParticlesSideWall << "\t"
                << pCounter[i]->escapeRateMomentum[0] << "\t"
                << pCounter[i]->escapeRateMomentum[1] << "\t"
                << pCounter[i]->escapeRateMomentum[2] << "\t"
                << pCounter[i]->escapeRateMomentumFrontWall[0] << "\t"
                << pCounter[i]->escapeRateMomentumFrontWall[1] << "\t"
                << pCounter[i]->escapeRateMomentumFrontWall[2] << "\t"
                << pCounter[i]->escapeRateMomentumBackWall[0] << "\t"
                << pCounter[i]->escapeRateMomentumBackWall[1] << "\t"
                << pCounter[i]->escapeRateMomentumBackWall[2] << "\t"
                << pCounter[i]->escapeRateMomentumSideWall[0] << "\t"
                << pCounter[i]->escapeRateMomentumSideWall[1] << "\t"
                << pCounter[i]->escapeRateMomentumSideWall[2] << "\t"
                << pCounter[i]->escapeRateKineticEnergy << "\t"
                << pCounter[i]->escapeRateKineticEnergyFrontWall << "\t"
                << pCounter[i]->escapeRateKineticEnergyBackWall << "\t"
                << pCounter[i]->escapeRateKineticEnergySideWall << "\t"
                << pCounter[i]->impactRateParticles << "\t"
                << pCounter[i]->impactRateMomentum[0] << "\t"
                << pCounter[i]->impactRateMomentum[1] << "\t"
                << pCounter[i]->impactRateMomentum[2] << "\t"
                << pCounter[i]->impactRateKineticEnergy << "\t"
                << pCounter[i]->injectRateParticles << "\t"
                << pCounter[i]->injectRateMomentum[0] << "\t"
                << pCounter[i]->injectRateMomentum[1] << "\t"
                << pCounter[i]->injectRateMomentum[2] << "\t"
                << pCounter[i]->injectRateKineticEnergy << "\t"
                << "\n" << flush;
    }
    // Reset counters
    for(unsigned int i = 0; i < pCounter.size(); ++i) {
        pCounter[i]->reset();
    }
}

//! Write field log file
void Diagnostics::logFields()
{
    static bool initDone = false;
    if(initDone == false) {
        flog.open("field.log");
        flog << scientific << showpos;
        flog.precision(10);
        flog
                << "% field\n"
                << "% columns = 14\n"
                << "% 01. Time [s]\n"
                << "% 02. avg(Bx) [T]\n"
                << "% 03. avg(By) [T]\n"
                << "% 04. avg(Bz) [T]\n"
                << "% 05. avg(B) [T]\n"
                << "% 06. max(B) [T]\n"
                << "% 07. avg(div(B)) [T/m]\n"
                << "% 08. max(div(B)) [T/m]\n"
                << "% 09. max(dx*div(B)/B) [-]\n"
                << "% 10. energy(sum(dV*B^2/2*mu0)) [J]\n"
                << "% 11. cutE rate [#/dt]\n"
                << "% 12. cutRhoQ rate [#/dt]\n"
                << "% 13. cutUe rate [#/dt]\n"
                << "% 14. maxVw rate [#/dt]\n"
                << flush;
        initDone = true;
    }
    // Finalize results
    Tgrid::fieldCounter.finalize();
    // Calculate instantenous field values
    MagneticLog magLog;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    g.calc_facediv(Tgrid::FACEDATA_B,magLog);
#else
    g.sph_calc_facediv(Tgrid::FACEDATA_B,Tgrid::CELLDATA_B,magLog);
#endif
    flog << Params::t << "\t";
    flog << magLog.avgBx << "\t";
    flog << magLog.avgBy << "\t";
    flog << magLog.avgBz << "\t";
    flog << magLog.avgB << "\t";
    flog << magLog.maxB << "\t";
    flog << magLog.avgDivB << "\t";
    flog << magLog.maxDivB << "\t";
    flog << magLog.maxDxDivBperB << "\t";
    flog << magLog.energyB << "\t";
    flog << Tgrid::fieldCounter.cutRateE << "\t";
    flog << Tgrid::fieldCounter.cutRateRhoQ << "\t";
    flog << Tgrid::fieldCounter.cutRateUe << "\t";
    flog << Tgrid::fieldCounter.cutRateMaxVw << "\t";
    flog << "\n" << flush;
    // Maximum B field reached => set program termination flag and save hc-files
    if (magLog.maxB >= Params::B_limit) {
        errorlog << "STOPPING: maximum B field reached (" << Params::B_limit*1e9 << " nT) at " << Tr3v(magLog.posMaxB).toString() << "\n";
        Params::stoppingPhase = true;
    }
    // Reset field counters
    Tgrid::fieldCounter.reset();
}

//! Constructor
ParticleCounter::ParticleCounter(const int populationid) : popid(populationid)
{
    reset();
}

//! Reset particle counters
void ParticleCounter::reset()
{
    totalWeight = 0.0;
    macroParticles = 0.0;
    avgVx = 0.0;
    avgVy = 0.0;
    avgVz = 0.0;
    avgV = 0.0;
    cutRateV = 0.0;
    kineticEnergy = 0.0;
    splittingRate = 0.0;
    joiningRate = 0.0;
    chargeExchangeRate = 0.0;
    escapeRateParticles = 0.0;
    escapeRateParticlesFrontWall = 0.0;
    escapeRateParticlesBackWall = 0.0;
    escapeRateParticlesSideWall = 0.0;
    for(int i=0; i<3; i++) {
        escapeRateMomentum[i] = 0.0;
        escapeRateMomentumFrontWall[i] = 0.0;
        escapeRateMomentumBackWall[i] = 0.0;
        escapeRateMomentumSideWall[i] = 0.0;
        impactRateMomentum[i] = 0.0;
        injectRateMomentum[i] = 0.0;
    }
    escapeRateKineticEnergy = 0.0;
    escapeRateKineticEnergyFrontWall = 0.0;
    escapeRateKineticEnergyBackWall = 0.0;
    escapeRateKineticEnergySideWall = 0.0;
    impactRateParticles = 0.0;
    impactRateKineticEnergy = 0.0;
    injectRateParticles = 0.0;
    injectRateKineticEnergy = 0.0;
    electronImpactIonizationRate = 0.0;
    resetTime = Params::t;
    resetTimestep = Params::cnt_dt;
}

//! Increase escape counter (front wall)
void ParticleCounter::increaseEscapeCountersFrontWall(TLinkedParticle& p)
{
    escapeRateParticlesFrontWall += p.w;
    escapeRateMomentumFrontWall[0] += p.w*p.vx;
    escapeRateMomentumFrontWall[1] += p.w*p.vy;
    escapeRateMomentumFrontWall[2] += p.w*p.vz;
    escapeRateKineticEnergyFrontWall += p.w*(sqr(p.vx) + sqr(p.vy) + sqr(p.vz));
}

//! Increase escape counter (back wall)
void ParticleCounter::increaseEscapeCountersBackWall(TLinkedParticle& p)
{
    escapeRateParticlesBackWall += p.w;
    escapeRateMomentumBackWall[0] += p.w*p.vx;
    escapeRateMomentumBackWall[1] += p.w*p.vy;
    escapeRateMomentumBackWall[2] += p.w*p.vz;
    escapeRateKineticEnergyBackWall += p.w*(sqr(p.vx) + sqr(p.vy) + sqr(p.vz));
}

//! Increase escape counter (side walls)
void ParticleCounter::increaseEscapeCountersSideWall(TLinkedParticle& p)
{
    escapeRateParticlesSideWall += p.w;
    escapeRateMomentumSideWall[0] += p.w*p.vx;
    escapeRateMomentumSideWall[1] += p.w*p.vy;
    escapeRateMomentumSideWall[2] += p.w*p.vz;
    escapeRateKineticEnergySideWall += p.w*(sqr(p.vx) + sqr(p.vy) + sqr(p.vz));
}

//! Increase impact counters (inner boundary)
void ParticleCounter::increaseImpactCounters(TLinkedParticle& p)
{
    impactRateParticles += p.w;
    impactRateMomentum[0] += p.w*p.vx;
    impactRateMomentum[1] += p.w*p.vy;
    impactRateMomentum[2] += p.w*p.vz;
    impactRateKineticEnergy += p.w*(sqr(p.vx) + sqr(p.vy) + sqr(p.vz));
}

//! Increase inject counters
void ParticleCounter::increaseInjectCounters(shortreal vx,shortreal vy,shortreal vz,shortreal w)
{
    injectRateParticles += w;
    injectRateMomentum[0] += w*vx;
    injectRateMomentum[1] += w*vy;
    injectRateMomentum[2] += w*vz;
    injectRateKineticEnergy += w*(sqr(vx) + sqr(vy) + sqr(vz));
}

//! Finalize counters
void ParticleCounter::finalizeCounters()
{
    // counting time & number of timesteps
    real Dt = Params::t - resetTime;
    real timeSteps = static_cast<real>(1+Params::cnt_dt - resetTimestep);
    // temporal averages
    if(Dt > 0) {
        escapeRateParticlesFrontWall /= Dt;
        escapeRateParticlesBackWall /= Dt;
        escapeRateParticlesSideWall /= Dt;
        for(int i=0; i<3; i++) {
            escapeRateMomentumFrontWall[i] *= Params::pops[popid]->m/Dt;
            escapeRateMomentumBackWall[i] *= Params::pops[popid]->m/Dt;
            escapeRateMomentumSideWall[i] *= Params::pops[popid]->m/Dt;
            impactRateMomentum[i] *= Params::pops[popid]->m/Dt;
            injectRateMomentum[i] *= Params::pops[popid]->m/Dt;
        }
        escapeRateKineticEnergyFrontWall *= 0.5*Params::pops[popid]->m/Dt;
        escapeRateKineticEnergyBackWall *= 0.5*Params::pops[popid]->m/Dt;
        escapeRateKineticEnergySideWall *= 0.5*Params::pops[popid]->m/Dt;
        impactRateParticles /= Dt;
        impactRateKineticEnergy *= 0.5*Params::pops[popid]->m/Dt;
        injectRateParticles /= Dt;
        injectRateKineticEnergy *= 0.5*Params::pops[popid]->m/Dt;
    } else {
        escapeRateParticlesFrontWall = 0.0;
        escapeRateParticlesBackWall = 0.0;
        escapeRateParticlesSideWall = 0.0;
        for(int i=0; i<3; i++) {
            escapeRateMomentumFrontWall[i] = 0.0;
            escapeRateMomentumBackWall[i] = 0.0;
            escapeRateMomentumSideWall[i] = 0.0;
            impactRateMomentum[i] = 0.0;
            injectRateMomentum[i] = 0.0;
        }
        escapeRateKineticEnergyFrontWall = 0.0;
        escapeRateKineticEnergyBackWall = 0.0;
        escapeRateKineticEnergySideWall = 0.0;
        impactRateParticles = 0.0;
        impactRateKineticEnergy = 0.0;
        injectRateParticles = 0.0;
        injectRateKineticEnergy = 0.0;
    }
    // sum total escape rates
    escapeRateParticles = escapeRateParticlesFrontWall +
                          escapeRateParticlesBackWall +
                          escapeRateParticlesSideWall;
    for(int i=0; i<3; i++) {
        escapeRateMomentum[i] = escapeRateMomentumFrontWall[i] +
                                escapeRateMomentumBackWall[i] +
                                escapeRateMomentumSideWall[i];
    }
    escapeRateKineticEnergy = escapeRateKineticEnergyFrontWall +
                              escapeRateKineticEnergyBackWall +
                              escapeRateKineticEnergySideWall;
    // timestep averages
    if(timeSteps > 0) {
        splittingRate /= timeSteps;
        joiningRate /= timeSteps;
        chargeExchangeRate /= timeSteps;
        electronImpactIonizationRate /= timeSteps;
        cutRateV /= timeSteps;
    } else {
        splittingRate = 0;
        joiningRate = 0;
        chargeExchangeRate = 0;
        electronImpactIonizationRate = 0;
        cutRateV = 0;
    }
    // averages over a whole population
    if (totalWeight > 0) {
        avgVx /= totalWeight;
        avgVy /= totalWeight;
        avgVz /= totalWeight;
        avgV /= totalWeight;
    } else {
        avgVx = 0;
        avgVy = 0;
        avgVz = 0;
        avgV = 0;
    }
}

