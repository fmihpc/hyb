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

#ifndef DIAGNOSTICS_H
#define DIAGNOSTICS_H

#include <fstream>
#include <vector>
#include "definitions.h"
#include "particle.h"

struct ParticleCounter;

//! Simulation diagnostics
class Diagnostics
{
public:
    Diagnostics();
    ~Diagnostics();
    void init();
    void run();
    static bool particleAnalyzeFunction(TLinkedParticle& p);
    std::vector<ParticleCounter*> pCounter; //!< Particle counters
private:
    std::vector<std::ofstream*> plog; //!< Particle population log files
    std::ofstream flog; //!< Field log file
    void logParticles();
    void logFields();
};

//! Particle counters
struct ParticleCounter {
    const int popid;
    // Snapshot counters
    real totalWeight;
    real macroParticles;
    real avgVx;
    real avgVy;
    real avgVz;
    real avgV;
    real cutRateV;
    real kineticEnergy;
    // Integrating counters
    real splittingRate;
    real joiningRate;
    real chargeExchangeRate;
    real escapeRateParticles;
    real escapeRateParticlesFrontWall;
    real escapeRateParticlesBackWall;
    real escapeRateParticlesSideWall;
    real escapeRateMomentum[3];
    real escapeRateMomentumFrontWall[3];
    real escapeRateMomentumBackWall[3];
    real escapeRateMomentumSideWall[3];
    real escapeRateKineticEnergy;
    real escapeRateKineticEnergyFrontWall;
    real escapeRateKineticEnergyBackWall;
    real escapeRateKineticEnergySideWall;
    real impactRateParticles;
    real impactRateMomentum[3];
    real impactRateKineticEnergy;
    real injectRateParticles;
    real injectRateMomentum[3];
    real injectRateKineticEnergy;
    real electronImpactIonizationRate;
    real resetTime;
    int resetTimestep;
    ParticleCounter(const int populationid);
    void reset();
    void increaseEscapeCountersFrontWall(TLinkedParticle& p);
    void increaseEscapeCountersBackWall(TLinkedParticle& p);
    void increaseEscapeCountersSideWall(TLinkedParticle& p);
    void increaseImpactCounters(TLinkedParticle& p);
    void increaseInjectCounters(shortreal vx,shortreal vy,shortreal vz,shortreal w);
    void finalizeCounters();
};

#endif

