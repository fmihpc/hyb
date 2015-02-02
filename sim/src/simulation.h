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

#ifndef SIMULATION_H
#define SIMULATION_H

#include "grid.h"
#include "timepool.h"
#include "params.h"
#include "logger.h"
#include "vis/vis_db.h"
#include "vis/vis_data_source_simulation.h"
#include "diagnostics.h"

extern Logger mainlog, errorlog, paramslog;

//! Simulation main class
class Simulation
{
public:
    Simulation();
    ~Simulation();
    void run();
    void runOnlyParameterDynamics();
    void updateParams();
    int finalize();
    static bool PropagateV(TLinkedParticle& part);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    static bool sph_PropagateV(TLinkedParticle& part);
#endif
private:
    Ttimepool timepool;
    GridRefinementProfile refineFunc;
    real macroParticlePropagations;
    SimulationVisDataSourceImpl* visDataSourceImpl;
    std::vector<VisDB*> visWriters;
    void initializeSimulation();
    void stepForward();
    bool finalizeTimestep(bool doBreakpointing = true);
    void saveStep();
    void dumpState(const char *fileName);
    void readState(const char *fileName);
    void saveVisualizationFiles();
    void saveExtraHcFiles();
    static bool output(TLinkedParticle&);
    static bool AlwaysTrue(TLinkedParticle&);
    static void BoundaryB(datareal celldata[Tgrid::NCELLDATA][3], int dim);
    static bool PropagateX(TLinkedParticle& part);
    void fieldpropagate(Tgrid::TFaceDataSelect fsBnew, Tgrid::TFaceDataSelect fsBold,
                        Tgrid::TFaceDataSelect fsBrhs,
                        real fp_dt, bool do_upwinding);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    void sph_stepForward();
    static void sph_BoundaryB(datareal celldata[Tgrid::NCELLDATA][3], gridreal sph_centroid[3], int dim);
    static bool sph_PropagateX(TLinkedParticle& part);
    void sph_fieldpropagate(Tgrid::TFaceDataSelect fsBnew, Tgrid::TFaceDataSelect fsBold,Tgrid::TFaceDataSelect fsBrhs,real fp_dt, bool do_upwinding);
#endif
#ifdef USE_PARTICLE_SUBCYCLING
    static bool PropagatePart1(TLinkedParticle&, ParticlePassArgs);
    static bool PropagatePart2(TLinkedParticle&);
#endif
};

#endif

