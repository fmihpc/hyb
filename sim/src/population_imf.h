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

#ifndef POPULATION_IMF_H
#define POPULATION_IMF_H

#include "definitions.h"

//! IMF-directed population
class PopulationIMF : public Population
{
public:
    PopulationIMF(PopulationArgs args);
    ~PopulationIMF();
    void initialize();
    void createParticles();
    void updateArgs();
    void writeExtraHcFile();
    bool updateDistributions();
    std::string configDump();
    std::string toString();
private:
    real n;
    real V;
    double V_avg;
    real A0, usw, Bangle, uswCrossBdirMagn, k, v1;
    Tr3v rotVec, vVec;
    real rotAngle;
    real backWallWeight;
    bool negativeV, conserveE;
    void newParticle();
    void writeLog();
    real t0;
    int E_distr_type;
    real xmin, xmax, ymin, ymax, zmin, zmax, box_x, box_y, box_z;
    int pitch_distr_type;
    real pitch_cutoff, pitch_width;
    real CartesianIntegrator(Tr3v uswVec,real parkerAngle,real clockAngle);
    real EnergyPitchIntegrator(Tr3v uswVec,real parkerAngle,real clockAngle), vf(real v, real pitch);
    real pitchf(real pitch);
    int dtlevel;
};

#endif

