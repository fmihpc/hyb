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

#ifndef POPULATION_EXOSPHERIC_H
#define POPULATION_EXOSPHERIC_H

//! Exospheric population
class PopulationExospheric : public Population
{
public:
    PopulationExospheric(PopulationArgs args);
    ~PopulationExospheric();
    void initialize();
    void createParticles();
    real getNeutralDensity(const gridreal[]);
    void updateArgs();
    void writeExtraHcFile();
    std::string configDump();
    std::string toString();
private:
    real R;
    real totalRate;
    MultipleProductDistribution distFunc;
    TPDF_ID distFuncId;
    void newParticle();
    void writeLog();
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    void sph_newParticle();
#endif
};

#endif

