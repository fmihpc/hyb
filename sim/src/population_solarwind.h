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

#ifndef POPULATION_SOLARWIND_H
#define POPULATION_SOLARWIND_H

//! Solar wind population
class PopulationSolarWind : public Population
{
public:
    PopulationSolarWind(PopulationArgs args);
    ~PopulationSolarWind();
    void initialize();
    void createParticles();
    void addParticle(shortreal x,shortreal y,shortreal z,real w);
    void updateArgs();
    void writeExtraHcFile();
    std::string configDump();
    std::string toString();
private:
    real n;
    real V;
    real backWallWeight;
    bool negativeV;
    void newParticle();
    void writeLog();
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    void sph_newParticle();
#endif
};

#endif

