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

#ifndef ATMOSPHERE_H
#define ATMOSPHERE_H

#include "definitions.h"
#include <string>
#include <vector>

//! Real valued scalar field
class ScalarField
{
public:
    ScalarField();
    virtual ~ScalarField();
    virtual real getValue(const gridreal[]); 
    bool isDefined();
    std::string toString();
protected:
    std::string name; //!< Name of the the scalar field
    std::vector<real> args; //!< Arguments of the scalar field
private:
};

//! Spatial distribution function of atmospheric particle populations
class SpatialDistribution : public ScalarField
{
public:
    SpatialDistribution();
    SpatialDistribution(std::string funcName,std::vector<real> args,unsigned int popid);
    ~SpatialDistribution();
    real getValue(const gridreal r[]);
private:
    unsigned int popid;
    real (SpatialDistribution::*ptr)(const gridreal[]);
    real defaultFunction(const gridreal[]);
    // SPATIAL DISTRIBUTION PROFILES
    real ionoConstantDayConstantNight(const gridreal[]);
    real ionoCosSzaDayConstantNight(const gridreal[]);
    real ionoCosSzaNoonToMidnight(const gridreal[]);
    real ionoLinearSzaDayConstantNight(const gridreal[]);
    real ionoLinearSzaNoonToMidnight(const gridreal[]);
    real ionoCosSzaConstNight_anySolarDirection(const gridreal[]);
    real ionoConstDayAndNight_anySolarDirection(const gridreal[]);
    real ionoCosSza_anySolarDirection(const gridreal[]);
    real ionoSmoothDistribution(const gridreal[]);
    real neutralDensityVenusHydrogen(const gridreal[]);
    real neutralDensityVenusOxygenHot(const gridreal[]);
    real photoionChamberlainTitan(const gridreal[]);
    real neutralDensityChamberlainT(const gridreal[]);
    real neutralDensityChamberlainH(const gridreal[]);
    real neutralDensityPowerLaw(const gridreal[]);
    real neutralDensityExponential(const gridreal[]);
    bool shadowTitan(const gridreal[], const real);
    real ionizationConstant(const gridreal[]);
    real shadow(const gridreal[]);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    real sph_ionoCosSzaDayConstantNight(const gridreal[]);
    real sph_neutralDensityVenusHydrogen(const gridreal[]);
    real sph_neutralDensityVenusOxygenHot(const gridreal[]);
    real sph_shadow(const gridreal[]);
#endif
};

//! Product of multiple spatial distribution functions
class MultipleProductDistribution : public ScalarField
{
public:
    MultipleProductDistribution();
    MultipleProductDistribution(std::vector<SpatialDistribution> distFuncs);
    ~MultipleProductDistribution();
    real getValue(const gridreal r[]);
    real getValue(const gridreal r[],const unsigned int n);
    std::string toString(std::string prefix=std::string(""),std::string delim=std::string("\n"));
private:
    real (MultipleProductDistribution::*ptr)(const gridreal[]);
    std::vector<SpatialDistribution> distFuncs;
    real defaultFunction(const gridreal[]);
    real productFunction(const gridreal[]);
};

real CosSZA(const gridreal[]);

#endif

