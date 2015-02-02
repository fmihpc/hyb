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

#ifndef BACKGROUNDCHARGE_H
#define BACKGROUNDCHARGE_H

#include <cmath>
#include <cstdio>
#include "definitions.h"
#include "atmosphere.h"

//! Background charge density profile
class BackgroundChargeDensityProfile : public ScalarField
{
public:
    BackgroundChargeDensityProfile();
    BackgroundChargeDensityProfile(std::string funcName,std::vector<real> args);
    ~BackgroundChargeDensityProfile();
    real getValue(const gridreal r[]);
private:
    real (BackgroundChargeDensityProfile::*ptr)(const gridreal[]);
    real defaultFunction(const gridreal[]);
    // BACKGROUND CHARGE DENSITY PROFILES
    real smoothObstacle(const gridreal[]);
    real sharpObstacle(const gridreal[]);
    real bgCosSzaDayConstantNight(const gridreal[]);
    real bgTitan(const gridreal[]);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    real sph_smoothObstacle(const gridreal[]);
#endif
};

void initializeBackgroundChargeDensity();

#endif

