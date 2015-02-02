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

#ifndef RESISTIVITY_H
#define RESISTIVITY_H

#include <cmath>
#include <cstdio>
#include "definitions.h"
#include "atmosphere.h"

//! Resistivity profiles
class ResistivityProfile : public ScalarField
{
public:
    ResistivityProfile();
    ResistivityProfile(std::string funcName,std::vector<real> args);
    ~ResistivityProfile();
    real getValue(const gridreal r[]);
private:
    real (ResistivityProfile::*ptr)(const gridreal[]);
    real defaultFunction(const gridreal r[]);
    // RESISTIVITY PROFILES
    real resistivityConstant(const gridreal r[]);
    void setArgs_resistivityConstant();
    real resistivityConstantOutsideR(const gridreal r[]);
    void setArgs_resistivityConstantOutsideR();
    real resistivityWithinObstacle(const gridreal r[]);
    void setArgs_resistivityWithinObstacle();
    real resistivitySpherical(const gridreal r[]);
    void setArgs_resistivitySpherical();
    real resistivityCartesian(const gridreal r[]);
    void setArgs_resistivityCartesian();
    real resistivitySphericalPolynomial(const gridreal r[]);
    void setArgs_resistivitySphericalPolynomial();
    // RESISTIVITY PARAMETERS
    void resetParameters();
    real eta,eta0_c,R,R2,dR_coef,dR,res1,res2,x1,y1,z1,x2,y2,z2;
    std::vector<real> etaArray;
    std::vector<real> rMinSqrArray;
    std::vector<real> rMaxSqrArray;
    std::vector<real> polyCoeffArray;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    real sph_resistivityConstantOutsideR(const gridreal r[]);
    void setArgs_sph_resistivityConstantOutsideR();
#endif
};

void setResistivity();
void initializeResistivity();

#endif

