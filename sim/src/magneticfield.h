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

#ifndef MAGNETICFIELD_H
#define MAGNETICFIELD_H

#include "definitions.h"
#include <string>
#include <vector>

//! Real valued vector field
class VectorField
{
public:
    VectorField();
    virtual ~VectorField();
    virtual void getValue(const gridreal r[],datareal V[]);
    bool isDefined();
protected:
    std::string name; //!< Name of the the vector field
    std::vector<real> args; //!< Arguments of the vector field
};

//! Magnetic field profiles
class MagneticFieldProfile : public VectorField
{
public:
    MagneticFieldProfile();
    MagneticFieldProfile(std::string funcName,std::vector<real> args);
    ~MagneticFieldProfile();
    void getValue(const gridreal r[], datareal B[]);
private:
    void (MagneticFieldProfile::*ptr)(const gridreal r[],datareal B[]);
    void defaultFunction(const gridreal r[],datareal B[]);
    // MAGNETIC FIELD PROFILES
    void constantB(const gridreal r[],datareal B[]);
    void setArgs_constantB();
    void constantBx(const gridreal r[],datareal B[]);
    void setArgs_constantBx();
    void laminarFlowAroundSphereB(const gridreal r[],datareal B[]);
    void setArgs_laminarFlowAroundSphereB();
    void laminarFlowAroundSphereBx(const gridreal r[],datareal B[]);
    void setArgs_laminarFlowAroundSphereBx();
    void dipoleB(const gridreal r[],datareal B[]);
    void setArgs_dipoleB();
    void translateDipoleB(const gridreal r[],datareal B[]);
    void setArgs_translateDipoleB();
    void generalDipoleB(const gridreal r[],datareal B[]);
    void setArgs_generalDipoleB();
    void hemisphericDipoleB(const gridreal r[],datareal B[]);
    void setArgs_hemisphericDipoleB();
    void dipoleCuspB(const gridreal r[],datareal B[]);
    void setArgs_dipoleCuspB();
    // MAGNETIC FIELD PARAMETERS
    void resetParameters();
    real Bx,By,Bz,Btot,R,R3;
    real sinTheta,cosTheta,sinPhi,cosPhi,xOrigin,yOrigin,zOrigin;
    real dipSurfB,dipSurfR,dipMomCoeff,dipRmin2,hemiCoeffDip,hemiCoeffQuad;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    int  sph_dipole_dir;
    void sph_laminarFlowAroundSphereBx(const gridreal r[],datareal B[]);
    void setArgs_sph_laminarFlowAroundSphereBx();
    void sph_laminarFlowAroundSphereBz(const gridreal r[],datareal B[]);
    void setArgs_sph_laminarFlowAroundSphereBz();
    void sph_dipoleB(const gridreal r[],datareal B[]);
    void setArgs_sph_dipoleB();
    void sph_Br(const gridreal r[],datareal B[]);
    void setArgs_sph_Br();
    real sph_B0, sph_R0, sph_d;
#endif
};

void setInitialMagneticField(const gridreal r[], datareal B[]);
void addConstantMagneticField(const gridreal r[], datareal B[]);
void initializeMagneticField();

#endif
