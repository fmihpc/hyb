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

/* \brief  Resistivity profiles
 *
 * Resitivity has the unit of [Ohm m] = Vs/(Am)m^2/s
 * 
 * The diffusion equation of the magnetic field is:
 *
 * dB/dt = eta/mu_0 * nabla^2(B)
 *
 * Dimensional analysis gives for eta:
 *
 * eta ~ mu_0*dx^2/dt
 *
 * Thus, when changing the grid cell size (dx) or the timestep (dt),
 * magnetic field diffusion changes as dx^2/dt. We parametrize eta in
 * these resitivity profiles as:
 *
 * eta = eta0_c * mu_0 * dx^2/dt,
 *
 * where eta0_c is a dimensionless parameter given in a config file and
 * dx is the base grid cell size and dt is the timestep.
 *
 */

#include <sstream>
#include "resistivity.h"
#include "params.h"
#include "grid.h"
#include "simulation.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Default constructor
ResistivityProfile::ResistivityProfile()
{
    this->ptr = &ResistivityProfile::defaultFunction;
    resetParameters();
}

#define ELSEIF_RESISTIVITY(func) else if(funcName.compare(#func) == 0) { this->ptr = &ResistivityProfile::func; ; setArgs_ ## func(); }

//! Constructor
ResistivityProfile::ResistivityProfile(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &ResistivityProfile::defaultFunction;
    this->args = args;
    resetParameters();
    // RESISTIVITY PROFILES
    if(funcName.compare("") == 0) {
        ERRORMSG("empty resistivity function name");
        doabort();
    }
    ELSEIF_RESISTIVITY(resistivityConstant)
    ELSEIF_RESISTIVITY(resistivityConstantOutsideR)
    ELSEIF_RESISTIVITY(resistivityWithinObstacle)
    ELSEIF_RESISTIVITY(resistivitySpherical)
    ELSEIF_RESISTIVITY(resistivityCartesian)
    ELSEIF_RESISTIVITY(resistivitySphericalPolynomial)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    ELSEIF_RESISTIVITY(sph_resistivityConstantOutsideR)
#endif
    else {
        ERRORMSG2("bad resistivity function name",funcName);
        doabort();
    }
}

//! Destructor
ResistivityProfile::~ResistivityProfile() { }

//! Returns the value at point r
real ResistivityProfile::getValue(const gridreal r[3])
{
    return (this->*ptr)(r);
}

real ResistivityProfile::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

//! Reset private class variables
void ResistivityProfile::resetParameters()
{
    // RESISTIVITY PARAMETERS
    eta = eta0_c = R = R2 = dR_coef = dR = res1 = res2 = x1 = y1 = z1 = x2 = y2 = z2 = 0.0;
    etaArray.clear();
    rMinSqrArray.clear();
    rMaxSqrArray.clear();
    polyCoeffArray.clear();
}

// RESISTIVITY PROFILES

//! Set function arguments
void ResistivityProfile::setArgs_resistivityConstant()
{
    if(args.size() != 1) {
        ERRORMSG2("function takes one argument",name);
        doabort();
    }
    eta0_c = args[0];
    eta = eta0_c*Params::mu_0*(sqr(Params::dx)/Params::dt);
}

//! Constant resistivity
real ResistivityProfile::resistivityConstant(const gridreal r[3])
{
    return eta;
}

//! Set function arguments
void ResistivityProfile::setArgs_resistivityConstantOutsideR()
{
    if(args.size() != 2) {
        ERRORMSG2("function takes two arguments",name);
        doabort();
    }
    eta0_c = args[0];
    eta = eta0_c*Params::mu_0*(sqr(Params::dx)/Params::dt);
    R = args[1];
    R2 = sqr(R);
}

//! Constant resistivity outside the obstacle, zero within
real ResistivityProfile::resistivityConstantOutsideR(const gridreal r[3])
{
    const real r2 = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    if(r2 <= R2) {
        return 0.0;
    }
    return eta;
}

//! Set function arguments
void ResistivityProfile::setArgs_resistivityWithinObstacle()
{
    if(args.size() != 3) {
        ERRORMSG2("function takes three arguments",name);
        doabort();
    }
    eta0_c = args[0];
    eta = eta0_c*Params::mu_0*(sqr(Params::dx)/Params::dt);
    R = args[1];
    dR_coef = args[2];
    R2 = sqr(R);
    dR = dR_coef*Params::dx/real(1 << Params::currentGridRefinementLevel);
}

//! sets resistivity within the obstacle to a set level with a smooth transition
real ResistivityProfile::resistivityWithinObstacle(const gridreal r[3])
{
    const real rr = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    if (rr > sqr(R + dR)) {
        return 0.0;
    } else if (rr < sqr(R - dR)) {
        return eta;
    } else {
        return eta * (R + dR - sqrt(rr))/(2*dR);
    }
}

//! Set function arguments
void ResistivityProfile::setArgs_resistivitySpherical()
{
    // Number of function arguments
    unsigned int nArgs = this->args.size();
    // FUNCTION ARGUMENTS
    if(nArgs % 3 != 0 || nArgs == 0) {
        ERRORMSG2("function takes N*3 arguments, N = 1,2..",name);
        doabort();
    }
    // Put argumets into an array
    for(unsigned int i=0; i < nArgs; i = i+3) {
        etaArray.push_back(args[i]*Params::mu_0*(sqr(Params::dx)/Params::dt));
        rMinSqrArray.push_back(sqr(args[i+1]));
        rMaxSqrArray.push_back(sqr(args[i+2]));
    }
}

//! Spherically symmetric resistivity shells
real ResistivityProfile::resistivitySpherical(const gridreal r[3])
{
    const real rr = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    for(unsigned int i = 0; i < etaArray.size(); i++) {
        if(rr > rMinSqrArray[i] && rr < rMaxSqrArray[i]) {
            return etaArray[i];
        }
    }
    return 0.0;
}

//! Set function arguments
void ResistivityProfile::setArgs_resistivityCartesian()
{
    // Number of function arguments
    unsigned int nArgs = this->args.size();
    // FUNCTION ARGUMENTS
    if(nArgs != 8) {
        ERRORMSG2("function takes 8 arguments, N = 1,2..",name);
        doabort();
    }
    res1 = args[0]*Params::mu_0*(sqr(Params::dx)/Params::dt);
    res2 = args[1]*Params::mu_0*(sqr(Params::dx)/Params::dt);
    x1 = args[2];
    y1 = args[3];
    z1 = args[4];
    x2 = args[5];
    y2 = args[6];
    z2 = args[7];
}

//! Cartesian resistivity for reconnection
real ResistivityProfile::resistivityCartesian(const gridreal r[3])
{
    static real rd = res1 - res2;
    real res_xyz = min(res1 - rd/(x2-x1)*(fabs(r[0])-x1), res1 - rd/(y2-y1)*(fabs(r[1])-y1));
    res_xyz = min( res_xyz, res1 - rd/(z2-z1)*(fabs(r[2])-z1) );
    res_xyz = max(res2, res_xyz);
    res_xyz = min(res1, res_xyz);
    return res_xyz;
}

//! Set function arguments
void ResistivityProfile::setArgs_resistivitySphericalPolynomial()
{
    // Number of function arguments
    unsigned int nArgs = this->args.size();
    // FUNCTION ARGUMENTS
    if(nArgs < 1) {
        ERRORMSG2("function takes >0 arguments",name);
        doabort();
    }
    // Put argumets into an array
    for(unsigned int i=0; i < nArgs; ++i) {
        polyCoeffArray.push_back(args[i]*Params::mu_0*(sqr(Params::dx)/Params::dt));
    }
}

//! Spherically symmetric polynomial resistivity
real ResistivityProfile::resistivitySphericalPolynomial(const gridreal r[3])
{
    const real r2 = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    if(r2 > Params::R_zeroFields2)  {
        return 0.0;
    }
    const real absr = sqrt(r2);
    real resultEta = 0.0;
    for(unsigned int i = 0; i < polyCoeffArray.size(); i++) {
        resultEta += polyCoeffArray[i]*pow(absr,static_cast<real>(i));
    }
    return resultEta;
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Set function arguments
void ResistivityProfile::setArgs_sph_resistivityConstantOutsideR()
{
    if(args.size() != 2) {
        ERRORMSG2("function takes two arguments",name);
        doabort();
    }
    eta0_c = args[0];
    eta = eta0_c*Params::mu_0*(sqr(Params::dx)/Params::dt);
    R = args[1];
}

//! (SPHERICAL) Spherical version of "ionoCosSzaDayConstantNight"
real ResistivityProfile::sph_resistivityConstantOutsideR(const gridreal r[3])
{
    const real rr = r[0];
    if(rr <= R) {
        return 0.0;
    }
    return eta;
}

#endif

//! Set resistivity
void setResistivity()
{
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("resistivityFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false || tempNames.size() < 1 || tempArgs.size() < 1) {
        if(Params::resistivityFunction.isDefined() == false) {
            mainlog << "none\n";
        }
        return;
    }
    if(tempNames.size() > 1 || tempArgs.size() > 1) {
        WARNINGMSG("multiple resistivity functions defined, using only the first one");
    }
    if(Params::resistivityFunction.isDefined() == false) {
        Params::resistivityFunction = ResistivityProfile(tempNames[0],tempArgs[0]);
        mainlog << "RESISTIVITY FUNCTION: ";
    } else {
        ResistivityProfile newResFunc = ResistivityProfile(tempNames[0],tempArgs[0]);
        if(newResFunc.toString().compare(Params::resistivityFunction.toString()) == 0) {
            // same function with same parameters, no update
            return;
        }
        Params::resistivityFunction = newResFunc;
        mainlog << "RESISTIVITY FUNCTION UPDATE: ";
    }
    mainlog << tempNames[0] << " ";
    for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
        mainlog << tempArgs[0][i] << " ";
    }
    mainlog << "\n";
    // Set resistivity in the grid
    g.set_resistivity(Params::resistivityFunction);
    // Save extra hc-file
    if(Params::saveExtraHcFiles == true) {
        string fn = "EXTRA_resistivity_" + Params::getSimuTimeStr() + ".hc";
        string hcHeader = "resistivity function: " + Params::resistivityFunction.toString();
        g.hcwrite_EXTRA(fn,&Params::resistivityFunction,hcHeader);
    }
}

//! Initialize resisitivity
void initializeResistivity()
{
    MSGFUNCTIONCALL("initializeResisitity");
    setResistivity();
    MSGFUNCTIONEND("initializeResistivity");
}

