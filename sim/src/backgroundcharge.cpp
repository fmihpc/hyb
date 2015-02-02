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

#include <sstream>
#include "backgroundcharge.h"
#include "params.h"
#include "grid.h"
#include "simulation.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Default constructor
BackgroundChargeDensityProfile::BackgroundChargeDensityProfile()
{
    // Make sure the class is not used without
    // proper initialization (this aborts execution
    // if getValue called).
    this->ptr = &BackgroundChargeDensityProfile::defaultFunction;
}

#define ELSEIF_BGCHARGEDENS_FUNC(func) else if(funcName.compare(#func) == 0) { this->ptr = &BackgroundChargeDensityProfile::func; }

//! Constructor
BackgroundChargeDensityProfile::BackgroundChargeDensityProfile(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &BackgroundChargeDensityProfile::defaultFunction;
    this->args = args;
    // BACKGROUND CHARGE DENSITY PROFILES
    if(funcName.compare("") == 0) {
        ERRORMSG("empty background charge density function name");
        doabort();
    }
    ELSEIF_BGCHARGEDENS_FUNC(smoothObstacle)
    ELSEIF_BGCHARGEDENS_FUNC(sharpObstacle)
    ELSEIF_BGCHARGEDENS_FUNC(bgCosSzaDayConstantNight)
    ELSEIF_BGCHARGEDENS_FUNC(bgTitan)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    ELSEIF_BGCHARGEDENS_FUNC(sph_smoothObstacle)
#endif
    else {
        ERRORMSG2("bad background charge density function name",funcName);
        doabort();
    }
}

BackgroundChargeDensityProfile::~BackgroundChargeDensityProfile() { }

//! Returns the density value at point r
real BackgroundChargeDensityProfile::getValue(const gridreal r[])
{
    return (this->*ptr)(r);
}

//! Default function, which aborts the program if called
real BackgroundChargeDensityProfile::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

// BACKGROUND CHARGE DENSITY PROFILES

//! Smooth obstacle: n0*exp(-(r-R)/r0)
real BackgroundChargeDensityProfile::smoothObstacle(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 3) {
        ERRORMSG2("function takes three arguments",this->name);
        doabort();
    }
    const real R = args[0];
    const real r0 = args[1];
    const real n0 = args[2];
    const fastreal absr = normvec(r);
    if(absr < R) {
        return n0;
    } else {
        return n0*exp(-(absr-R)/r0);
    }
}
//! Sharp obstacle: n0*exp(-(r-R)^2/r0^2)
real BackgroundChargeDensityProfile::sharpObstacle(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 3) {
        ERRORMSG2("function takes three arguments",this->name);
        doabort();
    }
    const real R = args[0];
    const real r0 = args[1];
    const real n0 = args[2];
    const fastreal absr = normvec(r)-R;
    if(absr < 0.) {
        return n0;
    } else {
        return n0*exp(-absr*absr/r0/r0);
    }
}
//! Smooth obstacle: Goes from SZA=0 to 90 as cos(SZA) and SZA>90 is constant
real BackgroundChargeDensityProfile::bgCosSzaDayConstantNight(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 5) {
        ERRORMSG2("function takes five arguments",this->name);
        doabort();
    }
    const real R = args[0];
    const real r0 = args[1];
    const real n0 = args[2];
    const real daySideFactor = args[3];
    const real nightSideFactor = args[4];
    if (nightSideFactor < 0) {
        ERRORMSG2("nightSideFactor < 0",this->name);
        doabort();
    }
    if (daySideFactor < 0) {
        ERRORMSG2("daySideFactor < 0",this->name);
        doabort();
    }
    const fastreal absr = normvec(r);
    const fastreal sza = acos(r[0]/absr);
    // No ionosphere inside the radius
    //   if(rr < R) { return 0.0; }
    fastreal effdens;
    if(sza < pi/2) {
        effdens = daySideFactor + (nightSideFactor - daySideFactor)*(1-cos(sza));
    } else {
        effdens = nightSideFactor;
    }
    if(absr < R) {
        return effdens*n0;
    } else {
        return effdens*n0*exp(-(absr-R)/r0);
    }
}
//! Bg function for Titan giving SZA dependance for NeMax and for Hmax and the altitude profile
real BackgroundChargeDensityProfile::bgTitan(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 2) {
        ERRORMSG2("function takes two arguments",this->name);
        doabort();
    }
    const real hscale_high = args[0]; //scale hight for above H_max (200km) [Withouth MP 4e5m]
    const real hscale_low = args[1]; //scale hight for bringing Ne to zero below H_max (600km)
    const real Ne_amb = 0; //Ne for high altitudes 0 [Withouth MP 200cm-3*e]
    const real cosSZA=CosSZA(r);
    const real SZAdeg=acos(cosSZA)*180/pi;
    //H_max and Ne_max SZA dependance from K.Agren's RPWS values Feb 2009
    const real ne0 = 5.6e8*Params::e; //Ne_max at night side (SZA>SZAboundary+delta)
    const real ne1 = 2.3e9*Params::e; //Ne_max at day side is ne1+ne0 (SZA<SZAboundary-delta)
    const real SZAboundary = 80;
    const real SZAdelta = 34;
    const real H_max = 950*SZAdeg + 1.03e6; //small dependance on SZA
    real Ne_max=ne0;
    if (SZAdeg<=SZAboundary-SZAdelta) {
        Ne_max+= ne1;
    } else if (SZAdeg<SZAboundary+SZAdelta) {
        Ne_max+= ne1*0.5*(1+cos(pi*0.5*(SZAdeg-SZAboundary+SZAdelta)/SZAdelta));
    }
    // const real Ne_max = (Ne_max-Ne_amb)*0.5*(1.0+cosSZA) + Ne_amb; // old
    // altitude dependance
    const real H = normvec(r) - Params::R_P; //altitude
    real dens;
    if (H>=H_max) {
        dens = (Ne_max-Ne_amb) * exp((H_max-H)/hscale_high) + Ne_amb; //ambient density not included
    } else if (H > H_max - hscale_low) {
        dens = Ne_max * 0.5*(1.0 + cos((H-H_max)/hscale_low *pi));
    } else {
        dens=0; //density inside
    }
    return (dens<0)? 0 : dens;  // additional check for not returning negative values
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "Smooth obstacle"
real BackgroundChargeDensityProfile::sph_smoothObstacle(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 3) {
        ERRORMSG2("function takes three arguments",this->name);
        doabort();
    }
    const real R = args[0];
    const real r0 = args[1];
    const real n0 = args[2];
    const fastreal absr = normvec(r);
    if(r[0] < R) {
        return n0;
    } else {
        return n0*exp(-(r[0]-R)/r0);
    }
}

#endif

//! Initialize background charge density
void initializeBackgroundChargeDensity()
{
    MSGFUNCTIONCALL("initializeBackgroundChargeDensity");
    // Set background charge density
    mainlog << "BACKGROUND CHARGE DENSITY FUNCTION: ";
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("bgChargeDensityFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false) {
        mainlog << "none\n";
    } else if (tempNames.size() != 1 || tempArgs.size() != 1) {
        ERRORMSG("only one background charge density function can be defined");
        doabort();
    } else {
        mainlog << tempNames[0] << " ";
        for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
            mainlog << tempArgs[0][i] << " ";
        }
        mainlog << "\n";
        Params::bgChargeDensityFunction = BackgroundChargeDensityProfile(tempNames[0],tempArgs[0]);
        // Set background charge density into the grid (cells)
        g.set_bgRhoQ(Params::bgChargeDensityFunction);
        // Save extra hc-file
        if(Params::saveExtraHcFiles == true) {
            string fileName = "EXTRA_backgroundChargeDensity.hc";
            string hcHeader = "background charge density function: " + Params::bgChargeDensityFunction.toString();
            g.hcwrite_EXTRA(fileName,&Params::bgChargeDensityFunction,hcHeader);
        }
    }
    MSGFUNCTIONEND("initializeBackgroundChargeDensity");
}
