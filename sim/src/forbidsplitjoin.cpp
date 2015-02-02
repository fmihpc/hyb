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
#include "forbidsplitjoin.h"
#include "params.h"
#include "grid.h"
#include "simulation.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Default constructor
ForbidSplitAndJoinProfile::ForbidSplitAndJoinProfile()
{
    this->ptr = &ForbidSplitAndJoinProfile::defaultFunction;
}

#define ELSEIF_SPLITJOIN(func) else if(funcName.compare(#func) == 0) { this->ptr = &ForbidSplitAndJoinProfile::func; }

//! Constructor
ForbidSplitAndJoinProfile::ForbidSplitAndJoinProfile(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &ForbidSplitAndJoinProfile::defaultFunction;
    this->args = args;
    // FORBID SPLITANDJOIN PROFILES
    if(funcName.compare("") == 0) {
        ERRORMSG("empty resistivity function name");
        doabort();
    }
    ELSEIF_SPLITJOIN(forbidSplitAndJoinDefault)
    ELSEIF_SPLITJOIN(forbidSplitAndJoinInsideSphere)
    ELSEIF_SPLITJOIN(forbidSplitAndJoinOutsideSphere)
    ELSEIF_SPLITJOIN(forbidSplitAndJoinInsideCuboid)
    ELSEIF_SPLITJOIN(forbidSplitAndJoinOutsideCuboid)
    ELSEIF_SPLITJOIN(forbidSplitAndJoinTitan)
    else {
        ERRORMSG2("bad resistivity function name",funcName);
        doabort();
    }
}

//! Destructor
ForbidSplitAndJoinProfile::~ForbidSplitAndJoinProfile() { }

//! Returns the value at point r
real ForbidSplitAndJoinProfile::getValue(const gridreal r[3])
{
    return (this->*ptr)(r);
}

//! Default function, which aborts the program if called
bool ForbidSplitAndJoinProfile::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

// FORBID SPLITANDJOIN PROFILES

//! Default function for (spatially) forbidding split and join
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinDefault(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 0) {
        ERRORMSG2("function takes no arguments",this->name);
        doabort();
    }
    return false;
}

//! Forbid split&join inside a sphere of radius r
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinInsideSphere(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes one argument",this->name);
        doabort();
    }
    gridreal r2 = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    gridreal R2 = sqr(args[0]);
    if ( r2 <=  R2 ) {
        return true;
    }
    return false;
}

//! Forbid split&join outside a sphere of radius r
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinOutsideSphere(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 1) {
        ERRORMSG2("function takes one argument",this->name);
        doabort();
    }
    gridreal r2 = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    gridreal R2 = sqr(args[0]);
    if ( r2 >  R2 ) {
        return true;
    }
    return false;
}

//! Forbid split&join inside of a cuboid xmin xmax ymin ymax zmin zmax
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinInsideCuboid(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 6) {
        ERRORMSG2("function takes six arguments",this->name);
        doabort();
    }
    const gridreal x = r[0];
    const gridreal y = r[1];
    const gridreal z = r[2];
    const gridreal xmin = args[0];
    const gridreal xmax = args[1];
    const gridreal ymin = args[2];
    const gridreal ymax = args[3];
    const gridreal zmin = args[4];
    const gridreal zmax = args[5];
    if (x >= xmin  && x <= xmax &&
        y >= ymin  && y <= ymax &&
        z >= zmin  && z <= zmax) {
        return true;
    }
    return false;
}

//! Forbid split&join outside of a cuboid xmin xmax ymin ymax zmin zmax
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinOutsideCuboid(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 6) {
        ERRORMSG2("function takes six arguments",this->name);
        doabort();
    }
    const gridreal x = r[0];
    const gridreal y = r[1];
    const gridreal z = r[2];
    const gridreal xmin = args[0];
    const gridreal xmax = args[1];
    const gridreal ymin = args[2];
    const gridreal ymax = args[3];
    const gridreal zmin = args[4];
    const gridreal zmax = args[5];
    if (x < xmin || x > xmax ||
        y < ymin || y > ymax ||
        z < zmin || z > zmax) {
        return true;
    }
    return false;
}

//! Forbid split&join profile for Titan
bool ForbidSplitAndJoinProfile::forbidSplitAndJoinTitan(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 2) {
        ERRORMSG2("function takes two arguments",this->name);
        doabort();
    }
    real D2 = args[0]; //scale for Rectangular part
    real D3 = args[1]; //front
    if ( Params::box_xmax - r[0] > D3*0.4*Params::R_P ) {
        return false;
    }
    if ( abs(r[1]) < (Params::box_Y/2.0 - D2*0.75*Params::R_P) ||
         abs(r[2]) < (Params::box_Z/2.0 - D2*0.75*Params::R_P) ) {
        return true;
    } else {
        return false;
    }
}

//! Initialize forbid split&join
void initializeForbidSplitJoin()
{
    MSGFUNCTIONCALL("initializeForbidSplitJoin");
    mainlog << "FORBID SPLIT&JOIN FUNCTION: ";
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("forbidSplitAndJoinFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false) {
        mainlog << "none\n";
    } else if (tempNames.size() != 1 || tempArgs.size() != 1) {
        ERRORMSG("only one forbid split&join function can be defined");
        doabort();
    } else {
        mainlog << tempNames[0] << " ";
        for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
            mainlog << tempArgs[0][i] << " ";
        }
        mainlog << "\n";
        Params::forbidSplitAndJoinFunction = ForbidSplitAndJoinProfile(tempNames[0],tempArgs[0]);
        // Set forbid split&join flags into the grid (cells)
        g.forbid_split_and_join(Params::forbidSplitAndJoinFunction);
        // Save extra hc-file
        if(Params::saveExtraHcFiles == true) {
            string fileName = "EXTRA_forbidSplitAndJoin.hc";
            string hcHeader = "forbid split&join function: " + Params::forbidSplitAndJoinFunction.toString();
            g.hcwrite_EXTRA(fileName,&Params::forbidSplitAndJoinFunction,hcHeader);
        }
    }
    MSGFUNCTIONEND("initializeForbidSplitJoin");
}

