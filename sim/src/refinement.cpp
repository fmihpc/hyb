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
#include "refinement.h"
#include "params.h"
#include "grid.h"
#include "simulation.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Constructor
GridRefinementProfile::GridRefinementProfile()
{
    this->ptr = &GridRefinementProfile::defaultFunction;
}

#define ELSEIF_REFINEMENT(func) else if(funcName.compare(#func) == 0) { this->ptr = &GridRefinementProfile::func; }

//! Constructor
GridRefinementProfile::GridRefinementProfile(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &GridRefinementProfile::defaultFunction;
    this->args = args;
    // GRID REFINEMENT PROFILES
    if(funcName.compare("") == 0) {
        ERRORMSG("empty grid refinement function name");
        doabort();
    }
    ELSEIF_REFINEMENT(refineSpherical)
    ELSEIF_REFINEMENT(refineCartesian)
    ELSEIF_REFINEMENT(refineSphericalAndCartesian)
    ELSEIF_REFINEMENT(refineTitan)
    else {
        ERRORMSG2("bad grid refinement function name",funcName);
        doabort();
    }
}

//! Returns the value at point r
GridRefinementProfile::~GridRefinementProfile() { }

real GridRefinementProfile::getValue(const gridreal r[3])
{
    return (this->*ptr)(r);
}

//! Default function, which aborts the program if called
gridreal GridRefinementProfile::defaultFunction(const gridreal r[3])
{
    ERRORMSG("function pointer not set");
    doabort();
    return -1;
}

// GRID REFINEMENT PROFILES

//! Spherical grid refinements - config file arguments: R1min R1max R2min R2max ...
gridreal GridRefinementProfile::refineSpherical(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() % 2 != 0 || this->args.size() == 0) {
        ERRORMSG2("function takes even amount arguments",this->name);
        doabort();
    }
    // Number of function arguments
    unsigned int nArgs = args.size();
    real rr = sqr(r[0]) + sqr(r[1]) + sqr(r[2]);
    real factor = 1.1;
    for(unsigned int i=0; i < nArgs; i=i+2) {
        real rMin = args[i];
        real rMax = args[i+1];
        if( rr > sqr(rMin)  && rr < sqr(rMax) ) {
            factor /= 2.0;
        }
    }
    return Params::dx*factor;
}

//! Cartesian grid refinements - config file arguments: xmin1 xmax1 ymin1 ymax1 zmin1 zmax1 xmin2 xmax2 ymin2 ymax2 zmin2 zmax2 ...
gridreal GridRefinementProfile::refineCartesian(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() % 6 != 0 || this->args.size() == 0) {
        ERRORMSG2("function takes N*6 arguments, N = 1,2..",this->name);
        doabort();
    }
    // Number of function arguments
    unsigned int nArgs = args.size();
    const real x = r[0];
    const real y = r[1];
    const real z = r[2];
    real factor = 1.1;
    for(unsigned int i=0; i < nArgs; i=i+6) {
        real xmin = args[i];
        real xmax = args[i+1];
        real ymin = args[i+2];
        real ymax = args[i+3];
        real zmin = args[i+4];
        real zmax = args[i+5];
        if(x >= xmin  && x <= xmax &&
	   y >= ymin  && y <= ymax &&
	   z >= zmin  && z <= zmax) {
            factor /= 2.0;
        }
    }
    return Params::dx*factor;
}

//! Spherical and Cartesian grid refinements - config file arguments: S1 R1min C1 xmin1 xmax1 ymin1 ymax1 zmin1 zmax1 C2 xmin2 xmax2 ymin2 ymax2 zmin2 zmax2 ...
gridreal GridRefinementProfile::refineSphericalAndCartesian(const gridreal r[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() < 2) {
        ERRORMSG2("function takes at least two arguments",this->name);
        doabort();
    }
    // Number of function arguments
    unsigned int nArgs = args.size();
    const real x = r[0];
    const real y = r[1];
    const real z = r[2];
    const real r2 = sqr(x) + sqr(y) + sqr(z);
    real factor = 1.1;
    for(unsigned int i=0; i < nArgs; i++) {
        const real C = args[i];
        //! Spherical refinement
        if(C == 1) {
            const real rMin2 = sqr(args[i+1]);
            const real rMax2 = sqr(args[i+2]);
            if( r2 > rMin2  && r2 < rMax2 ) {
                factor /= 2.0;
            }
            i += 2;
        }
        //! Cartesian refinement
        else if(C == 2) {
            real xmin = args[i+1];
            real xmax = args[i+2];
            real ymin = args[i+3];
            real ymax = args[i+4];
            real zmin = args[i+5];
            real zmax = args[i+6];
            if(x >= xmin  && x <= xmax &&
               y >= ymin  && y <= ymax &&
               z >= zmin  && z <= zmax) {
                factor /= 2.0;
            }
            i += 6;
        }
        else {
            ERRORMSG2("unrecognized refinement type, use 1 for spherical and 2 for Cartesian refinements",this->name);
            doabort();
        }
    }
    return Params::dx*factor;
}

//! Grid refinment for Titan
gridreal GridRefinementProfile::refineTitan(const gridreal x[3])
{
    // FUNCTION ARGUMENTS
    if(this->args.size() != 3) {
        ERRORMSG2("function takes three arguments",this->name);
        doabort();
    }
    real D1 = args[0]; //radius (in R_P) for inner half-sphere and cylinder
    real D2 = args[1]; //scale for Rectangular part
    real D3 = args[2]; //front
    const gridreal x_1 = Params::box_xmin, x_2 = Params::box_xmax,
                   y_1 = Params::box_ymin, y_2 = Params::box_ymax,
                   z_1 = Params::box_zmin, z_2 = Params::box_zmax;
    gridreal g_size  = 1.1;
    //Refinement at the distance of dG from the edges.
    gridreal dG = 1.4 * Params::dx;
    if ((x[0]<(x_2-dG)) && (x[0]>(x_1+dG))
        && (x[1]<(y_2-dG)) && (x[1]>(y_1+dG))
        && (x[2]<(z_2-dG)) && (x[2]>(z_1+dG))) {
        g_size = 0.55;
        // Rectangular refinement
        const gridreal dG2 = D2*0.7*Params::R_P, dG3 = D2*1.1*Params::R_P, dGx = D3*0.4*Params::R_P;
        if ( (x[0]<(x_2-dGx)) && (x[0]>(x_1+dG2))
             && (x[1]<(y_2-dG3)) && (x[1]>(y_1+dG3))
             && (x[2]<(z_2-dG3)) && (x[2]>(z_1+dG3)) ) {
            g_size = 0.30;
            // Refinement in shape of a half-sphere and a cylinder
            if ( (x[0]<(D1-1)*Params::R_P) && (x[0]>(-D1)*Params::R_P) ) { //in third this was abs(x[0])<3.0*R_P
                const gridreal yz2=sqr(x[1])+sqr(x[2]);
                const gridreal Rad2=sqr(D1*Params::R_P);
                if (x[0]<-Params::R_P) {
                    if (yz2<Rad2) {
                        g_size = 0.15;
                    }
                } else {
                    if (yz2+sqr(x[0]+Params::R_P)<Rad2) {
                        g_size = 0.15;
                    }
                }
            }
        }
    }
    return g_size*Params::dx;
}

//! Initialize grid refinement functions
void initializeGridRefinement()
{
    MSGFUNCTIONCALL("initializeGridRefinement");
    mainlog << "GRID REFINEMENT FUNCTION: ";
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("gridRefinementFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false) {
        mainlog << "none\n";
    } else if (tempNames.size() != 1 || tempArgs.size() != 1) {
        ERRORMSG("only one grid refinement function can be defined");
        doabort();
    } else {
        mainlog << tempNames[0] << " ";
        for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
            mainlog << tempArgs[0][i] << " ";
        }
        mainlog << "\n";
        Params::gridRefinementFunction = GridRefinementProfile(tempNames[0],tempArgs[0]);
        // Refine the grid
        g.Refine(Params::gridRefinementFunction);
        // Save extra hc-file
        if(Params::saveExtraHcFiles == true) {
            string fileName = "EXTRA_gridRefinement.hc";
            string hcHeader = "grid refinement function: " + Params::gridRefinementFunction.toString();
            g.hcwrite_EXTRA(fileName,&Params::gridRefinementFunction,hcHeader);
        }
    }
    MSGFUNCTIONEND("initializeGridRefinement");
}

