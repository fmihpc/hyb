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
#include "params.h"
#include "population.h"
#include "population_ionospheric.h"
#include "simulation.h"
#include "random.h"
#include "atmosphere.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Constructor
PopulationIonospheric::PopulationIonospheric(PopulationArgs args) : Population(args)
{
    // Check the necessary information is given
    if(args.R.given == false) {
        ERRORMSG2("radius R must be given",idStr);
        doabort();
    }
    if(args.n.given == true && args.totalRate.given == true) {
        ERRORMSG2("only number density n or total emission rate totalRate can be given",idStr);
        doabort();
    } else if(args.n.given == false && args.totalRate.given == false) {
        ERRORMSG2("number density n or total emission rate totalRate must be given",idStr);
        doabort();
    }
    if(args.distFunc.given == false) {
        ERRORMSG2("distribution function must be given",idStr);
        doabort();
    }

    // Initialize variable values
    n = 0;
    totalRate = 0;
    R = 0;
    //distFunc = NULL;
    distFuncId = -1;
}

//! Destructor
PopulationIonospheric::~PopulationIonospheric() { }

//! Initialize ionospheric population
void PopulationIonospheric::initialize()
{
    if(args.distFunc.name.size() <= 0) {
        ERRORMSG2("distribution function not found",idStr);
        doabort();
    }

    // Set distribution function
    vector<SpatialDistribution> tempFuncs;
    for(unsigned int i = 0; i < args.distFunc.name.size(); ++i) {
        if(args.distFunc.name[i].compare("") == 0) {
            ERRORMSG2("trying to set empty distribution function name",idStr);
            doabort();
        }
        tempFuncs.push_back( SpatialDistribution(args.distFunc.name[i],args.distFunc.funcArgs[i],popid) );
    }
    this->distFunc = MultipleProductDistribution(tempFuncs);

    // Prepare discretized distribution function
    real dummycumsum;
    g.prepare_PDF(&distFunc,distFuncId,dummycumsum);

    // Write parameter log
    if(logParams == true && Params::t <= 0) {
        writeLog();
    }
}

//! Create ionospheric population particles
void PopulationIonospheric::createParticles()
{
    for (int i = 0; i < probround(macroParticlesPerDt); ++i) {
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        newParticle();
#else
        sph_newParticle();
#endif
    }
}

//! Update ionospheric population arguments
void PopulationIonospheric::updateArgs()
{
    Population::updateBaseClassArgs(args);

    // Set radius
    if(args.R.given == true) {
        if(args.R.value >= 0) {
            this->R = args.R.value;
        } else {
            ERRORMSG2("trying to set R < 0",idStr);
            doabort();
        }
    }

    // Set n and totalRate
    if(args.n.given == true && args.totalRate.given == true) {
        ERRORMSG2("only number density n or total emission rate totalRate can be given",idStr);
        doabort();
    } else if(args.n.given == true) {
        if(args.n.value >= 0) {
            this->n = args.n.value;

            // Calculate total rate
            if(n >= 0 && R >= 0 && T >= 0) {
                this->totalRate = n*(4*pi*sqr(R)*sqrt(Params::k_B*T/(2*pi*m)));
            } else {
                ERRORMSG2("cannot calculate totalRate",idStr);
                doabort();
            }
        } else {
            ERRORMSG2("trying to set n < 0",idStr);
            doabort();
        }
    } else if(args.totalRate.given == true) {
        if(args.totalRate.value >= 0) {
            this->totalRate = args.totalRate.value;

            // Calculate number density
            if(R > 0 && T > 0) {
                this->n = totalRate/(4*pi*sqr(R)*sqrt(Params::k_B*T/(2*pi*m)));
            } else if(R == 0 || T == 0) {
                this->n = 0;
            } else {
                ERRORMSG2("cannot calculate n",idStr);
                doabort();
            }
        } else {
            ERRORMSG2("trying to set totalRate < 0",idStr);
            doabort();
        }
    }

    // Set macroParticlesPerDt
    int ionoPopsN = Params::popFactory.getNumberOfPopulations("ionospheric");
    if(args.macroParticlesPerDt.given == true) {
        this->macroParticlesPerDt = args.macroParticlesPerDt.value;
    } else {
        if (vth < 1e-9*Params::dx/(1<<Params::currentGridRefinementLevel)/Params::dt) {
            WARNINGMSG2("vth is too small to determine proper macroParticlesPerDt, set now to 100.",idStr);
            macroParticlesPerDt = 100;
        } else {
            macroParticlesPerDt = sqrt(3.0)*2/pi*vth/( Params::dx/(1<<Params::currentGridRefinementLevel) )*Params::dt;
            macroParticlesPerDt *= pi*sqr(R/(Params::dx/(1<<Params::currentGridRefinementLevel)))*Params::macroParticlesPerCell/ionoPopsN;
            // 2/pi*vth = velocity perpedicular to the surface
        }
    }

    // Calculate macroParticleStatisticalWeight
    if(macroParticlesPerDt <= 0) {
        macroParticleStatisticalWeight = 1;
    } else {
        macroParticleStatisticalWeight = totalRate*Params::dt/macroParticlesPerDt;
    }

    // Write parameter log
    if(logParams == true && Params::t > 0) {
        writeLog();
    }
}

//! New ionospheric population particle
void PopulationIonospheric::newParticle()
{
    gridreal r[3];
    g.generate_random_point(distFuncId,r);

    // Here we project the coordinates of the particle to the sphere of radius r0
    const fastreal r0 = R;
    const fastreal norm = r0/sqrt( sqr(r[0]) + sqr(r[1]) + sqr(r[2]) );

    const shortreal x = r[0]*norm;
    const shortreal y = r[1]*norm;
    const shortreal z = r[2]*norm;
    shortreal vx = vth*gaussrnd();
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();

    // Ensure that speed is upward
    if (vx*x + vy*y + vz*z < 0) {
        vx = -vx;
        vy = -vy;
        vz = -vz;
    }

    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
}

//! Write distribution function into hc-file
void PopulationIonospheric::writeExtraHcFile()
{
    string fileName = "EXTRA_" + idStr + ".hc";
    stringstream ss;
    ss << "ionospheric distribution (idStr=" << idStr
       << ",popid=" << popid << "): " << distFunc.toString(""," ");
    string hcHeader = ss.str();
    g.hcwrite_EXTRA(fileName,&distFunc,hcHeader);
}

//! Write population log
void PopulationIonospheric::writeLog()
{
    // Write particle log header lines
    if(logHeaderWritten == false) {
        // Create logger + logfile
        string logfile = "population_";
        logfile.append(idStr);
        logfile.append(".log");
        populationlog = Logger(logfile.c_str(),100000,false);
        // Initialize log (logfile is created at this point)
        populationlog.init();
        populationlog << "% particle population logfile (ionospheric)\n";
        populationlog << "% idStr = " << idStr << "\n";
        populationlog << "% popid = " <<  static_cast<int>(popid) << "\n";
        populationlog << "% m = " << m << "\n";
        populationlog << "% q = " << q << "\n";
        int nn = 1;
        string separator(1,'\t');
        populationlog << "% ";
        populationlog << int2string(nn++,2) << ". t [s]        " << separator;
        populationlog << int2string(nn++,2) << ". T [K]          " << separator;
        populationlog << int2string(nn++,2) << ". vth [m/s]      " << separator;
        populationlog << int2string(nn++,2) << ". mples/dt [#]   " << separator;
        populationlog << int2string(nn++,2) << ". statWeight [#] " << separator;
        populationlog << int2string(nn++,2) << ". split [-]      " << separator;
        populationlog << int2string(nn++,2) << ". join [-]       " << separator;
        populationlog << int2string(nn++,2) << ". n [#/m^3]      " << separator;
        populationlog << int2string(nn++,2) << ". totRate [#/s]  " << separator;
        populationlog << int2string(nn++,2) << ". R [m]          " << separator;
        populationlog << "\n";
        populationlog << scientific << showpos;
        populationlog.precision(10);
        logHeaderWritten = true;
    }
    string separator(2,' ');
    separator.append(1,'\t');
    populationlog << Params::t << separator;
    populationlog << T << separator;
    populationlog << vth << separator;
    populationlog << macroParticlesPerDt << separator;
    populationlog << macroParticleStatisticalWeight << separator;
    populationlog << static_cast<real>(split) << separator;
    populationlog << static_cast<real>(join) << separator;
    populationlog << n << separator;
    populationlog << totalRate << separator;
    populationlog << R << separator;
    populationlog << "\n";
}

//! Dump example population config (not implemented)
string PopulationIonospheric::configDump()
{
    return "";
}

//! Returns string representation of the population
string PopulationIonospheric::toString()
{
    stringstream ss;
    ss << "type = ionospheric\n";
    ss << Population::toStringGeneral();
    ss << "n(r=R) = " << n/1e6 << " #/cm^3\n";
    ss << "total emission rate = " << totalRate << " #/s\n";
    ss << "R = " << R/1e3  << " km\n";
    ss << "distribution function = \n{\n" << distFunc.toString(" ","\n") << "}\n";
    ss << "DPDF id = " << distFuncId << "\n";
    return ss.str();
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "newParticle"
void PopulationIonospheric::sph_newParticle()
{
    gridreal r[3];
    g.generate_random_point(distFuncId,r);
    // Here we project the coordinates of the particle to the sphere of radius r0
    const fastreal r0 = R;
    //const fastreal norm = r0/sqrt( sqr(r[0]) + sqr(r[1]) + sqr(r[2]) );
    const fastreal norm = r0/r[0];
    const shortreal x = r[0]*norm;
    const shortreal y = r[1];
    const shortreal z = r[2];
    //const shortreal y = r[1]*norm;
    //const shortreal z = r[2]*norm;
    shortreal vx = vth*gaussrnd();
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();
    // Ensure that speed is upward
    if (vx*x + vy*y + vz*z < 0) {
        vx = -vx;
        vy = -vy;
        vz = -vz;
    }
    sph_transf_H2S_R(r);
    gridreal theta = r[1];
    gridreal v[3] = {vx, vy, vz};
    sph_transf_H2S_V(v);
    sph_transf_S2C_V(r,v);
    vx = v[0];
    vy = v[1];
    vz = v[2];
    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*sin(theta),popid);
}

#endif

