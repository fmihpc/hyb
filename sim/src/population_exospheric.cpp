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
#include "population_exospheric.h"
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
PopulationExospheric::PopulationExospheric(PopulationArgs args) : Population(args)
{
    // Check the necessary information is given
    if(args.R.given == false) {
        ERRORMSG2("exobase R must be given",idStr);
        doabort();
    }
    if(args.distFunc.given == false) {
        ERRORMSG2("distribution function must be given",idStr);
        doabort();
    }
    // Initialize variable values
    R = 0;
    totalRate = 0;
    //distFunc = NULL;
    distFuncId = -1;
}

//! Destructor
PopulationExospheric::~PopulationExospheric() { }

//! Initialize exospheric population
void PopulationExospheric::initialize()
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
    g.prepare_PDF(&distFunc,distFuncId,totalRate);
    updateArgs();
    // Write parameter log
    if(logParams == true && Params::t <= 0) {
        writeLog();
    }
}

//! Create exospheric population particles
void PopulationExospheric::createParticles()
{
    for (int i = 0; i < probround(macroParticlesPerDt); ++i) {
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        newParticle();
#else
        sph_newParticle();
#endif
    }
}

//! Get the neutral density related to the population
real PopulationExospheric::getNeutralDensity(const gridreal r[3])
{
    // assume the first given function is the neutral density
    return distFunc.getValue(r,0);
}

//! Write population log
void PopulationExospheric::writeLog()
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
        populationlog << "% particle population logfile (exospheric)\n";
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
        populationlog << int2string(nn++,2) << ". R [#/m^3]      " << separator;
        populationlog << int2string(nn++,2) << ". totRate [#/s]  " << separator;
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
    populationlog << R << separator;
    populationlog << totalRate << separator;
    populationlog << "\n";
}

//! Update exospheric population arguments
void PopulationExospheric::updateArgs()
{
    Population::updateBaseClassArgs(args);
    // Set exobase
    if(args.R.given == true) {
        if(args.R.value >= 0) {
            this->R = args.R.value;
        } else {
            ERRORMSG2("trying to set R < 0",idStr);
            doabort();
        }
    }
    if(args.totalRate.given == true) {
        // Set total production rate
        if(args.totalRate.value >= 0) {
            this->totalRate = args.totalRate.value;
        } else {
            ERRORMSG2("trying to set totalRate < 0",idStr);
            doabort();
        }
    }
    // Set macroParticlesPerDt
    int exoPopsN = Params::popFactory.getNumberOfPopulations("exospheric");
    if(args.macroParticlesPerDt.given == true) {
        this->macroParticlesPerDt = args.macroParticlesPerDt.value;
    } else {
        if (vth < 1e-9*Params::dx/(1<<Params::currentGridRefinementLevel)/Params::dt) {
            WARNINGMSG2("vth is too small to determine proper macroParticlesPerDt, set now to 100",idStr);
            macroParticlesPerDt = 100;
        } else {
            macroParticlesPerDt = 2*sqrt(3.0)*vth/( Params::dx/(1<<Params::currentGridRefinementLevel) )*Params::dt;
            macroParticlesPerDt *= pi*sqr(R/(Params::dx/(1<<Params::currentGridRefinementLevel)))*Params::macroParticlesPerCell/exoPopsN;
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

//! New exospheric population particle
void PopulationExospheric::newParticle()
{
    gridreal r[3];
    // This loop is executed only once in 99.999.. % of the cases
    do {
        g.generate_random_point(distFuncId,r);
    } while ( sqr(r[0]) + sqr(r[1]) + sqr(r[2]) <= sqr(R) );
    const shortreal x = r[0];
    const shortreal y = r[1];
    const shortreal z = r[2];
    const shortreal vx = vth*gaussrnd();
    const shortreal vy = vth*gaussrnd();
    const shortreal vz = vth*gaussrnd();
    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
}

//! Write distribution function into hc-file
void PopulationExospheric::writeExtraHcFile()
{
    string fileName = "EXTRA_" + idStr + ".hc";
    stringstream ss;
    ss << "exospheric distribution (idStr=" << idStr
       << ",popid=" << popid << "): " << distFunc.toString(""," ");
    string hcHeader = ss.str();
    g.hcwrite_EXTRA(fileName,&distFunc,hcHeader);
}

//! Dump example population config (not implemented)
string PopulationExospheric::configDump()
{
    return "";
}

//! Returns string representation of the population
string PopulationExospheric::toString()
{
    stringstream ss;
    ss << "type = exospheric\n";
    ss << Population::toStringGeneral();
    ss << "total production rate = " << totalRate << " #/s\n";
    ss << "R = " << R/1e3 << " km\n";
    ss << "distribution function = \n{\n" << distFunc.toString(" ","\n") << "}\n";
    ss << "DPDF id = " << distFuncId << "\n";
    return ss.str();
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "newParticle"
void PopulationExospheric::sph_newParticle()
{
    gridreal r[3];
    // This loop is executed only once in 99.999.. % of the cases
    do {
        g.generate_random_point(distFuncId,r);
    } while ( r[0] <= R );
    const shortreal x = r[0];
    const shortreal y = r[1];
    const shortreal z = r[2];
    shortreal vx = vth*gaussrnd();
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();
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

