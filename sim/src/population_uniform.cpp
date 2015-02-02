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
#include "population_uniform.h"
#include "simulation.h"
#include "random.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Constructor
PopulationUniform::PopulationUniform(PopulationArgs args) : Population(args)
{
    // Initialize variable values
    n = 0;
    V = 0;
    if(args.n.given == false) {
        particleCreationDone = true;
    } else {
        particleCreationDone = false;
    }
}

//! Destructor
PopulationUniform::~PopulationUniform() { }

//! Initialize uniform population
void PopulationUniform::initialize()
{
    // If previous state has been loaded, particles are already created
    if (Params::wsFileGiven) {
        particleCreationDone = true;
    }
    if(logParams == true && Params::t <= 0) {
        writeLog();
    }
}

//! Create uniform population particles
void PopulationUniform::createParticles()
{
    if(particleCreationDone == false) {
        for (int i = 0; i < probround(macroParticlesPerDt); ++i) {
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
            newParticle();
#else
            sph_newParticle();
#endif
        }
        particleCreationDone = true;
    }
}

//! Update uniform population arguments
void PopulationUniform::updateArgs()
{
    Population::updateBaseClassArgs(args);
    // Set number density
    if(args.n.given == true) {
        if(args.n.value >= 0) {
            this->n = args.n.value;
        } else {
            ERRORMSG2("trying to set n < 0",idStr);
            doabort();
        }
    }
    // Set bulk speed
    if(args.V.given == true) {
        if(args.V.value >= 0) {
            this->V = args.V.value;
        } else {
            ERRORMSG2("trying to set V < 0",idStr);
            doabort();
        }
    }
    // Set macroParticlesPerDt
    if(args.macroParticlesPerDt.given == true) {
        this->macroParticlesPerDt = args.macroParticlesPerDt.value;
    } else if(n <= 0) {
        this->macroParticlesPerDt = 0.0;
    } else {
        int uniformPopsN = Params::popFactory.getNumberOfPopulations("uniform");
        macroParticlesPerDt = Params::nx*Params::ny*Params::nz*Params::macroParticlesPerCell/uniformPopsN;
        // usually good idea to have splitting off for initial pops
        // - therefore more particles of initial population should be created
        if(split == true) {
            macroParticlesPerDt *= 1 + 2*cube(Params::currentGridRefinementLevel);
        } else {
            macroParticlesPerDt = int( (1.2 + 4*cube(Params::currentGridRefinementLevel) )*macroParticlesPerDt);
        }
    }
    // Calc macroParticleStatisticalWeight
    if(macroParticlesPerDt <= 0) {
        macroParticleStatisticalWeight = 1;
    }
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    else {
        macroParticleStatisticalWeight = n*Params::box_V/macroParticlesPerDt;
    }
#else
    else {
        real sph_coeff = 5;
        macroParticleStatisticalWeight = n*Params::box_V/macroParticlesPerDt*sph_coeff;
    }
#endif
    // Write parameter log
    if(logParams == true && Params::t > 0) {
        writeLog();
    }
}

//! New uniform population particle
void PopulationUniform::newParticle()
{
    const shortreal x = Params::box_xmin_tight + uniformrnd()*(Params::box_xmax_tight - Params::box_xmin_tight);
    const shortreal y = Params::box_ymin_tight + uniformrnd()*(Params::box_ymax_tight - Params::box_ymin_tight);
    const shortreal z = Params::box_zmin_tight + uniformrnd()*(Params::box_zmax_tight - Params::box_zmin_tight);
    const shortreal vx = -V + vth*gaussrnd();
    const shortreal vy = vth*gaussrnd();
    const shortreal vz = vth*gaussrnd();
    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
}

//! Nothing to write
void PopulationUniform::writeExtraHcFile() { }

//! Write population log
void PopulationUniform::writeLog()
{
    // Write particle log header lines
    if(logHeaderWritten == false) {
        // Create logger + logfile
        string logfile = "population_";
        logfile.append(idStr);
        logfile.append(".log");
        populationlog = Logger(logfile.c_str(),1000000,false);
        // Initialize log (logfile is created at this point)
        populationlog.init();
        populationlog << "% particle population logfile (uniform)\n";
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
        populationlog << int2string(nn++,2) << ". V [m/s]        " << separator;
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
    populationlog << V << separator;
    populationlog << "\n";
}

//! Dump example population config (not implemented)
string PopulationUniform::configDump()
{
    return "";
}

//! Returns string representation of the population
string PopulationUniform::toString()
{
    stringstream ss;
    ss << "type = uniform\n";
    ss << Population::toStringGeneral();
    ss << "n = " << n/1e6 << " #/cm^3\n";
    ss << "V = " << V/1e3 << " km/s\n";
    return ss.str();
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "newParticle"
void PopulationUniform::sph_newParticle()
{
    const shortreal x = Params::box_xmin_tight + uniformrnd()*(Params::box_xmax_tight - Params::box_xmin_tight);
    const shortreal y = Params::box_ymin_tight + uniformrnd()*(Params::box_ymax_tight - Params::box_ymin_tight);
    const shortreal z = Params::box_zmin_tight + uniformrnd()*(Params::box_zmax_tight - Params::box_zmin_tight);
    gridreal r[3] = {x,y,z};
    // Positions transformation from hybrid into spherical coordinates
    sph_transf_H2S_R(r);
    gridreal R = r[0];
    gridreal theta = r[1];
    gridreal R_min = Params::box_xmin;
    //gridreal V_hyb = Params::box_V;
    //gridreal V_sph = Params::sph_box_V;
    shortreal vx = -V + vth*gaussrnd();
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();
    gridreal v[3] = {vx, vy, vz};
    sph_transf_H2S_V(v);
    sph_transf_S2C_V(r,v);
    vx = v[0];
    vy = v[1];
    vz = v[2];
    // we need to include coef sqr(R/R_min)*sin(theta) to keep uniform distribution
    // of uniform population in spherical coordinates. Othewise we will get
    // nonuniform population density: rho_uni = m_uni/V_sph, V_sph = r^2*sin(theta)*dr*dtheta*dphi
    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*sqr(R/R_min)*sin(theta),popid);
    //g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
    //g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*V_sph/V_hyb,popid);
}

#endif

