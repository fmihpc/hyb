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
#include "population_solarwind.h"
#include "simulation.h"
#include "random.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Constructor
PopulationSolarWind::PopulationSolarWind(PopulationArgs args) : Population(args)
{
    // Check the necessary information is given
    if(args.n.given == false) {
        ERRORMSG2("number density must be given",idStr);
        doabort();
    }
    if(args.V.given == false) {
        ERRORMSG2("bulk speed must be given",idStr);
        doabort();
    }
    // Initialize variable values
    n = 0;
    V = 0;
    backWallWeight = 0;
    negativeV = false;
}

//! Destructor
PopulationSolarWind::~PopulationSolarWind() { }

//! Initialize solarwind population
void PopulationSolarWind::initialize()
{
    if(logParams == true && Params::t <= 0) {
        writeLog();
    }
}

//! Create solar wind population particles
void PopulationSolarWind::createParticles()
{
    for (int i = 0; i < probround(macroParticlesPerDt); ++i) {
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        newParticle();
#else
        sph_newParticle();
#endif
    }
}

//! Example addParticle function, which can be called from the main code
void PopulationSolarWind::addParticle(shortreal x,shortreal y,shortreal z,real w)
{
    shortreal vx = -vth*derivgaussrnd(V/vth);
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();
    g.addparticle(x,y,z,vx,vy,vz,w,popid);
}

//! Update solarwind population arguments
void PopulationSolarWind::updateArgs()
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
            negativeV = false;
        } else {
            this->V = -args.V.value;
            negativeV = true;
        }
    }
    // Set backwall weight
    if(args.backWallWeight.given == true) {
        if(args.backWallWeight.value >= 0) {
            this->backWallWeight = args.backWallWeight.value;
        } else {
            ERRORMSG2("trying to set backWallWeight < 0",idStr);
            doabort();
        }
    }
    // Set macroParticlesPerDt
    int swPopsN = Params::popFactory.getNumberOfPopulations("solarwind");
    if(args.macroParticlesPerDt.given == true) {
        this->macroParticlesPerDt = args.macroParticlesPerDt.value;
    } else {
        macroParticlesPerDt = Params::nx*Params::ny*Params::nz*Params::macroParticlesPerCell/swPopsN;
        macroParticlesPerDt *= Params::dt/((Params::box_xmax-Params::box_xmin)/V);
    }
    // Set macroParticleStatisticalWeight
    if(macroParticlesPerDt <= 0) {
        macroParticleStatisticalWeight = 1;
    } else {
        real macroParticlesInFreeFlow = Params::nx*Params::ny*Params::nz*Params::macroParticlesPerCell;
        real temp_sw_macros = macroParticlesInFreeFlow/swPopsN * Params::dt/(Params::box_X/V);
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        macroParticleStatisticalWeight = n*Params::box_V*swPopsN/(macroParticlesInFreeFlow);
#else
        const real sph_coeff = 2*1.5;
        macroParticleStatisticalWeight = n*Params::sph_box_V*swPopsN/(macroParticlesInFreeFlow)*sph_coeff;
#endif
        // if macroParticlePerDt given in config
        if(macroParticlesPerDt != temp_sw_macros) {
            macroParticleStatisticalWeight *= temp_sw_macros/macroParticlesPerDt;
        }
    }
    // Write parameter log
    if(logParams == true && Params::t > 0) {
        writeLog();
    }
}

//! New solar wind population particle
void PopulationSolarWind::newParticle()
{
    shortreal y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
    shortreal z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
    shortreal x,vx;
    if(negativeV == false) {
        // from the front wall
        x = Params::box_xmax_tight - V*Params::dt;
        vx = -vth*derivgaussrnd(V/vth);
    } else {
        // from the back wall
        x = Params::box_xmin_tight + V*Params::dt;
        vx = vth*derivgaussrnd(V/vth);
    }
    shortreal vy = vth*gaussrnd();
    shortreal vz = vth*gaussrnd();
    g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
    // Back wall flow for cases with high thermal velocity.
    if (backWallWeight > 0) {
        // uniform distribution
        x = Params::box_xmin_tight + 0.05*Params::dx;
        y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
        z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
        // vx > 0
        vx = vth*derivgaussrnd(-V/vth);
        vy = vth*gaussrnd();
        vz = vth*gaussrnd();
        g.addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*backWallWeight,popid);
    }
}

//! Nothing to write
void PopulationSolarWind::writeExtraHcFile() { }

//! Write population log
void PopulationSolarWind::writeLog()
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
        populationlog << "% particle population logfile (solarwind)\n";
        populationlog << "% idStr = " << idStr << "\n";
        populationlog << "% popid = " <<  static_cast<int>(popid) << "\n";
        populationlog << "% m = " << m << "\n";
        populationlog << "% q = " << q << "\n";
        // File structure
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
        populationlog << int2string(nn++,2) << ". backWall [-]   " << separator;
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
    populationlog << backWallWeight << separator;
    populationlog << "\n";
}

//! Dump example population config (not implemented)
string PopulationSolarWind::configDump()
{
    return "";
}

//! Returns string representation of the population
string PopulationSolarWind::toString()
{
    stringstream ss;
    if(negativeV == false) {
        ss << "type = solarwind\n";
    } else {
        ss << "type = solarwind (from the back wall)\n";
    }
    ss << Population::toStringGeneral();
    ss << "n = " << n/1e6 << " #/cm^3\n";
    if(negativeV == false) {
        ss << "V = " << V/1e3 << " km/s\n";
    } else {
        ss << "V = " << -V/1e3 << " km/s\n";
    }
    ss << "backWallWeight = " << backWallWeight << "\n";
    return ss.str();
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "newParticle". Use For Radial expanding and collapsing. Velocity jump.
/*void PopulationSolarWind::sph_newParticle()
{
	shortreal y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
	shortreal z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
	shortreal x,vx;
	if(negativeV == false) {
		// from the front wall
		x = Params::box_xmax_tight - V*Params::dt;
		vx = -vth*derivgaussrnd(V/vth);
	} else {
		// from the back wall
		x = Params::box_xmin_tight + V*Params::dt;
		vx = vth*derivgaussrnd(V/vth);
	}
	shortreal vy = vth*gaussrnd();
	shortreal vz = vth*gaussrnd();
	// This block is created for velocity jump
	real t = Params::t;   // Current time
	//real t_start  = 220;  // Start of velocity jump (small box)
	//real t_finish = 260;  // Finish of velocity jump (small box)
	//real t_start  = 1.6e6;           // Start of velocity jump (large box)
	//real t_finish = 1.6e6 + 0.04e6;  // Finish of velocity jump (large box)
	real t_start  = 2.0e7;           // Start of velocity jump (extra large box)
	real t_finish = 2.0e7 + 0.04e7;  // Finish of velocity jump (extra large box)
	if (t >= t_start && t <= t_finish) {
		vx *= 2.0;
		g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
	} else {
		g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight,popid);
	}
	// Back wall flow for cases with high thermal velocity.
	if (backWallWeight > 0) {
		// uniform distribution
		x = Params::box_xmin_tight + 0.05*Params::dx;
		y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
		z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
		// vx > 0
		vx = vth*derivgaussrnd(-V/vth);
		vy = vth*gaussrnd();
		vz = vth*gaussrnd();
		g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*backWallWeight,popid);
	}
}*/

//! (SPHERICAL) Spherical version of "newParticle". Use For Flat front propagation
void PopulationSolarWind::sph_newParticle()
{
    // Spehrical propagation: initial velocities are along of the one of spherical axis
    if(Params::sph_propagation_type == 0) {
        shortreal y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
        shortreal z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
        shortreal x=0.0,vx;
        gridreal r[3] = {x,y,z};
        // Positions transformation from hybrid into hybrid spherical coordinates
        sph_transf_H2S_R(r);
        gridreal theta = r[1];
        if(negativeV == false) {
            // from the front wall
            x = Params::box_xmax_tight - V*Params::dt;
            vx = -vth*derivgaussrnd(V/vth);
        } else {
            // from the back wall
            x = Params::box_xmin_tight + V*Params::dt;
            vx = vth*derivgaussrnd(V/vth);
        }
        shortreal vy = vth*gaussrnd();
        shortreal vz = vth*gaussrnd();
        gridreal v[3] = {vx, vy, vz};
        sph_transf_S2C_V(r,v);
        vx = v[0];
        vy = v[1];
        vz = v[2];
        g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*sin(theta),popid);
    } else if(Params::sph_propagation_type == 1) { // Cartesian propagation: initial velocities are along of the one of Cartesian axis
        // Now we have flat front only from front wall
        shortreal y = Params::box_ymin_tight + uniformrnd()*Params::box_Y_tight;
        shortreal z = Params::box_zmin_tight + uniformrnd()*Params::box_Z_tight;
        shortreal x,vx,vy,vz;
        x = Params::box_xmax_tight - V*Params::dt;
        gridreal  r[3] = {x,y,z};
        // Velocities transformation, to get flat front of solar wind in shperical coordinates
        // Positions transformation from hybrid into hybrid spherical coordinates
        sph_transf_H2S_R(r);
        gridreal R     = r[0];
        gridreal theta = r[1];
        gridreal phi   = r[2];
        // There are two possibiliteis how to launch plasm particles: along x or z axis
        // 1. Particle flat front propagation along x axis
        if(Params::sph_propagation_dir == 0) {
            /// vx component of velocity (in Cartesian coordinates)
            vx = -vth*derivgaussrnd(V/vth);
            vy =  vth*gaussrnd();
            vz =  vth*gaussrnd();
            shortreal v[3] = {vx, vy, vz};
            // Positions and velocities transformation from spherical into hybrid coordinates
            vx = v[0];
            vy = v[1];
            vz = v[2];
            real t = Params::t;
            real x_max = Params::box_xmax_tight;
            real S = V*t;
            real x_launch = R*sin(theta)*cos(phi); // Launch particle position
            if (x_launch >= 0.0 && x_launch >= x_max - S) { //! (front to back)
                //flat front along x-axis
                g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*sin(theta)*sin(theta)*cos(phi),popid);
            } else {
                return;
            }
        } else if(Params::sph_propagation_dir == 1) { // 2. Particle flat front propagation along z axis
            //! vz component of velocity (in Cartesian coordinates)
            vx =  vth*gaussrnd();
            vy =  vth*gaussrnd();
            vz = -vth*derivgaussrnd(V/vth);
            shortreal v[3] = {vx, vy, vz};
            // Positions and velocities transformation from spherical into hybrid coordinates
            vx = v[0];
            vy = v[1];
            vz = v[2];
            real t = Params::t;
            real x_max = Params::box_xmax_tight;
            real S = V*t;
            real z_launch = R*cos(theta);
            if (z_launch >= 0.0 && z_launch >= x_max - S) { //! (front to back)
                //flat front along z-axis
                g.sph_addparticle(x,y,z,vx,vy,vz,macroParticleStatisticalWeight*sin(theta)*cos(theta),popid);
            } else {
                return;
            }
        }
    }
}

#endif

