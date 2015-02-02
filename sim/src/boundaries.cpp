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
#include "boundaries.h"
#include "simulation.h"
#include "random.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

//! Default constructor
ParticleBoundaryConditions::ParticleBoundaryConditions()
{
    this->ptrObstacle = &ParticleBoundaryConditions::defaultFunction;
    this->ptrSideWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrFrontWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrBackWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrExtra = &ParticleBoundaryConditions::defaultFunction;
    resetParameters();
}

#define ELSEIF_OBSTACLE(func) else if(funcName[i].compare(#func) == 0) { this->ptrObstacle = &ParticleBoundaryConditions::func; obstacleName = funcName[i]; obstacleArgs = args[i]; setArgs_ ## func(); }
#define ELSEIF_SIDEWALL(func) else if(funcName[i].compare(#func) == 0) { this->ptrSideWall = &ParticleBoundaryConditions::func; sideWallName = funcName[i]; sideWallArgs = args[i]; setArgs_ ## func(); }
#define ELSEIF_FRONTWALL(func) else if(funcName[i].compare(#func) == 0) { this->ptrFrontWall = &ParticleBoundaryConditions::func; frontWallName = funcName[i]; frontWallArgs = args[i]; setArgs_ ## func(); }
#define ELSEIF_BACKWALL(func) else if(funcName[i].compare(#func) == 0) { this->ptrBackWall = &ParticleBoundaryConditions::func; backWallName = funcName[i]; backWallArgs = args[i]; setArgs_ ## func(); }
#define ELSEIF_EXTRA(func) else if(funcName[i].compare(#func) == 0) { this->ptrExtra = &ParticleBoundaryConditions::func; extraName = funcName[i]; extraArgs = args[i]; setArgs_ ## func(); }

//! Constructor
ParticleBoundaryConditions::ParticleBoundaryConditions(vector<string> funcName,vector< vector<real> > args,unsigned int popid)
{
    if( funcName.size() != args.size() ) {
        ERRORMSG("internal error");
        doabort();
    }
    if( funcName.size() < 4 ) {
        ERRORMSG("four particle boundary condition functions required");
        doabort();
    }
    this->popid = popid;
    if(this->popid > Params::pops.size()) {
        ERRORMSG2("given population id > max(id)","ParticleBoundaryConditions");
        doabort();
    }
    this->ptrObstacle = &ParticleBoundaryConditions::defaultFunction;
    this->ptrSideWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrFrontWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrBackWall = &ParticleBoundaryConditions::defaultFunction;
    this->ptrExtra = &ParticleBoundaryConditions::defaultNullFunction;
    resetParameters();
    // Counters for the input boundary condition functions
    unsigned int obstacleN = 0;
    unsigned int sideWallN = 0;
    unsigned int frontWallN = 0;
    unsigned int backWallN = 0;
    unsigned int extraN = 0;
    // Go through the functions
    for(unsigned int i = 0; i < funcName.size(); ++i) {
        // The name should begin with a string "obstacle", "sideWall", "frontWal" or "backWall"
        if(funcName[i].length() < 9) {
            ERRORMSG2("bad particle boundary condition function name",funcName[i]);
            doabort();
        }
        // PARTICLE BOUNDARY CONDITIONS
        // Compare the first 5 characters of the name
        string cmp = funcName[i].substr(0,5);
        if(cmp.compare("obsta") == 0) {
            if( funcName[i].compare("") == 0 ) {
                ERRORMSG("empty particle boundary condition function name");
                doabort();
            }
            ELSEIF_OBSTACLE(obstacleAbsorb)
            ELSEIF_OBSTACLE(obstacleReflect)
            ELSEIF_OBSTACLE(obstacleNoObstacle)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
            ELSEIF_OBSTACLE(obstacleAbsorb_sph)
            ELSEIF_OBSTACLE(obstacleReflect_sph)
            ELSEIF_OBSTACLE(obstacleBackSideEmission_sph)
#endif
            else {
                ERRORMSG2("bad particle boundary condition function name",funcName[i]);
                doabort();
            }
            ++obstacleN;
        } else if(cmp.compare("sideW") == 0) {
            if( funcName[i].compare("") == 0 ) {
                ERRORMSG("empty particle boundary condition function name");
                doabort();
            }
            ELSEIF_SIDEWALL(sideWallAbsorb)
            ELSEIF_SIDEWALL(sideWallReflect)
            ELSEIF_SIDEWALL(sideWallPeriodic)
            ELSEIF_SIDEWALL(sideWallPeriodicYAbsorbZ)
            ELSEIF_SIDEWALL(sideWallAmbient)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
            ELSEIF_SIDEWALL(sideWallPeriodic_sph)
#endif
            else {
                ERRORMSG2("bad particle boundary condition function name",funcName[i]);
                doabort();
            }
            ++sideWallN;
        } else if(cmp.compare("front") == 0) {
            if( funcName[i].compare("") == 0 ) {
                ERRORMSG("empty particle boundary condition function name");
                doabort();
            }
            ELSEIF_FRONTWALL(frontWallAbsorb)
            ELSEIF_FRONTWALL(frontWallPeriodic)
            ELSEIF_FRONTWALL(frontWallReflect)
            ELSEIF_FRONTWALL(frontWallReflectAmbient)
            else {
                ERRORMSG2("bad particle boundary condition function name",funcName[i]);
                doabort();
            }
            ++frontWallN;
        } else if(cmp.compare("backW") == 0) {
            if( funcName[i].compare("") == 0 ) {
                ERRORMSG("empty particle boundary condition function name");
                doabort();
            }
            ELSEIF_BACKWALL(backWallAbsorb)
            ELSEIF_BACKWALL(backWallReflect)
            ELSEIF_BACKWALL(backWallPeriodic)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
            ELSEIF_BACKWALL(backWallBackSideEmission_sph)
#endif
            else {
                ERRORMSG2("bad particle boundary condition function name",funcName[i]);
                doabort();
            }
            ++backWallN;
        } else if(cmp.compare("extra") == 0) {
            if( funcName[i].compare("") == 0 ) {
                ERRORMSG("empty particle boundary condition function name");
                doabort();
            }
            ELSEIF_EXTRA(extraNeutralCollisions)
            else {
                ERRORMSG2("bad particle boundary condition function name",funcName[i]);
                doabort();
            }
            ++extraN;
        } else {
            ERRORMSG2("bad particle boundary condition function name",funcName[i]);
            doabort();
        }
    }
    // Check there is a function for every boundary - extra is optional
    if(obstacleN != 1 || sideWallN != 1 || frontWallN != 1 || backWallN != 1 || extraN>1) {
        ERRORMSG("need boundary condition functions for obstacle, sidewall, frontwall and backwall");
        doabort();
    }
}

//! Destructor
ParticleBoundaryConditions::~ParticleBoundaryConditions() { }

//! Check boundary conditions for a particle
bool ParticleBoundaryConditions::checkBoundaries(TLinkedParticle& p,fastreal rAverage[])
{
    bool rAverageFlag = true; //if velocity or position of particle is changed the raveflag needs to be false
    bool keepParticle = true;
    // Check obstacle
    keepParticle = (this->*ptrObstacle)(p,rAverageFlag);
    if(keepParticle == false) {
#ifndef NO_DIAGNOSTICS
        // Count impacting particles
        Params::diag.pCounter[p.popid]->increaseImpactCounters(p);
#endif
        return false;
    }
    // Check back wall
    keepParticle = (this->*ptrBackWall)(p,rAverageFlag);
    if(keepParticle == false) {
#ifndef NO_DIAGNOSTICS
        // Count escaping particles
        Params::diag.pCounter[p.popid]->increaseEscapeCountersBackWall(p);
#endif
        return false;
    }
    // Check sidewalls
    keepParticle = (this->*ptrSideWall)(p,rAverageFlag);
    if(keepParticle == false) {
#ifndef NO_DIAGNOSTICS
        // Count escaping particles
        Params::diag.pCounter[p.popid]->increaseEscapeCountersSideWall(p);
#endif
        return false;
    }
    // Check front wall
    keepParticle = (this->*ptrFrontWall)(p,rAverageFlag);
    if(keepParticle == false) {
#ifndef NO_DIAGNOSTICS
        // Count escaping particles
        Params::diag.pCounter[p.popid]->increaseEscapeCountersFrontWall(p);
#endif
        return false;
    }

    // Check extra - stopped particles are counted in ImpactCounters!
    keepParticle = (this->*ptrExtra)(p,rAverageFlag);
    if(keepParticle == false) {
#ifndef NO_DIAGNOSTICS
        // Count escaping particles
        Params::diag.pCounter[p.popid]->increaseImpactCounters(p);
#endif
        return false;
    }

    //check if inside the box (tight boundaries
    bool outside = false;
    if(Params::insideBoxTightFrontWall(&p) == false) {
#ifndef NO_DIAGNOSTICS
        Params::diag.pCounter[p.popid]->increaseEscapeCountersFrontWall(p);
#endif
        outside = true;
    } else if(Params::insideBoxTightBackWall(&p) == false) {
#ifndef NO_DIAGNOSTICS
        Params::diag.pCounter[p.popid]->increaseEscapeCountersBackWall(p);
#endif
        outside = true;
    } else if(Params::insideBoxTightSideWall(&p) == false) {
#ifndef NO_DIAGNOSTICS
        Params::diag.pCounter[p.popid]->increaseEscapeCountersSideWall(p);
#endif
        outside = true;
    }
    if(outside == true) {
        errorlog
                << "Boundaries returned particle that is outside the domain:  pop = " << Params::pops[p.popid]->getIdStr()
                << ", r = [" << (int(p.x/Params::R_P*100))/100.0 << ", "
                << (int(p.y/Params::R_P*100))/100.0 << ", "
                << (int(p.z/Params::R_P*100))/100.0
                << "], v = " << sqrt(sqr(p.vx)+sqr(p.vy)+sqr(p.vz)) << " m/s"
                << ", rave = " << rAverageFlag << "\n";
        return false;
    }
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    if (rAverageFlag) {
        //half step back
#ifndef USE_PARTICLE_SUBCYCLING
        rAverage[0] = p.x - 0.5*p.vx*Params::dt;
        rAverage[1] = p.y - 0.5*p.vy*Params::dt;
        rAverage[2] = p.z - 0.5*p.vz*Params::dt;
#else
        rAverage[0] = p.x - 0.5*p.vx*Params::dt_psub[p.dtlevel];
        rAverage[1] = p.y - 0.5*p.vy*Params::dt_psub[p.dtlevel];
        rAverage[2] = p.z - 0.5*p.vz*Params::dt_psub[p.dtlevel];
#endif
        if (rAverage[0] < Params::box_xmin_tight || rAverage[0] > Params::box_xmax_tight ||
            rAverage[1] < Params::box_ymin_tight || rAverage[1] > Params::box_ymax_tight ||
            rAverage[2] < Params::box_zmin_tight || rAverage[2] > Params::box_zmax_tight) {
            // outside
            errorlog << "Boundaries: rAverage outside the domain (may be raveflag issue) \n"
                     << " pop = " << Params::pops[p.popid]->getIdStr()
                     << ", rA = [" << (int(rAverage[0]/Params::R_P*100))/100.0 << ", "
                     << (int(rAverage[1]/Params::R_P*100))/100.0 << ", "
                     << (int(rAverage[2]/Params::R_P*100))/100.0
                     << "], v = [" << p.vx << ", " << p.vy << ", " << p.vz << "] m/s" << "\n";
            rAverage[0] = p.x;
            rAverage[1] = p.y;
            rAverage[2] = p.z;
        }
    } else {
        rAverage[0] = p.x;
        rAverage[1] = p.y;
        rAverage[2] = p.z;
    }
#else
    rAverage[0] = p.x;
    rAverage[1] = p.y;
    rAverage[2] = p.z;
#endif
    return true;
}

//! String summary of the object
string ParticleBoundaryConditions::toString(string delim)
{
    stringstream ss;
    ss << " " << obstacleName << " ";
    for(unsigned int ii = 0; ii < obstacleArgs.size(); ++ii) {
        ss << obstacleArgs[ii] << " ";
    }
    ss << delim;
    ss << " " << backWallName << " ";
    for(unsigned int ii = 0; ii < backWallArgs.size(); ++ii) {
        ss << backWallArgs[ii] << " ";
    }
    ss << delim;
    ss << " " << sideWallName << " ";
    for(unsigned int ii = 0; ii < sideWallArgs.size(); ++ii) {
        ss << sideWallArgs[ii] << " ";
    }
    ss << delim;
    ss << " " << frontWallName << " ";
    for(unsigned int ii = 0; ii < frontWallArgs.size(); ++ii) {
        ss << frontWallArgs[ii] << " ";
    }
    return ss.str();
}

//! Default function, which aborts the program if called
bool ParticleBoundaryConditions::defaultFunction(TLinkedParticle& p,bool& raveflag)
{
    ERRORMSG("function pointer not set");
    doabort();
    return true;
}

//! Default function for extra type
bool ParticleBoundaryConditions::defaultNullFunction(TLinkedParticle& p,bool& raveflag)
{
    return true;
}

void ParticleBoundaryConditions::resetParameters()
{
    // PARTICLE BOUNDARY CONDITION PARAMETERS
    reflectPopID = -1;
    R = 0.0;
    V = 0.0;
    vth = 0.0;
}

// PARTICLE BOUNDARY CONDITIONS

//! Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleAbsorb()
{
    if(obstacleArgs.size() != 1) {
        ERRORMSG2("function takes one argument",obstacleName);
        doabort();
    }
    R = obstacleArgs[0];
    if(R >= 0) {
        R2 = sqr(R);
    } else {
        R2 = -1.0*sqr(R);
    }
}

//! Absorbing obstacle
bool ParticleBoundaryConditions::obstacleAbsorb(TLinkedParticle& p,bool& raveflag)
{
    if (sqr(p.x) + sqr(p.y) + sqr(p.z) < R2) {
        return false;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleReflect()
{
    if(obstacleArgs.size() != 1) {
        ERRORMSG2("function takes one argument",obstacleName);
        doabort();
    }
    R = obstacleArgs[0];
    R2 = sqr(R);
}

//! Reflecting obstacle
bool ParticleBoundaryConditions::obstacleReflect(TLinkedParticle& p,bool& raveflag)
{
    if(sqr(p.x) + sqr(p.y) + sqr(p.z) < R2) {
        raveflag = false;
        p.vx = -p.vx;
        p.vy = -p.vy;
        p.vz = -p.vz;
        p.x += p.vx*Params::dt;
        p.y += p.vy*Params::dt;
        p.z += p.vz*Params::dt;
        // particle not moved during this time step - only velocity turned
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleNoObstacle()
{
    if(obstacleArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",obstacleName);
        doabort();
    }
}

//! No obstacle
bool ParticleBoundaryConditions::obstacleNoObstacle(TLinkedParticle& p,bool& raveflag)
{
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallAbsorb()
{
    if(sideWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",sideWallName);
        doabort();
    }
}

//! Absorbing side walls
bool ParticleBoundaryConditions::sideWallAbsorb(TLinkedParticle& p,bool& raveflag)
{
    if (p.y < Params::box_ymin_tight || p.y > Params::box_ymax_tight) {
        return false;
    }
    if (p.z < Params::box_zmin_tight || p.z > Params::box_zmax_tight) {
        return false;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallReflect()
{
    if(sideWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",sideWallName);
        doabort();
    }
}

//! Reflecting side walls
bool ParticleBoundaryConditions::sideWallReflect(TLinkedParticle& p,bool& raveflag)
{
    if (p.y < Params::box_ymin_tight) {
        raveflag = false;
        p.vy = -p.vy;
        p.y += p.vy*Params::dt;
    } else if (p.y > Params::box_ymax_tight) {
        raveflag = false;
        p.vy = -p.vy;
        p.y += p.vy*Params::dt;
    }
    if (p.z < Params::box_zmin_tight) {
        raveflag = false;
        p.vz = -p.vz;
        p.z += p.vz*Params::dt;
    } else if (p.z > Params::box_zmax_tight) {
        raveflag = false;
        p.vz = -p.vz;
        p.z += p.vz*Params::dt;
        // particle not moved during this time step - only velocity turned
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallPeriodic()
{
    if(sideWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",sideWallName);
        doabort();
    }
}

//! Periodic side walls
bool ParticleBoundaryConditions::sideWallPeriodic(TLinkedParticle& p,bool& raveflag)
{
    if (p.y < Params::box_ymin_tight) {
        raveflag = false;
        p.y += Params::box_Y_tight;
    } else if (p.y > Params::box_ymax_tight) {
        raveflag = false;
        p.y -= Params::box_Y_tight;
    }
    if (p.z < Params::box_zmin_tight) {
        raveflag = false;
        p.z += Params::box_Z_tight;
    } else if (p.z > Params::box_zmax_tight) {
        raveflag = false;
        p.z -= Params::box_Z_tight;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallPeriodicYAbsorbZ()
{
    if(sideWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",sideWallName);
        doabort();
    }
}

//! Periodic side walls
bool ParticleBoundaryConditions::sideWallPeriodicYAbsorbZ(TLinkedParticle& p,bool& raveflag)
{
    if (p.z < Params::box_zmin_tight || p.z > Params::box_zmax_tight) {
        return false;
    }

    if (p.y < Params::box_ymin_tight) {
        raveflag = false;
        p.y += Params::box_Y_tight;
    } else if (p.y > Params::box_ymax_tight) {
        raveflag = false;
        p.y -= Params::box_Y_tight;
    }

    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallAmbient()
{
    if(sideWallArgs.size() != 1) {
        ERRORMSG2("function takes one argument",sideWallName);
        doabort();
    }
    V = sideWallArgs[0];
    vth = Params::pops[this->popid]->getThermalSpeed();
}

//! Ambient side walls
bool ParticleBoundaryConditions::sideWallAmbient(TLinkedParticle& p,bool& raveflag)
{
    //Y walls
    if (p.y < Params::box_ymin_tight || p.y > Params::box_ymax_tight) {
        raveflag = false;
        p.x = Params::box_xmin_tight + uniformrnd()*(Params::box_X_tight);
        p.z = Params::box_zmin_tight + uniformrnd()*(Params::box_Z_tight);
        p.vx = -V + vth*gaussrnd();
        p.vz = vth*gaussrnd();
        if (uniformrnd() < 0.50) {
            p.y = Params::box_ymin_tight + 0.1*V*Params::dt;
            p.vy = vth*derivgaussrnd(0);
        } else {
            p.y = Params::box_ymax_tight - 0.1*V*Params::dt;
            p.vy = - vth*derivgaussrnd(0);
        }
    } else {
        if (p.z < Params::box_zmin_tight || p.z > Params::box_zmax_tight) {
            //Z walls
            raveflag = false;
            p.x = Params::box_xmin_tight + uniformrnd()*(Params::box_X_tight);
            p.y = Params::box_ymin_tight + uniformrnd()*(Params::box_Y_tight);
            p.vx = -V + vth*gaussrnd();
            p.vy = vth*gaussrnd();

            if (uniformrnd() < 0.50) {
                p.z = Params::box_zmin_tight + 0.1*V*Params::dt;
                p.vz = vth*derivgaussrnd(0);
            } else {
                p.z = Params::box_zmax_tight - 0.1*V*Params::dt;
                p.vz = -vth*derivgaussrnd(0);
            }
        }
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_frontWallAbsorb()
{
    if(frontWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",frontWallName);
        doabort();
    }
}

//! Absorbing front wall
bool ParticleBoundaryConditions::frontWallAbsorb(TLinkedParticle& p,bool& raveflag)
{
    if(p.x > Params::box_xmax_tight) {
        return false;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_frontWallPeriodic()
{
    if(frontWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",frontWallName);
        doabort();
    }
}

//! Periodic front->back wall
bool ParticleBoundaryConditions::frontWallPeriodic(TLinkedParticle& p,bool& raveflag)
{
    if(p.x > Params::box_xmax_tight) {
        raveflag = false;
        p.x -= Params::box_X_tight;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_frontWallReflect()
{
    if(frontWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",frontWallName);
        doabort();
    }
}

//! Periodic front->back wall
bool ParticleBoundaryConditions::frontWallReflect(TLinkedParticle& p,bool& raveflag)
{
    if(p.x > Params::box_xmax_tight) {
        raveflag = false;
        p.vx = -p.vx;
        p.x += p.vx*Params::dt;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_frontWallReflectAmbient()
{
    if(frontWallArgs.size() != 1) {
        ERRORMSG2("function takes one argument",frontWallName);
        doabort();
    }
    reflectPopID = static_cast<int>(frontWallArgs[0]);
    if(reflectPopID < 0) {
        ERRORMSG2("bad popid for reflected population",frontWallName);
        doabort();
    }
}

//! Reflect a particle by deleting the old particle and adding a new one using given popid
bool ParticleBoundaryConditions::frontWallReflectAmbient(TLinkedParticle& p,bool& raveflag)
{
    if(p.x > Params::box_xmax_tight) {
        shortreal x = p.x - p.vx*Params::dt;
        Params::pops[reflectPopID]->addParticle(x,p.y,p.z,p.w);
        return false;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_backWallAbsorb()
{
    if(backWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",backWallName);
        doabort();
    }
}

//! Absorbing back wall
bool ParticleBoundaryConditions::backWallAbsorb(TLinkedParticle& p,bool& raveflag)
{
    if (p.x < Params::box_xmin_tight) {
        return false;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_backWallReflect()
{
    if(backWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",backWallName);
        doabort();
    }
}

//! Reflecting back wall
bool ParticleBoundaryConditions::backWallReflect(TLinkedParticle& p,bool& raveflag)
{
    if(p.x < Params::box_xmin_tight) {
        raveflag = false;
        p.vx = -p.vx;
        p.x += p.vx*Params::dt;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_backWallPeriodic()
{
    if(backWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",backWallName);
        doabort();
    }
}

//! Periodic back->front wall
bool ParticleBoundaryConditions::backWallPeriodic(TLinkedParticle& p,bool& raveflag)
{
    if(p.x < Params::box_xmin_tight) {
        raveflag = false;
        p.x += Params::box_X_tight;
    }
    return true;
}

//! Set function arguments
void ParticleBoundaryConditions::setArgs_extraNeutralCollisions()
{
    if(extraArgs.size() != 2) {
        ERRORMSG2("function takes two arguments",extraName);
        doabort();
    }
    R = extraArgs[0]; //lower limit
    R2 = sqr(R);
    Rcoll2 = sqr(extraArgs[1]); // upper radius limit
}

//! Neutral collisions with absortion (stopping cross section)
bool ParticleBoundaryConditions::extraNeutralCollisions(TLinkedParticle& p,bool& raveflag)
{
    gridreal r2=sqr(p.x) + sqr(p.y) + sqr(p.z);
    if(r2 < R2) {
        return false;
    }
    if(r2 > Rcoll2) {
        return true;    //above collision region
    }
    //Calculate energy exchange and new velocity for the particle
    gridreal r[3]= {p.x,p.y,p.z};
    gridreal density = Params::pops[0]->getNeutralDensity(r); //number density
    //'NeutralDensity' (population id 0 exospheric - no need to create MPs - only densityprofile};
    if (density<=0) {
        return true;
    }
    gridreal v2=sqr(p.vx) + sqr(p.vy) +sqr(p.vz);
    gridreal newv, reduction, v=sqrt(v2);
    gridreal nccs=NCcross(v2,Params::pops[p.popid]->getIdStr());
    gridreal deltaE = nccs*v*Params::dt*density; //cross section * path length * neutral density
    if (deltaE<=0) {
        return true;
    }
    raveflag = false; //better would be to calculate rAverage before changes to v
    gridreal E = 0.5*Params::pops[p.popid]->m*v2;
    //the delta_energy/E fraction should be always small - stopping should not happen!!
    // we 'subcycle' by doing the process repeatedly (n times) in the same location using time step dt/n
    int iterate=1;
    const fastreal fr_limit=0.15; // 15% limit for energy loss
    if (fr_limit * E <= deltaE) {
        //const real dE=E*p.w;
#ifndef NO_DIAGNOSTICS
        Params::diag.pCounter[p.popid]->cutRateV += 1.0;//count subcycling rate
#endif
        const fastreal fr = deltaE/E;
        if (fr > 1.2) { //stop the particle
            p.vx=0;
            p.vy=0;
            p.vz=0;
            return true;
#ifndef NO_DIAGNOSTICS
            Params::diag.pCounter[p.popid]->electronImpactIonizationRate+=E*p.w; //dE=deltaE*p.w
#endif
        }
        iterate = min( int (26 * sqr(fr) + 11 * fr + 0.24), 25); // 'optimal' function
        //const fastreal dt_new = Params::dt/iterate;
        deltaE=deltaE/iterate; //only dt was changed
    }
    E-=deltaE;
    for (int i=0; i<iterate; i++) {
#ifndef NO_DIAGNOSTICS
        Params::diag.pCounter[p.popid]->electronImpactIonizationRate+=deltaE*p.w; //dE=deltaE*p.w //add momentum counter?
        //g.increaseQ2H(r,deltaE); //store in field quantity for hc (U1 in avehc)
#endif
        newv=sqrt(E*2/Params::pops[p.popid]->m);
        reduction=newv/v;
        p.vx *= reduction;
        p.vy *= reduction;
        p.vz *= reduction;
        if (iterate>1 && i<iterate-1) {
            v = newv;
            nccs=NCcross(sqr(v),Params::pops[p.popid]->getIdStr());
            deltaE = nccs*v*Params::dt/iterate*density;
            //if (fr_limit * E * 1.5 < deltaE){
            //errorlog << "NC subcycling gave dE/E fraction: " << deltaE/E << "\n";
            if (E < deltaE) {
                p.vx=0;
                p.vy=0;
                p.vz=0;
                return true;
            }
            //stopped but not removed (energy not counted)
            //}
            E -= deltaE;
        }
    }
    return true;
}

//! Cross section for neutral collisions
real ParticleBoundaryConditions::NCcross(real v2, string popstr)
{
    int npop=8;
    string NCpops[8]= {"flowH+","exoH+","flowH2+","exoH2+","flowO+","flowOH+","exoCH4+","ionoN2+"}; //8
    //check that the NC is defined for the particle (population)
    int NCcheck=0;
    int NCpop;
    for (int ni=0; ni<npop; ni++) {
        if (popstr.compare(NCpops[ni])==0) { //found
            NCcheck=1;
            NCpop=ni;
            break;
        }
    }
    if (NCcheck==0) {
        WARNINGMSG2("Neutral collisions Boundary condition ''extraNeutralCollisions'' not defined for population: ",popstr);
        return 0;
    }
    real cs; //O+ and H+ stopping cross sections (combined sn (nucleus) and se (electron)) in N gas
    // from a fit to SRIM parameters (O.J. Tucker & Robert Johnson)
    //cross section is sum of H+ and O+ crosssections (using the corresponding velocity i.e. E'=m_p*v2 or m_O*v2)
    int nH=0;
    int nO=0;
    if (NCpop==0 || NCpop==1) {
        nH=1;   //H+
    }
    if (NCpop==2 || NCpop==3) {
        nH=2;   //H2+
    }
    if (NCpop==4 || NCpop==5) {
        nO=1;   //O+  (OH+ is what INUM moments calls the watergroup ions, at Titan O+)
    }
    if (NCpop==6) {
        nH=4;    //CH4+, as C+ and 4H+
        nO=1;
    }
    if (NCpop==7) {
        nO=2;   //N2+, as 2N+
    }
    real csH=0;
    real csO=0;
    real E_eV,E1,E2,E3;
    if (nH>0) { // up to 8th potence of E needed.
        E_eV = 0.5*v2*Params::m_p/Params::e;
        if (E_eV>=3.0 || E_eV<=3.0e6) {
            E1=log(E_eV);
            E2=sqr(E1);
            E3=E2*E1;
            csH=exp(4.40261e-5*E3*E3-0.00207302*E3*E2+0.0373833*E3*E1-0.330106*E3+1.5043*E2-2.98386*E1-33.0139);
            if (csH<0) {
                csH=0;
            }
        }
    }
    if (nO>0) {
        E_eV = 0.5*v2*Params::m_O/Params::e;
        if (E_eV>=1.1 || E_eV<=6e4) {
            E1=log(E_eV);
            csO=exp(-0.0349419*E1*E1+0.709552*E1-34.4666);
        }
        if (csO<0) {
            csO=0;
        }
    }
    cs=(nH*csH+nO*csO)*Params::e*1e-4; //transform eV*cm2 to J*m2
    return 2*cs; //for N2 instead of N gas
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleAbsorb_sph()
{
    if(obstacleArgs.size() != 1) {
        ERRORMSG2("function takes one argument",obstacleName);
        doabort();
    }
    R = obstacleArgs[0];
    if(R >= 0) {
        R2 = sqr(R);
    } else {
        R2 = -1.0*sqr(R);
    }
}

//! (SPHERICAL) Spherical version of "obstacleAbsorb"
bool ParticleBoundaryConditions::obstacleAbsorb_sph(TLinkedParticle& p,bool& raveflag)
{
    if ((p.x) < R) {
        return false;
    }
    return true;
}

//! (SPHERICAL) Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleReflect_sph()
{
    if(obstacleArgs.size() != 1) {
        ERRORMSG2("function takes one argument",obstacleName);
        doabort();
    }
    R = obstacleArgs[0];
    R2 = sqr(R);
}

//! (SPHERICAL) Spherical version of "obstacleReflect"
bool ParticleBoundaryConditions::obstacleReflect_sph(TLinkedParticle& p,bool& raveflag)
{
    if (p.x < R) {
        gridreal r_old[3] = {p.x, p.y, p.z};
        gridreal v_old[3] = {p.vx, p.vy, p.vz};
        gridreal v_new[3] = {0, 0, 0};
        gridreal r_new[3] = {0, 0, 0};
        sph_transf_H2S_R(r_old);
        sph_transf_C2S_V(r_old, v_old);
        raveflag = false;
        v_new[0] = -v_old[0];
        v_new[1] =  v_old[1];
        v_new[2] =  v_old[2];
        sph_transf_S2H_V(v_new);
        p.x += v_new[0]*Params::dt;
        p.y += v_new[1]*Params::dt;
        p.z += v_new[2]*Params::dt;
        sph_transf_H2S_V(v_new);
        sph_transf_S2C_V(r_old, v_new);
        p.vx = v_new[0];
        p.vy = v_new[1];
        p.vz = v_new[2];
    }
    return true;
}

//! (SPHERICAL) Set function arguments
void ParticleBoundaryConditions::setArgs_obstacleBackSideEmission_sph()
{
    if(obstacleArgs.size() != 1) {
        ERRORMSG2("function takes one argument",obstacleName);
        doabort();
    }
    R = obstacleArgs[0];
}

//! (SPHERICAL) Particle emission from back side of obstacle
bool ParticleBoundaryConditions::obstacleBackSideEmission_sph(TLinkedParticle& p,bool& raveflag)
{
    gridreal V  = abs(Params::V);
    gridreal t  = Params::t;
    gridreal dt = Params::dt;
    gridreal R_max = Params::box_xmax;
    gridreal r_new[3];
    gridreal v_new[3];
    if (p.x < R) {
        raveflag = false;
        gridreal r_old[3] = {p.x, p.y, p.z};
        gridreal v_old[3] = {p.vx, p.vy, p.vz};
        sph_transf_H2S_R(r_old);
        sph_transf_S2C_R(r_old);
        // This algorithm is valid for particle propagation along z - axis. if somebody needs to push particle along another Cartesian axis, has to change index in r_old (0 - x, 1 - y, 2 - z)
        // Checking the possitive value of z to avoid particle translation from back to front side of obstacle
        if (r_old[2] > 0.0 && r_old[2] < V*t - R_max) {
            v_new[0] =  v_old[0];
            v_new[1] =  v_old[1];
            v_new[2] =  v_old[2];
            r_new[0] =  r_old[0] + v_new[0]*dt;
            r_new[1] =  r_old[1] + v_new[1]*dt;
            r_new[2] = -r_old[2] + v_new[2]*dt;
            sph_transf_C2S_r(r_new);
            sph_transf_S2H_R(r_new);
            p.x  = r_new[0];
            p.y  = r_new[1];
            p.z  = r_new[2];
            p.vx = v_new[0];
            p.vy = v_new[1];
            p.vz = v_new[2];
        }
    }
    return true;
}

//! (SPHERICAL) Set function arguments
void ParticleBoundaryConditions::setArgs_sideWallPeriodic_sph()
{
    if(sideWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",sideWallName);
        doabort();
    }
}

//! (SPHERICAL) Periodic side walls
bool ParticleBoundaryConditions::sideWallPeriodic_sph(TLinkedParticle& p,bool& raveflag)
{
// This is very simple method that is not optimal but the probabiliy that particle cross the pole line is very low
//particles boundary conditions for theta
    if (p.y < Params::box_ymin_tight) {
        raveflag = false;
        if (p.z <= 0.0) {
            p.y  = Params::box_ymin_tight;
            p.z += Params::box_Z_tight/2.0;
        } else {
            p.y  = Params::box_ymin_tight;
            p.z -= Params::box_Z_tight/2.0;
        }
    } else if (p.y > Params::box_ymax_tight) {
        raveflag = false;
        if (p.z <= 0.0) {
            p.y  = Params::box_ymax_tight;
            p.z += Params::box_Z_tight/2.0;
        } else {
            p.y  = Params::box_ymax_tight;
            p.z -= Params::box_Z_tight/2.0;
        }
    }
//particle boundary conditions for phi
    if (p.z < Params::box_zmin_tight) {
        raveflag = false;
        p.z += Params::box_Z_tight;
    } else if (p.z > Params::box_zmax_tight) {
        raveflag = false;
        p.z -= Params::box_Z_tight;
    }
    return true;
}

//! (SPHERICAL) Set function arguments
void ParticleBoundaryConditions::setArgs_backWallBackSideEmission_sph()
{
    if(backWallArgs.size() != 0) {
        ERRORMSG2("function takes no arguments",backWallName);
        doabort();
    }
}

//! (SPHERICAL) Particle emission from the back side of obstacle.
bool ParticleBoundaryConditions::backWallBackSideEmission_sph(TLinkedParticle& p,bool& raveflag)
{
    gridreal V  = abs(Params::V);
    gridreal t  = Params::t;
    gridreal dt = Params::dt;
    gridreal R_max = Params::box_xmax;
    gridreal r_new[3];
    gridreal v_new[3];
    if (p.x < Params::box_xmin + 10*V*dt) {
        raveflag = false;
        gridreal r_old[3] = {p.x, p.y, p.z};
        gridreal v_old[3] = {p.vx, p.vy, p.vz};
        sph_transf_H2S_R(r_old);
        sph_transf_H2S_V(v_old);
        sph_transf_S2C_V(r_old, v_old);
        sph_transf_S2C_R(r_old);
        if (r_old[2] > 0.0 && r_old[2] < V*t - R_max) {
            v_new[0] =  v_old[0];
            v_new[1] =  v_old[1];
            v_new[2] =  v_old[2];
            r_new[0] =  r_old[0] + v_new[0]*dt;
            r_new[1] =  r_old[1] + v_new[1]*dt;
            r_new[2] = -r_old[2] + v_new[2]*dt;
            sph_transf_C2S_r(r_new);
            sph_transf_C2S_V(r_new, v_new);
            sph_transf_S2H_R(r_new);
            sph_transf_S2H_V(v_new);
            p.x  = r_new[0];
            p.y  = r_new[1];
            p.z  = r_new[2];
            p.vx = v_new[0];
            p.vy = v_new[1];
            p.vz = v_new[2];
        }
    }
    return true;
}

#endif

