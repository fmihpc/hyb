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

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "definitions.h"
#include "particle.h"

//! Particle boundary conditions
class ParticleBoundaryConditions
{
public:
    ParticleBoundaryConditions();
    ParticleBoundaryConditions(std::vector<std::string> funcName,std::vector< std::vector<real> > args,unsigned int popid);
    ~ParticleBoundaryConditions();
    bool checkBoundaries(TLinkedParticle& p,fastreal rAverage[]);
    std::string toString(std::string delim=std::string("\n"));
private:
    unsigned int popid; //!< ID of the population for these boundary conditions
    std::string obstacleName; //!< Name of the inner boundary condition
    std::string sideWallName; //!< Name of the side wall boundary condition
    std::string frontWallName; //!< Name of the front wall boundary condition
    std::string backWallName; //!< Name of the back wall boundary condition
    std::string extraName; //!< Name of the neutral collision boundary condition
    std::vector<real> obstacleArgs; //!< Arguments of the inner boundary condition
    std::vector<real> sideWallArgs; //!< Arguments of the side wall boundary condition
    std::vector<real> frontWallArgs; //!< Arguments of the front wall boundary condition
    std::vector<real> backWallArgs; //!< Arguments of the back wall boundary condition
    std::vector<real> extraArgs; //!< Arguments of the neutral collision boundary condition
    bool (ParticleBoundaryConditions::*ptrObstacle)(TLinkedParticle& p,bool& raveflag); //!< Function pointer of the inner boundary condition
    bool (ParticleBoundaryConditions::*ptrSideWall)(TLinkedParticle& p,bool& raveflag); //!< Function pointer of the side wall boundary condition
    bool (ParticleBoundaryConditions::*ptrFrontWall)(TLinkedParticle& p,bool& raveflag); //!< Function pointer of the front wall boundary condition
    bool (ParticleBoundaryConditions::*ptrBackWall)(TLinkedParticle& p,bool& raveflag); //!< Function pointer of the back wall boundary condition
    bool (ParticleBoundaryConditions::*ptrExtra)(TLinkedParticle& p,bool& raveflag); //!< Function pointer of the neutral collision boundary condition
    bool defaultFunction(TLinkedParticle& p,bool& raveflag);
    bool defaultNullFunction(TLinkedParticle& p,bool& raveflag);
    // PARTICLE BOUNDARY CONDITIONS
    bool obstacleAbsorb(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleAbsorb();
    bool obstacleReflect(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleReflect();
    bool obstacleNoObstacle(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleNoObstacle();
    bool sideWallAbsorb(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallAbsorb();
    bool sideWallReflect(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallReflect();
    bool sideWallPeriodic(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallPeriodic();
    bool sideWallPeriodicYAbsorbZ(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallPeriodicYAbsorbZ();
    bool sideWallAmbient(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallAmbient();
    bool frontWallAbsorb(TLinkedParticle& p,bool& raveflag);
    void setArgs_frontWallAbsorb();
    bool frontWallPeriodic(TLinkedParticle& p,bool& raveflag);
    void setArgs_frontWallPeriodic();
    bool frontWallReflect(TLinkedParticle& p,bool& raveflag);
    void setArgs_frontWallReflect();
    bool frontWallReflectAmbient(TLinkedParticle& p,bool& raveflag);
    void setArgs_frontWallReflectAmbient();
    bool backWallAbsorb(TLinkedParticle& p,bool& raveflag);
    void setArgs_backWallAbsorb();
    bool backWallReflect(TLinkedParticle& p,bool& raveflag);
    void setArgs_backWallReflect();
    bool backWallPeriodic(TLinkedParticle& p,bool& raveflag);
    void setArgs_backWallPeriodic();
    bool extraNeutralCollisions(TLinkedParticle& p,bool& raveflag);
    void setArgs_extraNeutralCollisions();
    // PARTICLE BOUNDARY CONDITION PARAMETERS
    void resetParameters();
    int reflectPopID;
    real R,R2,V,vth,Rcoll2;
    real NCcross(real v2, std::string popIdStr);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    bool obstacleAbsorb_sph(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleAbsorb_sph();
    bool obstacleReflect_sph(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleReflect_sph();
    bool obstacleBackSideEmission_sph(TLinkedParticle& p,bool& raveflag);
    void setArgs_obstacleBackSideEmission_sph();
    bool sideWallPeriodic_sph(TLinkedParticle& p,bool& raveflag);
    void setArgs_sideWallPeriodic_sph();
    bool backWallBackSideEmission_sph(TLinkedParticle& p,bool& raveflag);
    void setArgs_backWallBackSideEmission_sph();
#endif
};

#endif

