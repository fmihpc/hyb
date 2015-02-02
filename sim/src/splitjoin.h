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

#ifndef SPLITJOIN_H
#define SPLITJOIN_H

#include <iostream>
#include <string>
#include <vector>
#include "definitions.h"
#include "particle.h"

//! Macro particle splitting
class Split
{
public:
    Split();
    Split(std::string funcName,std::vector<real> args);
    ~Split();
    int doSplitting(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist);
private:
    std::string name;
    std::vector<real> args;
    int (Split::*ptr)(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist);
    int defaultFunction(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist);
    bool split12(TLinkedParticle& p, gridreal split_distance);
    // SPLIT METHODS
    int splitDefault(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist);
    void setArgs_splitDefault();
    int splitOriginal(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist);
    void setArgs_splitOriginal();
    // SPLIT PARAMETERS
    void resetParameters();
    gridreal distanceFactor;
    bool HALFCELL_SPLIT_DISTANCE;
};

//! Macro particle joining
class Join
{
public:
    Join();
    Join(std::string funcName,std::vector<real> args);
    ~Join();
    int doJoining(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int fast);
private:
    std::string name;
    std::vector<real> args;
    int (Join::*ptr)(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int);
    int defaultFunction(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int);
    void join32(const TLinkedParticle& P1, const TLinkedParticle& P2, const TLinkedParticle& P3, TLinkedParticle& A, TLinkedParticle& B);
    // JOIN METHODS
    int joinDefault(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int fast);
    void setArgs_joinDefault();
    int joinOriginal(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int);
    void setArgs_joinOriginal();
    // JOIN PARAMETERS
    //none
    void resetParameters();
};

void initializeSplitJoin();

#endif

