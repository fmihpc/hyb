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

#ifndef POPULATION_H
#define POPULATION_H

#include "atmosphere.h"
#include "boundaries.h"
#include "logger.h"

int getPopId(std::string popIdStr);

struct stringArg {
    std::string value;
    bool given;
};
struct realArg {
    real value;
    bool given;
};
struct realArg2 {
    real value[2];
    bool given;
};
struct fastrealArg {
    fastreal value;
    bool given;
};
struct boolArg {
    bool value;
    bool given;
};
struct functionArg {
    std::string name;
    std::vector<real> funcArgs;
    bool given;
};
struct functionArg2 {
    std::vector<std::string> name;
    std::vector< std::vector<real> > funcArgs;
    bool given;
};

//! All particle population arguments (including population type specific ones)
struct PopulationArgs {
    // TEMPLATE POPULATION VARIABLES
    // general population parameters
    stringArg idStr;
    stringArg hcFilePrefix;
    functionArg2 boundaryFUNC;
    realArg m;
    realArg q;
    realArg T;
    realArg vth;
    realArg macroParticlesPerDt;
    boolArg propagateV;
    boolArg accumulate;
    boolArg split;
    boolArg join;
    boolArg logParams;
    // population type specific parameters
    realArg n;
    realArg V;
    realArg backWallWeight;
    realArg R;
    realArg totalRate;
    functionArg2 distFunc;
#ifdef USE_PARTICLE_SUBCYCLING
    realArg subcycleSteps;
#endif
    PopulationArgs();
    void clearArgs();
};

//! General particle population
class Population
{
public:
    Population();
    Population(PopulationArgs args);
    const real m;
    const real q;
    virtual ~Population();
    virtual void initialize();
    virtual void createParticles();
    virtual void addParticle(shortreal x,shortreal y,shortreal z,real w);
    virtual real getNeutralDensity(const gridreal[]);
    virtual void updateArgs();
    virtual void writeExtraHcFile();
    virtual std::string configDump();
    virtual std::string toString();
    void clearArgs() {
        args.clearArgs();
    };
    void saveArgs(PopulationArgs args);
    int getPopId() {
        return popid;
    }
    static void getHcFileConfigs(std::vector<std::string>& filePrefix, std::vector< std::vector<int> >& popId);
    std::string getIdStr() {
        return idStr;
    }
    bool checkBoundaries(TLinkedParticle& p,fastreal rAverage[]);
    real getThermalSpeed() {
        return vth;
    }
    bool getPropagateV() {
        return propagateV;
    }
    bool getAccumulate() {
        return accumulate;
    }
    bool getSplit() {
        return split;
    }
    bool getJoin() {
        return join;
    }
#ifdef USE_PARTICLE_SUBCYCLING
    real subcycleSteps;
#endif
private:
    static unsigned int idCnt;
    static std::vector<std::string> idStrTbl;
    static std::vector<std::string> hcFilePrefixTbl;
    ParticleBoundaryConditions boundaries;
    void checkMassAndCharge();
    void setIdStr(const std::string str);
    void setHcFilePrefix(const std::string prefix);
protected:
    unsigned int popid;
    PopulationArgs args;
    std::string idStr;
    std::string hcFilePrefix;
    std::vector<std::string> boundaryFuncNames;
    std::vector< std::vector<real> > boundaryFuncArgs;
    real T;
    real vth;
    real macroParticleStatisticalWeight;
    real macroParticlesPerDt;
    bool propagateV;
    bool accumulate;
    bool split;
    bool join;
    bool logParams;
    bool logHeaderWritten;
    std::string configDumpGeneral();
    std::string toStringGeneral();
    Logger populationlog;
    void updateBaseClassArgs(PopulationArgs args);
};

//! Definition for population dummy construction functions
typedef Population* (*PopulationNewFunc)(PopulationArgs args);

//! Creates population objects
class PopulationFactory
{
public:
    PopulationFactory();
    ~PopulationFactory();
    void init();
    static Population* createPopulation(const std::string typeStr,PopulationArgs args);
    static int getNumberOfPopulations(const std::string typeStr);
    static std::vector<int> getPopulationIds(const std::string typeStr);
private:
    static std::vector<std::string> typeStrTbl;
    static std::vector<PopulationNewFunc> newFuncTbl;
    static std::vector<int> counterTbl;
    static std::vector< std::vector<int> > popIdTbl;
    void registerPopulation(const std::string typeStr,PopulationNewFunc newFunc);
    static void checkVectorSizes();
};

#endif

