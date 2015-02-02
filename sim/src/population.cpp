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
#include "population_solarwind.h"
#include "population_ionospheric.h"
#include "population_exospheric.h"
#include "population_imf.h"
#include "simulation.h"
#include "random.h"
#include "atmosphere.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Find the population id of a population with a given idStr. Returns -1 if no population was found.
int getPopId(string popIdStr)
{
    unsigned int result = -1;
    for(unsigned int i=0; i<Params::pops.size(); ++i) {
        if(Params::pops[i]->getIdStr() == popIdStr) {
            result = Params::pops[i]->getPopId();
        }
    }
    if(result < 0 || result >= Params::pops.size()) {
        result = -1;
    }
    return result;
}

//! Constructor for population argument struct
PopulationArgs::PopulationArgs()
{
    clearArgs();
}

//! Clear population arguments inside PopulationArgs struct
void PopulationArgs::clearArgs()
{
    // TEMPLATE POPULATION VARIABLES
    idStr.value = "";
    idStr.given = false;
    hcFilePrefix.value = "";
    hcFilePrefix.given = false;
    boundaryFUNC.name.clear();
    boundaryFUNC.funcArgs.clear();
    boundaryFUNC.given = false;
    m.value = 0.0;
    m.given = false;
    q.value = 0.0;
    q.given = false;
    T.value = 0.0;
    T.given = false;
    vth.value = 0.0;
    vth.given = false;
    macroParticlesPerDt.value = 0.0;
    macroParticlesPerDt.given = false;
    propagateV.value = 1;
    propagateV.given = false;
    accumulate.value = 1;
    accumulate.given = false;
    split.value = 0.0;
    split.given = false;
    join.value = 0.0;
    join.given = false;
    logParams.value = 0.0;
    logParams.given = false;
    n.value = 0.0;
    n.given = false;
    V.value = 0.0;
    V.given = false;
    backWallWeight.value = 0.0;
    backWallWeight.given = false;
    R.value = 0.0;
    R.given = false;
    totalRate.value = 0.0;
    totalRate.given = false;
    distFunc.name.clear();
    distFunc.funcArgs.clear();
    distFunc.given = false;
#ifdef USE_PARTICLE_SUBCYCLING
    subcycleSteps.value = 0.0;
    subcycleSteps.given = false;
#endif
}

unsigned int Population::idCnt = 0;
vector<string> Population::idStrTbl;
vector<string> Population::hcFilePrefixTbl;

//! Constructor to create particle population object
Population::Population(PopulationArgs args) : m(args.m.value), q(args.q.value)
{
    if(Params::POPULATIONS >= Params::MAX_POPULATIONS) {
        ERRORMSG2("too many particle populations",idStr);
        doabort();
    }
    Params::POPULATIONS++;
    this->popid = idCnt++;
    // Check the necessary information is given
    if(args.idStr.given == false) {
        ERRORMSG2("identification string must be given",idStr);
        doabort();
    }
    if(args.hcFilePrefix.given == false) {
        ERRORMSG2("hcFilePrefix must be given",idStr);
        doabort();
    }
    if(args.boundaryFUNC.given == false) {
        ERRORMSG2("boundaryFUNC must be given",idStr);
        doabort();
    }
    if(args.m.given == false) {
        ERRORMSG2("particle m must be given",idStr);
        doabort();
    }
    if(args.q.given == false) {
        ERRORMSG2("particle q must be given",idStr);
        doabort();
    }
    if(args.T.given == true && args.vth.given == true) {
        ERRORMSG2("only population T or vth can be given",idStr);
        doabort();
    }
    if(args.split.given == false) {
        ERRORMSG2("population split boolean must be given",idStr);
        doabort();
    }
    if(args.join.given == false) {
        ERRORMSG2("population join boolean must be given",idStr);
        doabort();
    }
    saveArgs(args);
    setIdStr(args.idStr.value);
    checkMassAndCharge();
    // Initialize variable values
    hcFilePrefix = "";
    boundaryFuncNames.clear();
    boundaryFuncArgs.clear();
    T = 0;
    vth = 0;
    macroParticleStatisticalWeight = 0;
    macroParticlesPerDt = 0;
    propagateV = true;
    accumulate = true;
    split = false;
    join = false;
    logParams = 0;
    logHeaderWritten = false;
}

//! Dummy constructor
Population::Population() : m(-1.0), q(-1.0) { }

//! Dummy virtual destructor
Population::~Population() { }

//! Dummy implementation for a virtual interface function
void Population::initialize()
{
    WARNINGMSG("dummy implementation function called");
}

//! Dummy implementation for a virtual interface function
void Population::createParticles()
{
    static long cnt = 0;
    if(cnt < 100) {
        WARNINGMSG2("dummy implementation function called",int2string(cnt,2));
        cnt++;
    }
}

//! Dummy implementation for a virtual interface function
void Population::addParticle(shortreal x,shortreal y,shortreal z,real w)
{
    static long cnt = 0;
    if(cnt < 100) {
        WARNINGMSG2("dummy implementation function called",int2string(cnt,2));
        cnt++;
    }
}

//! Dummy implementation for a virtual interface function
real Population::getNeutralDensity(const gridreal[])
{
    static long cnt = 0;
    if(cnt < 100) {
        WARNINGMSG2("dummy implementation function called",int2string(cnt,2));
        cnt++;
    }
    return -1.0;
}

//! Dummy implementation for a virtual interface function
void Population::updateArgs()
{
    WARNINGMSG("dummy implementation function called");
}

//! Dummy implementation for a virtual interface function
void Population::writeExtraHcFile()
{
    WARNINGMSG("dummy implementation function called");
}

//! Dummy implementation for a virtual interface function
string Population::configDump()
{
    WARNINGMSG("dummy implementation function called");
    return "";
}

//! Dummy implementation for a virtual interface function
string Population::toString()
{
    WARNINGMSG("dummy implementation function called");
    return "";
}

void Population::saveArgs(PopulationArgs args)
{
    this->args.clearArgs();
    this->args = args;
}

//! Set population parameters
void Population::updateBaseClassArgs(PopulationArgs args)
{
    this->args = args;
    if(args.hcFilePrefix.given == true) {
        setHcFilePrefix(args.hcFilePrefix.value);
    }
    // Set T and vth
    if(args.vth.given == true && args.T.given == true) {
        ERRORMSG2("only population T or vth can be given",idStr);
        doabort();
    } else if(args.vth.given == false && args.T.given == false) {
        this->T = 0.0;
        this->vth = 0.0;
    } else if(args.T.given == true) {
        if(args.T.value >= 0) {
            this->T = args.T.value;
            // Calculate vth
            if(T >= 0) {
                this->vth = sqrt(Params::k_B*T/m);
            } else {
                ERRORMSG2("cannot calculate vth",idStr);
                doabort();
            }
        } else {
            ERRORMSG2("trying to set T < 0",idStr);
            doabort();
        }
    } else if(args.vth.given == true) {
        if(args.vth.value >= 0) {
            this->vth = args.vth.value;
            // Calculate T
            if(vth >= 0) {
                this->T = sqr(vth)*m/Params::k_B;
            } else {
                ERRORMSG2("cannot calculate T",idStr);
                doabort();
            }
        } else {
            ERRORMSG2("trying to set v_th < 0", idStr);
            doabort();
        }
    }
    if(args.propagateV.given == true) {
        this->propagateV = args.propagateV.value;
    }
    if(args.accumulate.given == true) {
        this->accumulate = args.accumulate.value;
    }
    if(args.split.given == true) {
        this->split = args.split.value;
    }
    if(args.join.given == true) {
        this->join = args.join.value;
    }
    if(args.logParams.given == true) {
        this->logParams = args.logParams.value;
    }
    // Set boundary condition functions
    if(args.boundaryFUNC.given == true) {
        // Check if the functions or parameters have changed
        bool differenceFound = false;
        if( ( boundaryFuncNames.size() == args.boundaryFUNC.name.size() ) &&
            ( boundaryFuncArgs.size() == args.boundaryFUNC.funcArgs.size() ) &&
            (  boundaryFuncNames.size() == boundaryFuncArgs.size() ) ) {
            for(unsigned int i = 0; i < boundaryFuncNames.size(); ++i) {
                // Function names
                if(boundaryFuncNames[i].compare(args.boundaryFUNC.name[i]) != 0) {
                    differenceFound = true;
                }
                for(unsigned int j = 0; j < boundaryFuncArgs[i].size(); ++j) {
                    // Function parameters
                    if(boundaryFuncArgs[i][j] != args.boundaryFUNC.funcArgs[i][j]) {
                        differenceFound = true;
                    }
                }
            }
        } else {
            differenceFound = true;
        }
        // Do possible update
        if(differenceFound == true) {
            this->boundaryFuncNames = args.boundaryFUNC.name;
            this->boundaryFuncArgs = args.boundaryFUNC.funcArgs;
            boundaries = ParticleBoundaryConditions(boundaryFuncNames,boundaryFuncArgs,popid);
        }
    }
#ifdef USE_PARTICLE_SUBCYCLING
    if(args.subcycleSteps.given == true) {
        this->subcycleSteps = args.subcycleSteps.value;
    } else {
        this->subcycleSteps = 10;
    }
#endif
}

//! Check mass and charge of a population
void Population::checkMassAndCharge()
{
    if(m <= 0) {
        ERRORMSG2("trying to set m <= 0",idStr);
        doabort();
    }
    if(q < 0) {
        ERRORMSG2("trying to set q < 0",idStr);
        doabort();
    }
}

//! Set population identification stirng
void Population::setIdStr(const string str)
{
    if(str.compare("") == 0) {
        ERRORMSG("empty population id string");
        doabort();
    }
    // Check whether the given population id string already exists
    for(unsigned int i = 0; i < idStrTbl.size(); ++i) {
        if(idStrTbl[i].compare(str) == 0) {
            ERRORMSG2("given population id string not unique",str);
            doabort();
        }
    }
    this->idStr = str;
    idStrTbl.push_back(str);
}

//! Set hc-file name prefix for the population
void Population::setHcFilePrefix(string str)
{
    if(str.compare("") != 0) {
        this->hcFilePrefix = str;
        // Check if this is initial call
        if( (hcFilePrefixTbl.size() ) < idCnt ) {
            hcFilePrefixTbl.push_back(str);
        } else {
            hcFilePrefixTbl[this->popid] = str;
        }
    } else {
        ERRORMSG2("trying to set empty hc-file prefix",idStr);
        doabort();
    }
}

//! Check particle boundary conditions
bool Population::checkBoundaries(TLinkedParticle& p,fastreal rAverage[3])
{
    return boundaries.checkBoundaries(p,rAverage);
}

//! Return hc-file population configurations
void Population::getHcFileConfigs(vector<string>& filePrefix, vector< vector<int> >& popId)
{
    filePrefix.clear();
    popId.clear();
    // Check the sizes are consistent
    if(hcFilePrefixTbl.size() != idCnt) {
        ERRORMSG("internal error in the class");
        doabort();
    }
    // Construct the prefix table and the corresponding 2d id table
    for(unsigned int i = 0; i < hcFilePrefixTbl.size(); ++i) {
        bool prefixFound = false;
        for(unsigned int j = 0; j < filePrefix.size(); ++j) {
            if(filePrefix[j].compare(hcFilePrefixTbl[i]) == 0) {
                popId[j].push_back(i);
                prefixFound = true;
                break;
            }
        }
        if(prefixFound == false) {
            filePrefix.push_back(hcFilePrefixTbl[i]);
            vector<int> tempA;
            tempA.push_back(i);
            popId.push_back(tempA);
        }
    }
    if( filePrefix.size() != popId.size() ) {
        ERRORMSG("internal error in the function");
        doabort();
    }
}

//! Dump example population config (not implemented)
string Population::configDumpGeneral()
{
    return "";
}

//! Returns string representation of the population
string Population::toStringGeneral()
{
    stringstream ss;
    ss << "popid = " << popid  << "\n";
    ss << "idStr = " << idStr << "\n";
    ss << "hc = " << hcFilePrefix << "\n";
    ss << "logParams = " << logParams << "\n";
    ss << "m = " << m/Params::amu << " amu (" << m << " kg)\n";
    ss << "q = " << q/Params::e << " e (" << q << " C)\n";
    ss << "T = " << T <<" K\n";
    ss << "vth = " << vth/1e3 << " km/s (=sqrt(kB*T/m))\n";
    ss << "boundary conditions =\n";
    ss << "{\n";
    ss << boundaries.toString() << "\n";
    ss << "}\n";
    ss << "macros/dt = " << macroParticlesPerDt  << " #/dt\n";
    ss << "macroweight = " << macroParticleStatisticalWeight << " #\n";
    ss << "propagateV = " << propagateV << "\n";
    ss << "accumulate = " << accumulate << "\n";
    ss << "split = " << split << "\n";
    ss << "join = " << join << "\n";
#ifdef USE_PARTICLE_SUBCYCLING
    ss << "Substepping target steps/cell: " << subcycleSteps << "\n";
#endif
    return ss.str();
}

//! Population factory construction function: solarwind
Population* newPopulationSolarWind(PopulationArgs args)
{
    return new PopulationSolarWind(args);
}

//! Population factory construction function: uniform
Population* newPopulationUniform(PopulationArgs args)
{
    return new PopulationUniform(args);
}

//! Population factory construction function: ionospheric
Population* newPopulationIonospheric(PopulationArgs args)
{
    return new PopulationIonospheric(args);
}

//! Population factory construction function: exospheric
Population* newPopulationExospheric(PopulationArgs args)
{
    return new PopulationExospheric(args);
}

//! Population factory construction function: imf
Population* newPopulationIMF(PopulationArgs args)
{
    return new PopulationIMF(args);
}

vector<string> PopulationFactory::typeStrTbl;
vector<PopulationNewFunc> PopulationFactory::newFuncTbl;
vector<int> PopulationFactory::counterTbl;
vector< vector<int> > PopulationFactory::popIdTbl;

//! Constructor
PopulationFactory::PopulationFactory() { }

//! Initialize population factory
void PopulationFactory::init()
{
    static bool initDone = false;
    if(initDone == false) {
        registerPopulation("uniform",newPopulationUniform);
        registerPopulation("solarwind",newPopulationSolarWind);
        registerPopulation("ionospheric",newPopulationIonospheric);
        registerPopulation("exospheric",newPopulationExospheric);
        registerPopulation("imf",newPopulationIMF);
        initDone = true;
    }
}

//! Destructor
PopulationFactory::~PopulationFactory() { }

//! Creates a new population object
Population* PopulationFactory::createPopulation(const string typeStr,PopulationArgs args)
{
    mainlog << "CREATING A POPULATION [PopulationFactory::createPopulation]: " << args.idStr.value << "\n";
    checkVectorSizes();
    Population* temp = NULL;
    bool typeFound = false;
    for(unsigned int i = 0; i < typeStrTbl.size(); ++i) {
        if(typeStrTbl[i].compare(typeStr) == 0) {
            temp = (newFuncTbl[i])(args);
            typeFound = true;
            counterTbl[i]++;
            popIdTbl[i].push_back( temp->getPopId() );
            break;
        }
    }
    if(typeFound == false) {
        ERRORMSG2("cannot find population type",typeStr);
        doabort();
    }
    return temp;
}

//! Return number of populations of given type
int PopulationFactory::getNumberOfPopulations(const string typeStr)
{
    checkVectorSizes();
    for(unsigned int i = 0; i < typeStrTbl.size(); ++i) {
        if(typeStrTbl[i].compare(typeStr) == 0) {
            return counterTbl[i];
            break;
        }
    }
    return -1;
}

//! Returns a vector of population ids of "typeStr" type populations
vector<int> PopulationFactory::getPopulationIds(const string typeStr)
{
    checkVectorSizes();
    vector<int> tempA;
    for(unsigned int i = 0; i < typeStrTbl.size(); ++i) {
        if(typeStrTbl[i].compare(typeStr) == 0) {
            tempA = popIdTbl[i];
            break;
        }
    }
    return tempA;
}

//! Registers population to the factory
void PopulationFactory::registerPopulation(const string typeStr,PopulationNewFunc newFunc)
{
    checkVectorSizes();
    for(unsigned int i = 0; i < typeStrTbl.size(); ++i) {
        if(typeStrTbl[i].compare(typeStr) == 0) {
            ERRORMSG2("cannot register same population type twice",typeStr);
            doabort();
        }
    }
    typeStrTbl.push_back(typeStr);
    newFuncTbl.push_back(newFunc);
    counterTbl.push_back(0);
    vector<int> tempA;
    popIdTbl.push_back(tempA);
}

//! Check internal consistency in the class
void PopulationFactory::checkVectorSizes()
{
    if( typeStrTbl.size() != newFuncTbl.size() || typeStrTbl.size() != counterTbl.size() || newFuncTbl.size() != counterTbl.size() || counterTbl.size() != popIdTbl.size() ) {
        ERRORMSG("internal error in the class");
        doabort();
    }
}

