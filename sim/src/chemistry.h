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

#ifndef CHEMISTRY_H
#define CHEMISTRY_H

#include "simulation.h"

//! All process arguments
struct ProcessArgs {
    // TEMPLATE PROCESS VARIABLES
    stringArg procIdStr;
    stringArg incidentIonIdStr;
    stringArg exoNeutralCoronaIdStr;
    stringArg ENAIdStr;
    stringArg slowIonIdStr;
    fastrealArg crossSection;
    fastrealArg weightFactor;
    fastrealArg slowIonVth;
    realArg N_limitHeavyReactions;
    realArg probLimitHeavyReactions;
    ProcessArgs();
    void clearArgs();
};

//! Process variables
struct Processes {
    std::vector<int> incidentIonPopId;
    std::vector< std::vector<int> > exoNeutralCoronaPopId, ENAPopId, slowIonPopId;
    std::vector< std::vector<fastreal> > crossSection, macroParticleFactor, weightFactorA, weightFactorB, slowIonVth;
    std::vector< std::vector<real> > probLimitHeavyReactions,N_limitHeavyReactions,heavyReactionCounter;
    std::vector< std::vector<bool> > injectENA,injectSlowIon;
    std::vector< std::vector<std::string> > processIdStr;
};

//! Particle processes such as charge exchange and electron impact ionization
class ParticleProcesses
{
public:
    ParticleProcesses();
    ~ParticleProcesses();
    static void writeLog();
    static void createReaction(std::string processType,ProcessArgs args);
    static void updateReaction(std::string processType,ProcessArgs args);
    static bool run(TLinkedParticle& part,ParticlePassArgs a);
    static bool checkProcIdStrExists(std::string str);
    static bool isInitialized() {
        return initializedFlag;
    }
private:
    static Processes CX; //!< Charge exchange processes
    static Processes EI; //!< Electron impact ionization processes
    static bool initializedFlag; //!< If the class is initialized
    static void findProc(std::string str,int& iProc, int& jProc);
    static void initializeReactionChargeExchange(ProcessArgs args);
    static void updateReactionChargeExchange(ProcessArgs args,const int iProc,const int jProc);
    static void doChargeExchange(TLinkedParticle& part);
    static void initializeReactionElectronImpactIonization(ProcessArgs args);
    static void updateReactionElectronImpactIonization(ProcessArgs args,const int iProc,const int jProc);
    static void doElectronImpactIonization(TLinkedParticle& part);
};

#endif
