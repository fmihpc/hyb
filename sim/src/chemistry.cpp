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

#include "chemistry.h"
#include "definitions.h"

using namespace std;

bool ParticleProcesses::initializedFlag = false;
Processes ParticleProcesses::CX;
Processes ParticleProcesses::EI;

//! Constructor for process argument struct
ProcessArgs::ProcessArgs()
{
    clearArgs();
}

//! Clear population arguments inside PopulationArgs struct
void ProcessArgs::clearArgs()
{
    // TEMPLATE PROCESS VARIABLES
    procIdStr.value = "";
    procIdStr.given = false;
    incidentIonIdStr.value = "";
    incidentIonIdStr.given = false;
    exoNeutralCoronaIdStr.value = "";
    exoNeutralCoronaIdStr.given = false;
    ENAIdStr.value = "";
    ENAIdStr.given = false;
    slowIonIdStr.value = "";
    slowIonIdStr.given = false;
    crossSection.value = 0.0;
    crossSection.given = false;
    weightFactor.value = 0.0;
    weightFactor.given = false;
    slowIonVth.value = 0.0;
    slowIonVth.given = false;
    N_limitHeavyReactions.value = 0.0;
    N_limitHeavyReactions.given = false;
    probLimitHeavyReactions.value = 0.0;
    probLimitHeavyReactions.given = false;
}

//! Constructor
ParticleProcesses::ParticleProcesses() { }

//! Destructor
ParticleProcesses::~ParticleProcesses() { }

//! Create a particle reaction
void ParticleProcesses::createReaction(string processType,ProcessArgs args)
{
    if(processType.compare("ChargeExchange") == 0) {
        initializeReactionChargeExchange(args);
    } else if(processType.compare("ElectronImpactIonization") == 0) {
        initializeReactionElectronImpactIonization(args);
    } else {
        ERRORMSG2("unknown process type",processType);
        doabort();
    }
}

//! Update reaction arguments
void ParticleProcesses::updateReaction(string processType,ProcessArgs args)
{
    if(args.procIdStr.given == false) {
        ERRORMSG("cannot update process, procIdStr not given");
        return;
    }
    // find process id
    int iProc=-1,jProc=-1;
    findProc(args.procIdStr.value,iProc,jProc);
    if(iProc < 0 || jProc < 0) {
        ERRORMSG2("cannot update process, procIdStr not found",args.procIdStr.value);
        return;
    }
    // Update reaction
    if(processType.compare("ChargeExchange") == 0) {
        updateReactionChargeExchange(args,iProc,jProc);
    } else if(processType.compare("ElectronImpactIonization") == 0) {
        updateReactionElectronImpactIonization(args,iProc,jProc);
    } else {
        ERRORMSG2("unknown process type",processType);
        doabort();
    }
}

//! Write particle process log entry
void ParticleProcesses::writeLog()
{
    mainlog
            << "|----------------------------------------------------------------------------- PARTICLE REACTIONS -----------------------------------------------------------------------------|\n"
            << "| Initialized charge exchange reactions:\n";
    if(CX.incidentIonPopId.size() > 0) {
        mainlog << "| Fast ion + neutral -> ENA + slow ion, cross section [m^2], weight factor, slow ion vth [m/s], prob. limit heavy reaction, N_limit heavy reactions, inject ENA, inject slow ion\n";
        for(unsigned int i=0; i<CX.incidentIonPopId.size(); i++) {
            for(unsigned int j=0; j<CX.exoNeutralCoronaPopId[i].size(); j++) {
                mainlog
                        << "| " << CX.processIdStr[i][j] << ": fast(" <<  Params::pops[CX.incidentIonPopId[i]]->getIdStr()
                        << ") + neutral(" << Params::pops[CX.exoNeutralCoronaPopId[i][j]]->getIdStr()
                        << ") -> H-ENA(" << (CX.ENAPopId[i][j] > 0 ? Params::pops[CX.ENAPopId[i][j]]->getIdStr() : "-")
                        << ") + slow(" << (CX.slowIonPopId[i][j] > 0 ? Params::pops[CX.slowIonPopId[i][j]]->getIdStr() : "-")
                        << "), " << CX.crossSection[i][j] << ", "
                        << CX.macroParticleFactor[i][j] << ", "
                        << CX.slowIonVth[i][j] << ", "
                        << CX.probLimitHeavyReactions[i][j] << ", "
                        << CX.N_limitHeavyReactions[i][j] << ", "
                        << CX.injectENA[i][j] << ", "
                        << CX.injectSlowIon[i][j] << "\n";
            }
        }
    }
    mainlog << "| Initialized electron impact ionization reactions:\n";
    if(EI.incidentIonPopId.size() > 0) {
        mainlog << "| Fast electron + neutral -> slow ion, cross section [m^2], weight factor, slow ion vth [m/s], prob. limit heavy reaction, N_limit heavy reactions, inject slow ion\n";
        for(unsigned int i=0; i<EI.incidentIonPopId.size(); i++) {
            for(unsigned int j=0; j<EI.exoNeutralCoronaPopId[i].size(); j++) {
                mainlog
                        << "| " << EI.processIdStr[i][j] << ": fast(" <<  Params::pops[EI.incidentIonPopId[i]]->getIdStr()
                        << ") + neutral(" << Params::pops[EI.exoNeutralCoronaPopId[i][j]]->getIdStr()
                        << ") -> slow(" << (EI.slowIonPopId[i][j] > 0 ? Params::pops[EI.slowIonPopId[i][j]]->getIdStr() : "-")
                        << "), " << EI.crossSection[i][j] << ", "
                        << EI.macroParticleFactor[i][j] << ", "
                        << EI.slowIonVth[i][j] << ", "
                        << EI.probLimitHeavyReactions[i][j] << ", "
                        << EI.N_limitHeavyReactions[i][j] << ", "
                        << EI.injectSlowIon[i][j] << "\n";
            }
        }
    } else {
        mainlog << "| none\n";
    }
    mainlog << "|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|\n";
}

//! Do particle reactions
bool ParticleProcesses::run(TLinkedParticle& part,ParticlePassArgs a)
{
    doChargeExchange(part);
    doElectronImpactIonization(part);
    return true;
}

//! Find if the reaction with a given process id string already exists
bool ParticleProcesses::checkProcIdStrExists(string str)
{
    for(unsigned int i=0; i<CX.processIdStr.size(); i++) {
        for(unsigned int j=0; j<CX.processIdStr[i].size(); j++) {
            if(CX.processIdStr[i][j].compare(str) == 0) {
                return true;
            }
        }
    }
    for(unsigned int i=0; i<EI.processIdStr.size(); i++) {
        for(unsigned int j=0; j<EI.processIdStr[i].size(); j++) {
            if(EI.processIdStr[i][j].compare(str) == 0) {
                return true;
            }
        }
    }
    return false;
}

//! Find the process with given name
void ParticleProcesses::findProc(string str,int& iProc, int& jProc)
{
    for(unsigned int i=0; i<CX.processIdStr.size(); i++) {
        for(unsigned int j=0; j<CX.processIdStr[i].size(); j++) {
            if(CX.processIdStr[i][j].compare(str) == 0) {
                iProc=i;
                jProc=j;
                return;
            }
        }
    }
    for(unsigned int i=0; i<EI.processIdStr.size(); i++) {
        for(unsigned int j=0; j<EI.processIdStr[i].size(); j++) {
            if(EI.processIdStr[i][j].compare(str) == 0) {
                iProc=i;
                jProc=j;
                return;
            }
        }
    }
}

//! Initialize a charge exchange reaction: fast(X+) + neutral(Y) -> ENA(X) + slow(Y+)
void ParticleProcesses::initializeReactionChargeExchange(ProcessArgs args)
{
    if(args.procIdStr.given == false) {
        ERRORMSG("process idStr must be given");
        doabort();
    }
    if(args.incidentIonIdStr.given == false) {
        ERRORMSG2("incident ion idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.exoNeutralCoronaIdStr.given == false) {
        ERRORMSG2("neutral profile idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.ENAIdStr.given == false) {
        ERRORMSG2("ENA idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.slowIonIdStr.given == false) {
        ERRORMSG2("slow ion idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.crossSection.given == false) {
        ERRORMSG2("cross section must be given",args.procIdStr.value);
        doabort();
    }
    if(args.weightFactor.given == false) {
        ERRORMSG2("weight factor must be given",args.procIdStr.value);
        doabort();
    }
    if(args.slowIonVth.given == false) {
        ERRORMSG2("slow ion vth must be given",args.procIdStr.value);
        doabort();
    }
    if(args.N_limitHeavyReactions.given == false) {
        ERRORMSG2("N_limit for heavy reactions must be given",args.procIdStr.value);
        doabort();
    }
    if(args.probLimitHeavyReactions.given == false) {
        ERRORMSG2("heavy reaction probablity limit must be given",args.procIdStr.value);
        doabort();
    }
    mainlog << "CREATING A REACTION [ParticleProcesses::initializeReactionChargeExchange]: " << args.procIdStr.value << "\n";
    bool popIdFound = false;
    int iFound = -1;
    // incident ion
    int popId = getPopId(args.incidentIonIdStr.value);
    if(popId >= 0) {
        for(unsigned int i=0; i<CX.incidentIonPopId.size(); ++i) {
            // popid already in vector
            if(CX.incidentIonPopId[i] == popId) {
                iFound=i;
                popIdFound = true;
                break;
            }
        }
        // popid not yet in vector
        if(popIdFound == false) {
            CX.incidentIonPopId.push_back(popId);
            CX.exoNeutralCoronaPopId.push_back(vector<int>());
            CX.ENAPopId.push_back(vector<int>());
            CX.slowIonPopId.push_back(vector<int>());
            CX.crossSection.push_back(vector<fastreal>());
            CX.macroParticleFactor.push_back(vector<fastreal>());
            CX.weightFactorA.push_back(vector<fastreal>());
            CX.weightFactorB.push_back(vector<fastreal>());
            CX.slowIonVth.push_back(vector<fastreal>());
            CX.N_limitHeavyReactions.push_back(vector<real>());
            CX.probLimitHeavyReactions.push_back(vector<real>());
            CX.heavyReactionCounter.push_back(vector<real>());
            CX.injectENA.push_back(vector<bool>());
            CX.injectSlowIon.push_back(vector<bool>());
            CX.processIdStr.push_back(vector<string>());
            iFound = CX.incidentIonPopId.size()-1;
        }
    } else {
        ERRORMSG2("ParticleProcesses::initializeReaction: bad incident ion idStr",args.incidentIonIdStr.value);
        doabort();
    }
    // neutral profile
    popId = getPopId(args.exoNeutralCoronaIdStr.value);
    if(popId >= 0) {
        CX.exoNeutralCoronaPopId[iFound].push_back(popId);
    } else {
        ERRORMSG2("ParticleProcesses::initializeReaction: bad exospheric neutral corona idStr",args.exoNeutralCoronaIdStr.value);
        doabort();
    }
    // ENA
    if(args.ENAIdStr.value.compare("-") == 0) {
        CX.injectENA[iFound].push_back(false);
        CX.ENAPopId[iFound].push_back(-1);
    } else {
        popId = getPopId(args.ENAIdStr.value);
        if(popId >= 0) {
            CX.ENAPopId[iFound].push_back(popId);
            CX.injectENA[iFound].push_back(true);
        } else {
            ERRORMSG2("ParticleProcesses::initializeReaction: bad ENA idStr",args.ENAIdStr.value);
            doabort();
        }
    }
    // slow ion
    if(args.slowIonIdStr.value.compare("-") == 0) {
        CX.injectSlowIon[iFound].push_back(false);
        CX.slowIonPopId[iFound].push_back(-1);
    } else {
        popId = getPopId(args.slowIonIdStr.value);
        if(popId >= 0) {
            CX.slowIonPopId[iFound].push_back(popId);
            CX.injectSlowIon[iFound].push_back(true);
        } else {
            ERRORMSG2("ParticleProcesses::initializeReaction: bad slow ion idStr",args.slowIonIdStr.value);
            doabort();
        }
    }
    // cross section
    if(args.crossSection.value > 0) {
        CX.crossSection[iFound].push_back(args.crossSection.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: cross section < 0");
        doabort();
    }
    // macroparticle factor
    if(args.weightFactor.value > 0) {
        CX.macroParticleFactor[iFound].push_back(args.weightFactor.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: macroparticle factor < 0");
        doabort();
    }
    if(args.slowIonVth.value > 0) {
        CX.slowIonVth[iFound].push_back(args.slowIonVth.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: slowIonVth < 0");
        doabort();
    }
    CX.weightFactorA[iFound].push_back(1.0/args.weightFactor.value);
    CX.weightFactorB[iFound].push_back(1.0 - 1.0/args.weightFactor.value);
    // heavy reaction limit
    CX.N_limitHeavyReactions[iFound].push_back(args.N_limitHeavyReactions.value);
    CX.probLimitHeavyReactions[iFound].push_back(args.probLimitHeavyReactions.value);
    CX.heavyReactionCounter[iFound].push_back(0.0);
    CX.processIdStr[iFound].push_back(args.procIdStr.value);
    initializedFlag = true;
}

//! Update a charge exchange reaction parameters
void ParticleProcesses::updateReactionChargeExchange(ProcessArgs args,const int iProc,const int jProc)
{
    if(args.incidentIonIdStr.given == true) {
        int popId = getPopId(args.incidentIonIdStr.value);
        if(popId >= 0) {
            CX.incidentIonPopId[iProc] = popId;
        } else {
            ERRORMSG2("popid not found",CX.processIdStr[iProc][jProc]);
        }
    }
    if(args.exoNeutralCoronaIdStr.given == true) {
        int popId = getPopId(args.exoNeutralCoronaIdStr.value);
        if(popId >= 0) {
            CX.exoNeutralCoronaPopId[iProc][jProc] = popId;
        } else {
            ERRORMSG2("popid not found",CX.processIdStr[iProc][jProc]);
        }
    }
    if(args.ENAIdStr.given == true) {
        if(args.ENAIdStr.value.compare("-") == 0) {
            CX.injectENA[iProc][jProc] = false;
            CX.ENAPopId[iProc][jProc] = -1;
        } else {
            int popId = getPopId(args.ENAIdStr.value);
            if(popId >= 0) {
                CX.ENAPopId[iProc][jProc] = popId;
                CX.injectENA[iProc][jProc] = true;
            } else {
                ERRORMSG2("popid not found",CX.processIdStr[iProc][jProc]);
            }
        }
    }
    if(args.slowIonIdStr.given == true) {
        if(args.slowIonIdStr.value.compare("-") == 0) {
            CX.injectSlowIon[iProc][jProc] = false;
            CX.slowIonPopId[iProc][jProc] = -1;
        } else {
            int popId = getPopId(args.slowIonIdStr.value);
            if(popId >= 0) {
                CX.slowIonPopId[iProc][jProc] = popId;
                CX.injectSlowIon[iProc][jProc] = true;
            } else {
                ERRORMSG2("popid not found",CX.processIdStr[iProc][jProc]);
            }
        }
    }
    if(args.crossSection.given == true) {
        if(args.crossSection.value > 0) {
            CX.crossSection[iProc][jProc] = args.crossSection.value;
        } else {
            ERRORMSG2("cross section < 0, not updating value",CX.processIdStr[iProc][jProc]);
        }
    }
    if(args.weightFactor.given == true) {
        if(args.weightFactor.value > 0) {
            CX.macroParticleFactor[iProc][jProc] = args.weightFactor.value;
        } else {
            ERRORMSG2("weight factor < 0, not updating value",CX.processIdStr[iProc][jProc]);
        }
    }
    if(args.slowIonVth.given == true) {
        if(args.slowIonVth.value > 0) {
            CX.slowIonVth[iProc][jProc] = args.slowIonVth.value;
        } else {
            ERRORMSG2("slow ion vth < 0, not updating value",CX.processIdStr[iProc][jProc]);
        }
    }
    if(args.N_limitHeavyReactions.given == true) {
        CX.N_limitHeavyReactions[iProc][jProc] = args.N_limitHeavyReactions.value;
    }
    if(args.probLimitHeavyReactions.given == true) {
        CX.probLimitHeavyReactions[iProc][jProc] = args.probLimitHeavyReactions.value;
    }
}

//! Do charge exchange
void ParticleProcesses::doChargeExchange(TLinkedParticle& part)
{
    const gridreal r[3] = {part.x, part.y, part.z};
    const fastreal vdt = sqrt(sqr(part.vx) + sqr(part.vy) + sqr(part.vz))*Params::dt;
    for(unsigned int i=0; i<CX.incidentIonPopId.size(); i++) {
        //check if the particle participates in the charge exchange
        if(part.popid == CX.incidentIonPopId[i]) {
            for(unsigned int j=0; j<CX.exoNeutralCoronaPopId[i].size(); j++) {
                fastreal neutralDensity=Params::pops[CX.exoNeutralCoronaPopId[i][j]]->getNeutralDensity(r);
                fastreal prob = neutralDensity*CX.crossSection[i][j]*vdt;
                // prob cannot be large, otherwise the probability of next loop is not correct (P2*(1-P1)~P2 when P1<<1)
                if(prob > CX.probLimitHeavyReactions[i][j]) {
                    CX.heavyReactionCounter[i][j] += 1.0;
                    errorlog << "ChargeExchange " << CX.processIdStr[i][j] << ": heavy reaction happens, counter = " << CX.heavyReactionCounter[i][j] << endl;
                    if(CX.heavyReactionCounter[i][j] > CX.N_limitHeavyReactions[i][j]) {
                        ERRORMSG2("ChargeExchange: too many heavy reactions happened. Reduce the time step!",CX.processIdStr[i][j]);
                        doabort();
                    }
                }
                // Charge exchange happens
                if(uniformrnd() < prob*CX.macroParticleFactor[i][j]) {
                    // inject ENA
                    if(CX.injectENA[i][j] == true) {
                        g.addparticle(part.x,part.y,part.z,part.vx,part.vy,part.vz,part.w*CX.weightFactorA[i][j],CX.ENAPopId[i][j]);
                    }
                    // inject slow ion
                    if(CX.injectSlowIon[i][j] == true) {
                        const fastreal vx = CX.slowIonVth[i][j]*gaussrnd();
                        const fastreal vy = CX.slowIonVth[i][j]*gaussrnd();
                        const fastreal vz = CX.slowIonVth[i][j]*gaussrnd();
                        g.addparticle(part.x,part.y,part.z,vx,vy,vz,part.w*CX.weightFactorA[i][j],CX.slowIonPopId[i][j]);
                    }
                    // update weight of the original ion
                    part.w *= CX.weightFactorB[i][j];
#ifndef NO_DIAGNOSTICS
                    Params::diag.pCounter[part.popid]->chargeExchangeRate += 1.0;
#endif
                }
            }
            return;
        }
    }
}

//! Initialize an electron impact ionization reaction: fast(X+) + neutral(Y) -> slow(Y+)
void ParticleProcesses::initializeReactionElectronImpactIonization(ProcessArgs args)
{
    if(args.procIdStr.given == false) {
        ERRORMSG("process idStr must be given");
        doabort();
    }
    if(args.incidentIonIdStr.given == false) {
        ERRORMSG2("incident ion idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.exoNeutralCoronaIdStr.given == false) {
        ERRORMSG2("neutral profile idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.ENAIdStr.given == true) {
        ERRORMSG2("ENA idStr cannot be given",args.procIdStr.value);
        doabort();
    }
    if(args.slowIonIdStr.given == false) {
        ERRORMSG2("slow ion idStr must be given",args.procIdStr.value);
        doabort();
    }
    if(args.crossSection.given == false) {
        ERRORMSG2("cross section must be given",args.procIdStr.value);
        doabort();
    }
    if(args.weightFactor.given == false) {
        ERRORMSG2("weight factor must be given",args.procIdStr.value);
        doabort();
    }
    if(args.slowIonVth.given == false) {
        ERRORMSG2("slow ion vth must be given",args.procIdStr.value);
        doabort();
    }
    if(args.N_limitHeavyReactions.given == false) {
        ERRORMSG2("N_limit for heavy reactions must be given",args.procIdStr.value);
        doabort();
    }
    if(args.probLimitHeavyReactions.given == false) {
        ERRORMSG2("heavy reaction probablity limit must be given",args.procIdStr.value);
        doabort();
    }
    mainlog << "CREATING A REACTION [ParticleProcesses::initializeReactionElectronImpactIonization]: " << args.procIdStr.value << "\n";
    bool popIdFound = false;
    int iFound = -1;
    // incident ion
    int popId = getPopId(args.incidentIonIdStr.value);
    if(popId >= 0) {
        for(unsigned int i=0; i<EI.incidentIonPopId.size(); ++i) {
            // popid already in vector
            if(EI.incidentIonPopId[i] == popId) {
                iFound=i;
                popIdFound = true;
                break;
            }
        }
        // popid not yet in vector
        if(popIdFound == false) {
            EI.incidentIonPopId.push_back(popId);
            EI.exoNeutralCoronaPopId.push_back(vector<int>());
            EI.slowIonPopId.push_back(vector<int>());
            EI.crossSection.push_back(vector<fastreal>());
            EI.macroParticleFactor.push_back(vector<fastreal>());
            EI.weightFactorA.push_back(vector<fastreal>());
            EI.weightFactorB.push_back(vector<fastreal>());
            EI.slowIonVth.push_back(vector<fastreal>());
            EI.N_limitHeavyReactions.push_back(vector<real>());
            EI.probLimitHeavyReactions.push_back(vector<real>());
            EI.heavyReactionCounter.push_back(vector<real>());
            EI.injectSlowIon.push_back(vector<bool>());
            EI.processIdStr.push_back(vector<string>());
            iFound = EI.incidentIonPopId.size()-1;
        }
    } else {
        ERRORMSG2("ParticleProcesses::initializeReaction: bad incident ion idStr",args.incidentIonIdStr.value);
        doabort();
    }
    // neutral profile
    popId = getPopId(args.exoNeutralCoronaIdStr.value);
    if(popId >= 0) {
        EI.exoNeutralCoronaPopId[iFound].push_back(popId);
    } else {
        ERRORMSG2("ParticleProcesses::initializeReaction: bad exospheric neutral corona idStr",args.exoNeutralCoronaIdStr.value);
        doabort();
    }
    // slow ion
    if(args.slowIonIdStr.value.compare("-") == 0) {
        EI.injectSlowIon[iFound].push_back(false);
        EI.slowIonPopId[iFound].push_back(-1);
    } else {
        popId = getPopId(args.slowIonIdStr.value);
        if(popId >= 0) {
            EI.slowIonPopId[iFound].push_back(popId);
            EI.injectSlowIon[iFound].push_back(true);
        } else {
            ERRORMSG2("ParticleProcesses::initializeReaction: bad slow ion idStr",args.slowIonIdStr.value);
            doabort();
        }
    }
    // cross section
    if(args.crossSection.value > 0) {
        EI.crossSection[iFound].push_back(args.crossSection.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: cross section < 0");
        doabort();
    }
    // macroparticle factor
    if(args.weightFactor.value > 0) {
        EI.macroParticleFactor[iFound].push_back(args.weightFactor.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: macroparticle factor < 0");
        doabort();
    }
    if(args.slowIonVth.value > 0) {
        EI.slowIonVth[iFound].push_back(args.slowIonVth.value);
    } else {
        ERRORMSG("ParticleProcesses::initializeReaction: slowIonVth < 0");
        doabort();
    }
    EI.weightFactorA[iFound].push_back(1.0/args.weightFactor.value);
    EI.weightFactorB[iFound].push_back(1.0 - 1.0/args.weightFactor.value);
    // heavy reaction limit
    EI.N_limitHeavyReactions[iFound].push_back(args.N_limitHeavyReactions.value);
    EI.probLimitHeavyReactions[iFound].push_back(args.probLimitHeavyReactions.value);
    EI.heavyReactionCounter[iFound].push_back(0.0);
    EI.processIdStr[iFound].push_back(args.procIdStr.value);
    initializedFlag = true;
}

//! Update an eletron impact ionization reaction parameters
void ParticleProcesses::updateReactionElectronImpactIonization(ProcessArgs args,const int iProc,const int jProc)
{
    if(args.incidentIonIdStr.given == true) {
        int popId = getPopId(args.incidentIonIdStr.value);
        if(popId >= 0) {
            EI.incidentIonPopId[iProc] = popId;
        } else {
            ERRORMSG2("popid not found",EI.processIdStr[iProc][jProc]);
        }
    }
    if(args.exoNeutralCoronaIdStr.given == true) {
        int popId = getPopId(args.exoNeutralCoronaIdStr.value);
        if(popId >= 0) {
            EI.exoNeutralCoronaPopId[iProc][jProc] = popId;
        } else {
            ERRORMSG2("popid not found",EI.processIdStr[iProc][jProc]);
        }
    }
    if(args.ENAIdStr.given == true) {
        ERRORMSG2("ENA idStr cannot be given, not updating",EI.processIdStr[iProc][jProc]);
    }
    if(args.slowIonIdStr.given == true) {
        if(args.slowIonIdStr.value.compare("-") == 0) {
            EI.injectSlowIon[iProc][jProc] = false;
            EI.slowIonPopId[iProc][jProc] = -1;
        } else {
            int popId = getPopId(args.slowIonIdStr.value);
            if(popId >= 0) {
                EI.slowIonPopId[iProc][jProc] = popId;
                EI.injectSlowIon[iProc][jProc] = true;
            } else {
                ERRORMSG2("popid not found",EI.processIdStr[iProc][jProc]);
            }
        }
    }
    if(args.crossSection.given == true) {
        if(args.crossSection.value > 0) {
            EI.crossSection[iProc][jProc] = args.crossSection.value;
        } else {
            ERRORMSG2("cross section < 0, not updating value",EI.processIdStr[iProc][jProc]);
        }
    }
    if(args.weightFactor.given == true) {
        if(args.weightFactor.value > 0) {
            EI.macroParticleFactor[iProc][jProc] = args.weightFactor.value;
        } else {
            ERRORMSG2("weight factor < 0, not updating value",EI.processIdStr[iProc][jProc]);
        }
    }
    if(args.slowIonVth.given == true) {
        if(args.slowIonVth.value > 0) {
            EI.slowIonVth[iProc][jProc] = args.slowIonVth.value;
        } else {
            ERRORMSG2("slow ion vth < 0, not updating value",EI.processIdStr[iProc][jProc]);
        }
    }
    if(args.N_limitHeavyReactions.given == true) {
        EI.N_limitHeavyReactions[iProc][jProc] = args.N_limitHeavyReactions.value;
    }
    if(args.probLimitHeavyReactions.given == true) {
        EI.probLimitHeavyReactions[iProc][jProc] = args.probLimitHeavyReactions.value;
    }
}

//! Do electron impact ionization
void ParticleProcesses::doElectronImpactIonization(TLinkedParticle& part)
{
    const gridreal r[3] = {part.x, part.y, part.z};
    for(unsigned int i=0; i<EI.incidentIonPopId.size(); i++) {
        //check if the particle participates in the electron impact ionization
        if(part.popid == EI.incidentIonPopId[i]) {
            for(unsigned int j=0; j<EI.exoNeutralCoronaPopId[i].size(); j++) {
                fastreal neutralDensity=Params::pops[EI.exoNeutralCoronaPopId[i][j]]->getNeutralDensity(r);
                fastreal prob = neutralDensity*EI.crossSection[i][j]*Params::dt;
                // prob cannot be large, otherwise the probability of next loop is not correct (P2*(1-P1)~P2 when P1<<1)
                if(prob > EI.probLimitHeavyReactions[i][j]) {
                    EI.heavyReactionCounter[i][j] += 1.0;
                    errorlog << "ElectronImpactIonization " << EI.processIdStr[i][j] << ": heavy reaction, counter = " << EI.heavyReactionCounter[i][j] << endl;
                    if(EI.heavyReactionCounter[i][j] > EI.N_limitHeavyReactions[i][j]) {
                        ERRORMSG2("ElectronImpactIonization: too many heavy reactions happened. Reduce the time step!",EI.processIdStr[i][j]);
                        doabort();
                    }
                }
                // Electron impact ionization happens
                if(uniformrnd() < prob*EI.macroParticleFactor[i][j]) {
                    if(EI.injectSlowIon[i][j] == true) {
                        const fastreal vx = EI.slowIonVth[i][j]*gaussrnd();
                        const fastreal vy = EI.slowIonVth[i][j]*gaussrnd();
                        const fastreal vz = EI.slowIonVth[i][j]*gaussrnd();
                        g.addparticle(part.x,part.y,part.z,vx,vy,vz,part.w*EI.weightFactorA[i][j],EI.slowIonPopId[i][j]);
                    }
#ifndef NO_DIAGNOSTICS
                    Params::diag.pCounter[part.popid]->electronImpactIonizationRate += 1.0;
#endif
                }
            }
            return;
        }
    }
}

