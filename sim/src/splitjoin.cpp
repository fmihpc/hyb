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

#include <cmath>
#include <cstdlib>
#include <sstream>
#include "splitjoin.h"
#include "simulation.h"
#include "params.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Params simuConfig;

//! Default constructor
Split::Split()
{
    // Make sure the class is not used without proper initialization
    this->ptr = &Split::defaultFunction;
    resetParameters();
}

#define ELSEIF_SPLIT(func) else if(funcName.compare(#func) == 0) { this->ptr = &Split::func; setArgs_ ## func(); }

//! Constructor
Split::Split(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &Split::defaultFunction;
    this->args = args;
    resetParameters();
    // SPLIT METHODS
    if(funcName.compare("") == 0)  {
        ERRORMSG("empty split function name");
        doabort();
    }
    ELSEIF_SPLIT(splitDefault)
    ELSEIF_SPLIT(splitOriginal)
    else {
        ERRORMSG2("bad split function name",funcName);
        doabort();
    }
}

//! Destructor
Split::~Split() { }

//! Split macroparticles in a given particle list
int Split::doSplitting(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist)
{
    return (this->*ptr)(boxmin,boxmax,tplist);
}

//! Default function, which aborts the program if called
int Split::defaultFunction(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist)
{
    ERRORMSG("function pointer not set");
    doabort();
    return 0;
}

/** \brief Split macroparticles 1->2
 *
 * Split particle old to produce particles A and B. The daughter
 * particles A,B have the same velocity. They are set apart by a
 * randomly directed distance whose length is split_distance/2 from
 * the old position. The weight is distributed evenly among A and B,
 * and (of course) A and B both have the same mass m1.
 */
bool Split::split12(TLinkedParticle& p, gridreal split_distance)
{
    // produce a random vector dx, whose length is |dx| = split_distance/2
    // and which is perpendicular to old.v
    const Tgr3v v(p.vx,p.vy,p.vz);
    Tgr3v dx = RandomPerpendicularUnitVector(v)*(0.5*split_distance);
    shortreal rNewA[3],rNewB[3];
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    rNewA[0] = p.x - dx[0];
    rNewA[1] = p.y - dx[1];
    rNewA[2] = p.z - dx[2];
    rNewB[0] = p.x + dx[0];
    rNewB[1] = p.y + dx[1];
    rNewB[2] = p.z + dx[2];
#else
    gridreal r_old[3] = {p.x,p.y,p.z};
    // Here is very important the sequence of transformation
    sph_transf_H2S_R(r_old);
    sph_transf_S2C_R(r_old);
    for (int d=0; d<3; d++) {
        rNewA[d] = r_old[d] - dx[d];
        rNewB[d] = r_old[d] + dx[d];
    }
    // Here is very important the sequence of transformation
    sph_transf_C2S_r(rNewA);
    sph_transf_C2S_r(rNewB);
    sph_transf_S2H_R(rNewA);
    sph_transf_S2H_R(rNewB);
#endif

    if(Params::insideBoxTight(rNewA) == false ||
       Params::insideBoxTight(rNewB) == false) {
        return false;
    }
    p.x = rNewA[0];
    p.y = rNewA[1];
    p.z = rNewA[2];
    p.w = 0.5*p.w;
    g.addparticle(rNewB[0],rNewB[1],rNewB[2],p.vx,p.vy,p.vz,p.w,p.popid,false);
#ifndef NO_DIAGNOSTICS
    // Increase counter
    Params::diag.pCounter[p.popid]->splittingRate += 1;
#endif
    return true;
}

//! Reset split&join parameters
void Split::resetParameters()
{
    // SPLIT PARAMETERS
    distanceFactor = 0.0;
    HALFCELL_SPLIT_DISTANCE = true;
}

// SPLIT METHODS

//! Set function arguments
void Split::setArgs_splitDefault()
{
    if(args.size() != 1) {
        ERRORMSG2("function takes one argument",name);
        doabort();
    }
    distanceFactor = args[0];
}

//! Split: choose a random population and find max(part.w) for that population
int Split::splitDefault(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist)
{
    //choose the population at random and
    // select the "weightiest" particle for splitting from the chosen population
    // now done in one loop over macro particles
    int poppi, kk=0, kp, choose[Params::POPULATIONS];
    fastreal maxweight[Params::POPULATIONS];
    TLinkedParticle *p, *maxweight_p[Params::POPULATIONS];
    for (poppi=0; poppi<Params::POPULATIONS; poppi++) {
        maxweight_p[poppi] = NULL;
        maxweight[poppi] = 0;
    }
    for (p=tplist.first; p; p=p->next) {
        if (Params::pops[p->popid]->getSplit() && p->w > maxweight[p->popid]) {
            maxweight_p[p->popid] = p;
            maxweight[p->popid] = p->w;
        }
    }
    //choose an available population at random
    for (poppi=0; poppi<Params::POPULATIONS; poppi++) {
        if (maxweight[poppi] > 0) {
            choose[kk++] = poppi;
        }
    }
    if (kk==0) {
        return 0;
    }
    kp = int(uniformrnd() * kk);
    if (kp==kk) {
        kp--;
        errorlog << "SplitDefault: uniformrnd() gave 1.00000" << endl;
    }
    // just to make sure in the (no idea how) unlikely event: uniformrnd()==1
    const int chosen_pop = choose[kp];
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    // Original spliting in hybrid coordinates (see Split12 in Grid.cpp)
    const gridreal size = boxmax[0] - boxmin[0];
#else
    // Splitting in Cartesian coordinates (see sph_Split12 in Grid.cpp)
    // To split particle we need to know minimal line element of the cell
    gridreal min_size[3];
    for (int d=0; d<3; d++) min_size[d] = boxmax[d] - boxmin[d];
    const gridreal size = min3(min_size[0], min_size[1], min_size[2]);
#endif
    const bool splitOk = split12(*maxweight_p[chosen_pop],distanceFactor*size);
    if(splitOk == true) {
        return 1;
    } else {
        return 0;
    }
}

//! Set function arguments
void Split::setArgs_splitOriginal()
{
    if(args.size() != 1) {
        ERRORMSG2("function takes one argument",name);
        doabort();
    }
    if(args[0] == 0) {
        HALFCELL_SPLIT_DISTANCE = false;
        distanceFactor = 0.1;
    } else if(args[0] == 1) {
        HALFCELL_SPLIT_DISTANCE = true;
        distanceFactor = 0.5;
    } else {
        ERRORMSG2("function takes only argument 0 or 1",name);
        doabort();
    }
}

//! Original Split function - chooses the particle with the Highest Absolute Weight!
int Split::splitOriginal(const gridreal minbox[3], const gridreal maxbox[3], TParticleList& tplist)
{
    TLinkedParticle *p;
    // select the "weightiest" particle for splitting (note heaviest is
    // the same as weightiest only if one species)
    shortreal maxweight = 0;
    TLinkedParticle *maxweight_p = 0;
    for (p=tplist.first; p; p=p->next) {
        if (!Params::insideBoxTight(p)) {
            continue;
        }
        if (Params::pops[p->popid]->getSplit() && p->w > maxweight) {
            maxweight_p = p;
            maxweight = p->w;
        }
    }
    if (maxweight_p == 0) {
        // cancel split if none was found in tightbox (rather unlikely, but can occur)
        return 0;
    }
    const gridreal size = maxbox[0] - minbox[0];
    const bool splitOk = split12(*maxweight_p,distanceFactor*size);
    if(splitOk == true) {
        return 1;
    } else {
        return 0;
    }
}

//! Default constructor
Join::Join()
{
    // Make sure the class is not used without proper initialization
    this->ptr = &Join::defaultFunction;
    resetParameters();
}

#define ELSEIF_JOIN(func) else if(funcName.compare(#func) == 0) { this->ptr = &Join::func; setArgs_ ## func(); }

//! Constructor
Join::Join(string funcName,vector<real> args)
{
    this->name = funcName;
    this->ptr = &Join::defaultFunction;
    this->args = args;
    resetParameters();
    // JOIN METHODS
    if(funcName.compare("") == 0)  {
        ERRORMSG("empty join function name");
        doabort();
    }
    ELSEIF_JOIN(joinDefault)
    ELSEIF_JOIN(joinOriginal)
    else {
        ERRORMSG2("bad join function name",funcName);
        doabort();
    }
}

//! Destructor
Join::~Join() { }

//! Join macroparticles
int Join::doJoining(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int fast)
{
    return (this->*ptr)(boxmin,boxmax,tplist,fast);
}

//! Default function, which aborts the program if called
int Join::defaultFunction(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int dummy)
{
    ERRORMSG("function pointer not set");
    doabort();
    return 0;
}

/** \brief Join macroparticles 3->2
 *
 * After splitting, the daughter particles A and B have exactly the
 * same velocity. If multiple splitting occur, we may have a group of
 * particles with the same velocity. If such a group is joined later
 * on, the angular momentum L and kinetic energy Wk relative to CM
 * are zero, and the joined particles fall exactly at the same
 * coordinates, the velocities also being the same. This can happen
 * easily in a region where external forces vanish, because then the
 * velocities do not get modified by the dynamics. This looks bizarre
 * but should not do any harm. The worst that happens is that we have
 * some degenerate cloned  particles. Once external forces start
 * acting, the formation of these stops, and subsequent join
 * operations find these and coalesce into one.
 */
void Join::join32(const TLinkedParticle& P1, const TLinkedParticle& P2, const TLinkedParticle& P3, TLinkedParticle& A, TLinkedParticle& B)
{
    const fastreal w1 = P1.w, w2 = P2.w, w3 = P3.w;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    Tgr3v r1(P1.x,P1.y,P1.z);
    Tgr3v r2(P2.x,P2.y,P2.z);
    Tgr3v r3(P3.x,P3.y,P3.z);
    Tgr3v v1(P1.vx,P1.vy,P1.vz);
    Tgr3v v2(P2.vx,P2.vy,P2.vz);
    Tgr3v v3(P3.vx,P3.vy,P3.vz);
#else
    gridreal r1_old[3] = {P1.x, P1.y, P1.z};
    gridreal r2_old[3] = {P2.x, P2.y, P2.z};
    gridreal r3_old[3] = {P3.x, P3.y, P3.z};
    gridreal v1_old[3] = {P1.vx, P1.vy, P1.vz};
    gridreal v2_old[3] = {P2.vx, P2.vy, P2.vz};
    gridreal v3_old[3] = {P3.vx, P3.vy, P3.vz};
    gridreal r1_new[3], r2_new[3];
    gridreal v1_new[3], v2_new[3];
    sph_transf_H2S_R(r1_old);
    sph_transf_H2S_R(r2_old);
    sph_transf_H2S_R(r3_old);
    sph_transf_S2C_R(r1_old);
    sph_transf_S2C_R(r2_old);
    sph_transf_S2C_R(r3_old);
    Tgr3v r1(r1_old[0], r1_old[1], r1_old[2]);
    Tgr3v r2(r2_old[0], r2_old[1], r2_old[2]);
    Tgr3v r3(r3_old[0], r3_old[1], r3_old[2]);
    Tgr3v v1(v1_old[0], v1_old[1], v1_old[2]);
    Tgr3v v2(v2_old[0], v2_old[1], v2_old[2]);
    Tgr3v v3(v3_old[0], v3_old[1], v3_old[2]);
#endif
    // Total weight of the input particles. Total mass is Wtot*pop_m[id].
    const fastreal Wtot = w1 + w2 + w3;
    // Total "momentum" of the input particles. Momentum is actually
    // this multiplied by the (physical) mass of the particles.
    const Tgr3v Ptot = w1*v1 + w2*v2 + w3*v3;
    // Center of mass velocity of the input particles
    const Tgr3v vCM = Ptot/Wtot;
    // Center of mass position of the input particles
    const Tgr3v rCM = (w1*r1 + w2*r2 + w3*r3)/Wtot;
    // Velocity and coordinate transform to CM
    r1 -= rCM;
    r2-= rCM;
    r3-= rCM;
    v1 -= vCM;
    v2-= vCM;
    v3-= vCM;
    // Total kinetic "energy" of the input particles in CM (actually
    // should be multiplied by the physical mass to the correct
    // dimension)
    const fastreal Wk = 0.5*(w1*magn2(v1) + w2*magn2(v2) + w3*magn2(v3));
    // Total angular "momentum" of the input particles in CM (actually
    // should be multiplied by the physical mass to the correct
    // dimension)
    const Tgr3v L = w1*Cross(r1,v1) + w2*Cross(r2,v2) + w3*Cross(r3,v3);
    // Both output particles have the same weight
    const fastreal w = 0.5*Wtot;
    // n is the normal unit vector of the plane formed by the three
    // input particles
    const Tgr3v nvec = Cross(r2-r1,r3-r1);
    const Tgr3v n = UnitVector(nvec);
    // Select v so that |v|=sqrt(Wk/w), v.L == 0 and v is as much
    // n-directed as possible
    const Tgr3v uL = UnitVector(L);
    const fastreal vmagn = sqrt(Wk/w);
    // Now v is projection of n to plane perpendicular to L
    Tgr3v v = vmagn*UnitVector_robust(n - uL*dot(uL,n),uL);
    fastreal v_dot_n = dot(v,n);
    // Ensure that v**n >= 0 ==> eases division by it below
    if (v_dot_n < 0) {
        v = -v;
        v_dot_n = -v_dot_n;
    }
    fastreal W_epsilon = 1e-15*P1.w*(sqr(P1.vx) + sqr(P1.vy) + sqr(P1.vz));
    if (W_epsilon == 0) {
        W_epsilon = 1e-30;
    }
    // Now 2*m*Cross(r,v) == L, i.e. angular momentum is conserved
    Tgr3v r = Cross(v,L)/(W_epsilon + 2*Wk);
    fastreal v_epsilon = vmagn*1e-7;
    if (v_epsilon == 0) {
        v_epsilon = 1e-30;
    }
    r -= v*dot(n,r)/(v_epsilon + v_dot_n);
    //sets new particles in the plane of the original particles
    // by adding to r a vector in v direction
    //the angular momentum is still conserved!
    const Tgr3v rA = rCM + r;
    const Tgr3v rB = rCM - r;
    const Tgr3v vA = vCM + v;
    const Tgr3v vB = vCM - v;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    A.x = rA(0);
    A.y = rA(1);
    A.z = rA(2);
    B.x = rB(0);
    B.y = rB(1);
    B.z = rB(2);
    A.vx = vA(0);
    A.vy = vA(1);
    A.vz = vA(2);
    B.vx = vB(0);
    B.vy = vB(1);
    B.vz = vB(2);
#else
    for(int d=0; d<3; d++) {
        r1_new[d] = rA(d);
        r2_new[d] = rB(d);
        v1_new[d] = vA(d);
        v2_new[d] = vB(d);
    }
    // Here is very important the sequence of transformation
    sph_transf_C2S_r(r1_new);
    sph_transf_C2S_r(r2_new);
    sph_transf_S2H_R(r1_new);
    sph_transf_S2H_R(r2_new);
    A.x = r1_new[0];
    A.y = r1_new[1];
    A.z = r1_new[2];
    B.x = r2_new[0];
    B.y = r2_new[1];
    B.z = r2_new[2];
    A.vx = v1_new[0];
    A.vy = v1_new[1];
    A.vz = v1_new[2];
    B.vx = v2_new[0];
    B.vy = v2_new[1];
    B.vz = v2_new[2];
#endif
    A.w = B.w = w;
    A.popid = B.popid = P1.popid;
}

void Join::resetParameters()
{
    // JOIN PARAMETERS
    //none
}

// JOIN METHODS

//! Set function arguments
void Join::setArgs_joinDefault()
{
    if(args.size() != 0) {
        ERRORMSG2("function takes no arguments",name);
        doabort();
    }
}

/** \brief Combine particles
 *
 * Particles are combined conserving momentum, kinetic energy, angular
 * momentum and the position of the center of mass. Only particles in
 * the same population are combined.
 *
 * First, find out how many different particle populations have
 * at least MIN_COUNT(=6) macroparticles in this cell.
 * Second, choose a suitable population with prob: N/Ntot.
 * Third, find two particles in that population
 * which minimize (w1+w2)|v1 - v2|.
 * Forth, find third particle in that population
 * which minimizes w3(v3-v_2cm), v_2cm=0.5(v1+v2);
 * Finally, call the join32 function with these 3 particle pointers.
 *
 * if the optional input is 'false', use fast joining:
 * in Third step: find the smallest weight w1,
 *   then other that minizes the same (w1+w2)|v1-v2|
 * otherwise the same - makes time consumption scale as N and not N*N
 */
int Join::joinDefault(const gridreal boxmin[3], const gridreal boxmax[3], TParticleList& tplist, int fast)
{
    TLinkedParticle *p,*q,*r;
    //choose the population
    const int MIN_COUNT = 6;
    int poppi, kp, ncount, kk=0, chosen_pop=-1;
    int pop_choose[Params::POPULATIONS];
    int pop_count[Params::POPULATIONS], totcount=0;
    for (poppi=0; poppi<Params::POPULATIONS; poppi++) {
        pop_count[poppi] = 0;
    }
    for (p=tplist.first; p; p=p->next) {
        pop_count[p->popid]++;
    }
    //choose an available population p_i with probability N_p_i/N_tot
    for (poppi=0; poppi<Params::POPULATIONS; poppi++) {
        if (Params::pops[poppi]->getJoin() && pop_count[poppi] > MIN_COUNT) {
            pop_choose[kk] = poppi;
            totcount+=pop_count[poppi];
            kk++;
        }
    }
    if (kk==0) {
        return 0;
    }
    int rand_pop = int(uniformrnd() * totcount);
    // just to make sure in the (no idea how) unlikely event: uniformrnd()==1
    if (rand_pop==totcount) {
        rand_pop--;
        errorlog << "JoinDefault: uniformrnd() gave 1.00000" << endl;
    }
    for (kp=0, ncount=0; kp<kk; kp++) {
        ncount+= pop_count[pop_choose[kp]];
        if (rand_pop<ncount) {
            chosen_pop = pop_choose[kp];
            break;
        }
    }
    // Find the first two particles
    real merit=0, bestmerit=-1;
    TLinkedParticle *pp[4] = {NULL,NULL,NULL,NULL};
    if (fast==0) {
        for (p = tplist.first->next; p; p = p->next) {
            if (p->popid!=chosen_pop) {
                continue;
            }
            for (q = p->next; q; q = q->next) {
                if (q->popid!=chosen_pop) {
                    continue;
                }
                merit = real(p->w + q->w)*(p->w * q->w) *
                        (sqr(p->vx - q->vx) + sqr(p->vy - q->vy) + sqr(p->vz - q->vz));
                if (bestmerit < 0 || merit < bestmerit) {
                    bestmerit = merit;
                    pp[0] = p;
                    pp[1] = q;
                }
            }
        }
    } else { //fast join
        gridreal w_low=-1;
        for (p = tplist.first->next; p; p = p->next) {
            if (p->popid!=chosen_pop) {
                continue;
            }
            if (w_low<0||p->w<w_low) {
                w_low = p->w;
                pp[0] = p;
            }
        }
        for (q = tplist.first->next; q; q = q->next) {
            if (q->popid!=chosen_pop || pp[0]==q) {
                continue;
            }
            merit = real(pp[0]->w + q->w)*(pp[0]->w * q->w) *
                    (sqr(pp[0]->vx - q->vx) + sqr(pp[0]->vy - q->vy) + sqr(pp[0]->vz - q->vz));
            if (bestmerit < 0 || merit < bestmerit) {
                bestmerit = merit;
                pp[1] = q;
            }
        }
    }
    const fastreal v1[3]= {pp[0]->w * pp[0]->vx, pp[0]->w * pp[0]->vy, pp[0]->w * pp[0]->vz};
    const fastreal v2[3]= {pp[1]->w * pp[1]->vx, pp[1]->w * pp[1]->vy, pp[1]->w * pp[1]->vz};
    const fastreal invw = 1.0 / (pp[0]->w + pp[1]->w);
    const fastreal v_2cm[3] = {invw*(v1[0] + v2[0]),invw*(v1[1] + v2[1]),invw*(v1[2] + v2[2])};
    // Find the third particle
    bestmerit=-1;
    for (r = tplist.first->next; r; r = r->next) {
        if (r->popid!=chosen_pop || r==pp[0] || r==pp[1]) {
            continue;
        }
        merit = real(r->w) * r->w *
                (sqr(r->vx - v_2cm[0]) + sqr(r->vy - v_2cm[1]) + sqr(r->vz - v_2cm[2]));
        if (bestmerit < 0 || merit < bestmerit) {
            bestmerit = merit;
            pp[2] = r;
        }
    }
    // Must copy *p1,*p2 so that link field isn't destroyed
    TLinkedParticle *const p1=pp[0];
    TLinkedParticle *const p2=pp[1];
    TLinkedParticle *const p3=pp[2];
    const TLinkedParticle P1=*p1;
    const TLinkedParticle P2=*p2;
    const TLinkedParticle P3=*p3;
    TLinkedParticle PA = *pp[0], PB = *pp[1];
    join32(P1,P2,P3, PA,PB);
    if (!Params::insideBoxTight(&PA) || !Params::insideBoxTight(&PB))  {
        return 0;
    }
    // Situation ok, copy PA,PB to data structure
    // Link fields are inherited from original *p1,*p2
    *pp[0] = PA;
    *pp[1] = PB;
    q = pp[2];
    // Now delete p3 (and p4)
    if (pp[2] == tplist.first) {
        tplist.first = pp[2]->next;
    } else {
        // Find predecessor
        TLinkedParticle *pred;
        for (pred=tplist.first; pred->next!=q; pred=pred->next);
        pred->next = pp[2]->next;
    }
    delete q;
#ifndef NO_DIAGNOSTICS
    // Increase counter
    Params::diag.pCounter[PA.popid]->joiningRate += 1;
#endif
    tplist.n_part--;
    return -1;
}

//! Set function arguments
void Join::setArgs_joinOriginal()
{
    if(args.size() != 0) {
        ERRORMSG2("function takes no arguments",name);
        doabort();
    }
}

//! Original join function (coalesce)
int Join::joinOriginal(const gridreal minbox[3], const gridreal maxbox[3], TParticleList& tplist, int)
{
    TLinkedParticle *p,*q;
    // select three particles for coalescing
    // first find out how many different mass-species occur, all coalesced particles
    // must belong in the same species.
    int poppi,i,selected_species = -1;
    int pop_count[Params::POPULATIONS];
    for (poppi=0; poppi<Params::POPULATIONS; poppi++) {
        pop_count[poppi] = 0;
    }
    for (p=tplist.first; p; p=p->next) {
        pop_count[p->popid]++;
    }
    // now mass[0..nspecies-1] contains the mass spectrum of this plist
    // and N_of_mass[0..nspecies-1] the corresponding frequencies (numbers)
    //   - now pop_count[] replaces N_of_mass [OLD mass window stuff]
    TLinkedParticle *p1s[Params::POPULATIONS], *p2s[Params::POPULATIONS], *p3s[Params::POPULATIONS];
    real vCM[3], bestmerit, best_species_merit = -1;
    for (i=0; i<Params::POPULATIONS; i++) {
        if (pop_count[i] < 3) continue;
        // select the best pair from species i
        // the best pair has the smallest (w1+w2)*(v1-v2)^2
        bestmerit = -1;
        for (p=tplist.first; p; p=p->next) for (q=p->next; q; q=q->next) {
                if (p->popid != i || q->popid != i) continue;	// pass through only species i
                const real merit = real(p->w + q->w)*real(sqr(p->vx-q->vx) + sqr(p->vy-q->vy) + sqr(p->vz-q->vz));
                // CAVEAT! If shortreal=float, merit can overflow! Thus weight and v^2 must be casted to real explicitly!
                if (bestmerit < 0 || merit < bestmerit) {
                    bestmerit = merit;
                    p1s[i] = p;
                    p2s[i] = q;
                }
            }
        // now p1s[i],p2s[i] contains the best pair
        // find the third member so that w1*(v1-vCM)^2 + w2*(v2-vCM)^2 + w3*(v3-vCM)^2 is minimum
        // where vCM is the center of mass velocity of all three particles
        bestmerit = -1;
        for (p=tplist.first; p; p=p->next) {
            if (p->popid != i || p == p1s[i] || p == p2s[i]) continue;
            const fastreal invwsum = p1s[i]->w + p2s[i]->w + p->w;
            vCM[0] = (real(p1s[i]->w)*p1s[i]->vx + real(p2s[i]->w)*p2s[i]->vx + real(p->w)*p->vx)*invwsum;
            vCM[1] = (real(p1s[i]->w)*p1s[i]->vy + real(p2s[i]->w)*p2s[i]->vy + real(p->w)*p->vy)*invwsum;
            vCM[2] = (real(p1s[i]->w)*p1s[i]->vz + real(p2s[i]->w)*p2s[i]->vz + real(p->w)*p->vz)*invwsum;
            // CAVEAT! Also here it is safest to cast the multiplicants to real
            const real merit =
                real(p1s[i]->w)*real(sqr(p1s[i]->vx-vCM[0]) + sqr(p1s[i]->vy-vCM[1]) + sqr(p1s[i]->vz-vCM[2])) +
                real(p2s[i]->w)*real(sqr(p2s[i]->vx-vCM[0]) + sqr(p2s[i]->vy-vCM[1]) + sqr(p2s[i]->vz-vCM[2])) +
                real(p->w)*real(sqr(p->vx-vCM[0]) + sqr(p->vy-vCM[1]) + sqr(p->vz-vCM[2]));
            // CAVEAT! And here
            if (bestmerit < 0 || merit < bestmerit) {
                bestmerit = merit;
                p3s[i] = p;
            }
        }
        const real speciesmerit = Params::pops[i]->m*bestmerit;
        if (best_species_merit < 0 || speciesmerit < best_species_merit) {
            best_species_merit = speciesmerit;
            selected_species = i;
        }
    }
    TLinkedParticle *const p1 = p1s[selected_species];
    TLinkedParticle *const p2 = p2s[selected_species];
    TLinkedParticle *const p3 = p3s[selected_species];
    const TLinkedParticle P1 = *p1;
    const TLinkedParticle P2 = *p2;
    const TLinkedParticle P3 = *p3;
    TLinkedParticle PA=*p1,PB=*p2;		// must copy *p1,*p2 so that link field isn't destroyed
    join32(P1,P2,P3, PA,PB);
    if (!Params::insideBoxTight(&PA) || !Params::insideBoxTight(&PB))  {
        return 0;
    }
    // situation ok, copy PA,PB to data structure
    *p1 = PA;		// link fields are inherited from original *p1,*p2
    *p2 = PB;
    // now delete p3
    if (p3 == tplist.first) {
        q = p3;
        tplist.first = p3->next;
        delete q;
    } else {
        q = p3;
        // find predecessor
        TLinkedParticle *pred;
        for (pred=tplist.first; pred->next!=q; pred=pred->next);
        pred->next = p3->next;
        delete q;
    }
#ifndef NO_DIAGNOSTICS
    // Increase counter
    Params::diag.pCounter[PA.popid]->joiningRate += 1;
#endif
    tplist.n_part--;
    return -1;
}

//! Initialize macro particle split&join
void initializeSplitJoin()
{
    MSGFUNCTIONCALL("initializeSplitJoin");
    mainlog << "SPLIT FUNCTION: ";
    vector<string> tempNames;
    vector< vector<real> > tempArgs;
    bool nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("splitFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false) {
        ERRORMSG("splitting function has to be defined");
        doabort();
    } else if (tempNames.size() != 1 || tempArgs.size() != 1) {
        ERRORMSG("only one splitting function can be defined");
        doabort();
    } else {
        mainlog << tempNames[0] << " ";
        for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
            mainlog << tempArgs[0][i] << " ";
        }
        mainlog << "\n";
        Params::splittingFunction = Split(tempNames[0], tempArgs[0]);
    }
    mainlog << "JOIN FUNCTION: ";
    tempNames.clear();
    tempArgs.clear();
    nonZeroFuncs = simuConfig.getFunctionNamesAndArgs("joinFUNC",tempNames,tempArgs);
    if(nonZeroFuncs == false) {
        ERRORMSG("joining function has to be defined");
        doabort();
    } else if (tempNames.size() != 1 || tempArgs.size() != 1) {
        ERRORMSG("only one joining function can be defined");
        doabort();
    } else {
        mainlog << tempNames[0] << " ";
        for(unsigned int i=0; i < tempArgs[0].size(); ++i) {
            mainlog << tempArgs[0][i] << " ";
        }
        mainlog << "\n";
        Params::joiningFunction = Join(tempNames[0],tempArgs[0]);
    }
    MSGFUNCTIONEND("initializeSplitJoin");
}

