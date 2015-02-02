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
#include "particle.h"
#include "random.h"
#include "simulation.h"
#include "params.h"

using namespace std;

extern Tgrid g;

//! Add one new particle with given parameters.
void TParticleList::add(shortreal x, shortreal y, shortreal z, shortreal vx, shortreal vy, shortreal vz, shortreal w, int popid)
{
    TLinkedParticle *const p = new TLinkedParticle;
    p->x = x;
    p->y = y;
    p->z = z;
    p->vx = vx;
    p->vy = vy;
    p->vz = vz;
    p->w = w;
    p->popid = popid;
    p->next = first;
    first = p;
    n_part++;
#ifdef USE_PARTICLE_SUBCYCLING
    p->dtlevel = 0;
    p->accumed = 1;
#endif
}

/** \brief Call op for all particles
 *
 * If op returns false, delete the particle afterwards. Returns
 * number of deletions.
 */
int TParticleList::pass(bool (*op)(TLinkedParticle& part,ParticlePassArgs a), ParticlePassArgs a)
{
    TLinkedParticle *p,*prev,*q;
    int ndel = 0;
    for (p=first,prev=0; p;) {
        if ((*op)(*p,a)) {
            prev = p;
            p = p->next;
        } else {
            q = p;
            if (prev) {
                prev->next = p->next;
            } else {
                first = p->next;
            }
            p = p->next;
            delete q;
            ndel++;
            n_part--;
        }
    }
    return ndel;
}

//! Pass thru the particles in the list with relocation
int TParticleList::pass_with_relocate(bool (*op)(TLinkedParticle& p,ParticlePassArgs a), ParticlePassArgs a)
{
    TLinkedParticle *p,*prev,*q;
    int ndel = 0;
    for (p=first,prev=0; p;) {
        if ((*op)(*p,a)) {
            TParticleList *newplist = g.find_plist(*p);
            if (newplist != NULL && newplist != this) {
                // particle p needs to be moved from *this to *newplist
                q = p;
                if (prev) prev->next = p->next;
                else first = p->next;
                p = p->next;
                q->next = newplist->first;
                newplist->first = q;
                n_part--;
                newplist->n_part++;
            }
            // Remove the particle if no particle list found (=out of box)
            else if(newplist == NULL) {
                ERRORMSG("no particle list found, removing particle");
                q = p;
                if (prev) prev->next = p->next;
                else first = p->next;
                p = p->next;
                delete q;
                ndel++;
                n_part--;
            } else {
                prev = p;
                p = p->next;
            }
        } else {
            q = p;
            if (prev) prev->next = p->next;
            else first = p->next;
            p = p->next;
            delete q;
            ndel++;
            n_part--;
        }
    }
    return ndel;
}

//! Check if the particle belongs into any of the populations in popId. If popId.size() <= 0, return true.
bool TParticleList::particleInPop(const TLinkedParticle& P, const vector<int> popId) const
{
    unsigned int popsN = popId.size();
    if (popsN <= 0) {
        return true;
    } else {
        for(unsigned int i = 0; i < popsN; i++) {
            if (P.popid == popId[i]) {
                return true;
            }
        }
        return false;
    }
}

//! Calculate sum(w) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
real TParticleList::calc_weight(vector<int> popId) const
{
    real result = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return 0;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            result += static_cast<real>(p->w);
        }
    }
    return result;
}

//! Calculate sum(w*m) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
real TParticleList::calc_mass(vector<int> popId) const
{
    real result = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return 0;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            result += static_cast<real>(p->w)*real(Params::pops[p->popid]->m);
        }
    }
    return result;
}

//! Calculate sum(w*m) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
real TParticleList::calc_charge(vector<int> popId) const
{
    real result = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return 0;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            result += static_cast<real>(p->w)*real(Params::pops[p->popid]->q);
        }
    }
    return result;
}

//! Calculate sum(w*v)/sum(w) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
void TParticleList::calc_avev(real& vx0, real& vy0, real& vz0, vector<int> popId) const
{
    vx0 = vy0 = vz0 = 0;
    real wsum = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            vx0 += static_cast<real>(p->w)*real(p->vx);
            vy0 += static_cast<real>(p->w)*real(p->vy);
            vz0 += static_cast<real>(p->w)*real(p->vz);
            wsum += static_cast<real>(p->w);
        }
    }
    if (wsum == 0) {
        return;
    }
    const real invwsum = 1.0/wsum;
    vx0 *= invwsum;
    vy0 *= invwsum;
    vz0 *= invwsum;
}

//! Calculate U = sum(m*w*v)/sum(m*w) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
void TParticleList::calc_U(real& Ux0, real& Uy0, real& Uz0, vector<int> popId) const
{
    Ux0 = Uy0 = Uz0 = 0;
    real wsum = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            const real wfac = static_cast<real>(p->w) * static_cast<real>(Params::pops[p->popid]->m);
            Ux0 += wfac*static_cast<real>(p->vx);
            Uy0 += wfac*static_cast<real>(p->vy);
            Uz0 += wfac*static_cast<real>(p->vz);
            wsum += wfac;
        }
    }
    if (wsum == 0) {
        return;
    }
    const real invwsum = 1.0/wsum;
    Ux0 *= invwsum;
    Uy0 *= invwsum;
    Uz0 *= invwsum;
}

//! Calculate sum(w*m*(v-v0)^2)/sum(w) over pops listed in popId. If popId.size() <= 0, calculate over all populations.
real TParticleList::calc_avemv2(real vx0, real vy0, real vz0, vector<int> popId) const
{
    real mv2 = 0, denom = 0;
    TLinkedParticle *p;
    if (first == 0) {
        return 0;
    }
    for (p=first; p; p=p->next) {
        if (particleInPop(*p, popId) == true) {
            mv2 += real(p->w)*static_cast<real>(Params::pops[p->popid]->m)
                   * (sqr(p->vx-vx0) + sqr(p->vy-vy0) + sqr(p->vz-vz0));
            denom += p->w;
        }
    }
    if (denom == 0) {
        return 0;
    }
    return mv2/denom;
}

//! Stream operator
ostream& operator<<(ostream& o, const TParticleList& pl)
{
    TLinkedParticle *p;
    o << '(';
    for (p=pl.first; p; p=p->next) {
        o << "x=" << p->x << ",y=" << p->y << ",z=" << p->z
          << ",vx=" << p->vx << ",vy=" << p->vy << ",vz=" << p->vz << "; ";
    }
    o << ')';
    o << ": mv2=" << pl.calc_avemv2(0, 0, 0) << ",npart=" << pl.n_part;
    return o;
}

//! Particle list to string
string TParticleList::toString() const
{
    stringstream ss;
    ss << "(";
    TLinkedParticle *p;
    for (p = first; p; p = p->next) {
        ss << "x=" << p->x << ",y=" << p->y << ",z=" << p->z
           << ",vx=" << p->vx << ",vy=" << p->vy << ",vz=" << p->vz << "; ";
    }
    ss << "): mv2=" << calc_avemv2(0,0,0) << ", npart=" << n_part;
    return ss.str();
}

//! Destructor
TParticleList::~TParticleList()
{
    TLinkedParticle *p=first,*q;
    while (p) {
        q = p->next;
        delete p;
        p = q;
    }
}

