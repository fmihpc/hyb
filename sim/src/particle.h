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

#ifndef PLIST_H
#define PLIST_H

#include <iostream>
#include <vector>
#include "definitions.h"
#include <stdint.h>

//! Linked simulation macroparticle
struct TLinkedParticle {
    shortreal x,y,z,vx,vy,vz; //!< Position and velocity of a macro particle
    shortreal w; //!< Statistical weight of a macro particle (=how many real particle a macro particle represents)
    int popid; //!< Particle population ID
    TLinkedParticle *next; //!< Next particle in linked list
#ifdef USE_PARTICLE_SUBCYCLING
    uint8_t dtlevel;
    real accumed; //!< debug var
#endif
};

//! Arguments from the grid to the particle pass function
struct ParticlePassArgs {
    datareal rho_q;
    gridreal size;
};

//! Linked particle list (unidirectional), for storing to grid cells. For all functions having popID[], pops: take only particles in the specified populations.
class TParticleList
{
private:
    TLinkedParticle *first;
    int n_part;
    bool particleInPop(const TLinkedParticle& P, const std::vector<int> popId) const;
    friend class Split;
    friend class Join;
    template <class Ret, class Func> friend class ParticleMapper;
public:
    void init() {
        first = 0;
        n_part=0;
    }
    TParticleList() {
        init();
    }
    void add(shortreal x, shortreal y, shortreal z, shortreal vx, shortreal vy, shortreal vz, shortreal w, int popid);
    template <class Func> int pass(Func& op);
    template <class Func> void pass(Func& op) const;
    int pass(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), ParticlePassArgs a);
    template <class Func> int pass_with_relocate(Func& op);
    int pass_with_relocate(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), ParticlePassArgs a);
    int Nparticles() const;
    real calc_weight(std::vector<int> popId = std::vector<int>()) const;
    real calc_mass(std::vector<int> popId = std::vector<int>()) const;
    real calc_charge(std::vector<int> popId = std::vector<int>()) const;
    void calc_avev(real& vx0, real& vy0, real& vz0, std::vector<int> popId = std::vector<int>()) const;
    void calc_U(real& Ux0, real& Uy0, real& Uz0, std::vector<int> popId = std::vector<int>()) const;
    real calc_avemv2(real vx0, real vy0, real vz0, std::vector<int> popId = std::vector<int>()) const;
    friend std::ostream& operator<<(std::ostream& o, const TParticleList& pl);
    std::string toString() const;
    ~TParticleList();
};

//! Return the number of particles in a particle list
inline int TParticleList::Nparticles() const
{
    return n_part;
}

#endif

