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

#ifndef TEMPLATES_H
#define TEMPLATES_H

extern Tgrid g;
extern Logger mainlog, errorlog, paramslog;

#define ForAll(i,j,k) for (i=0; i<nx; i++) for (j=0; j<ny; j++) for (k=0; k<nz; k++)
#define ForAllxyz(i,j,k) for (k=0; k<nz; k++) for (j=0; j<ny; j++) for (i=0; i<nx; i++)
#define ForInterior(i,j,k) for (i=1; i<nx-1; i++) for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)

//! Returns the number of particles deleted
template <class Func>
int Tgrid::Tcell::particle_pass_recursive(Func& op, bool relocate)
{
    int ndel = 0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) ndel+= child[0][0][ch]->particle_pass_recursive(op,relocate);
    } else {
        if (relocate)
            ndel = plist.pass_with_relocate(op);
        else
            ndel = plist.pass(op);
    }
    return ndel;
}

//! Pass all particles in the list to the function op
template <class Func>
int Tgrid::particle_pass(Func op, bool relocate)
{
    int i,j,k,ndel=0;
    TCellPtr c;
    ForAll(i,j,k) {
        c = cells[flatindex(i,j,k)];
        ndel+= c->particle_pass_recursive(op,relocate);
    }
    n_particles-= ndel;
    return ndel;
}

//! Pass cells (recursive)
template <class Func>
void Tgrid::Tcell::cellPassRecursive(Func& op)
{
    if (haschildren) {
        const int childCount = 8;
        for (int ch = 0; ch < childCount; ++ch)
            child[0][0][ch]->cellPassRecursive(op);
    } else {
        op(*this);
    }
}

//! Pass cells
template <class Func>
void Tgrid::cellPass(Func op)
{
    int i, j, k;
    ForAllxyz(i, j, k)
    cells[flatindex(i, j, k)]->cellPassRecursive(op);
}

/** \brief Call op for all particles
 *
 * If op returns false, delete the particle afterwards. Returns
 * number of deletions.
 */
template <class Func>
int TParticleList::pass(Func& op)
{
    TLinkedParticle *p,*prev,*q;
    int ndel = 0;
    for (p=first,prev=0; p;) {
        if (op(*p)) {
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

/** \brief Call op for all particles
 *
 * Const pass. Like pass, but doesn't change anything.
 */
template <class Func>
void TParticleList::pass(Func& op) const
{
    for (TLinkedParticle* p=first; p; p=p->next) {
        op(*p);
    }
}

//! Pass thru the particles in the list with relocation
template <class Func>
int TParticleList::pass_with_relocate(Func& op)
{
    TLinkedParticle *p,*prev,*q;
    int ndel = 0;
    for (p=first,prev=0; p;) {
        if (op(*p)) {
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

#endif

