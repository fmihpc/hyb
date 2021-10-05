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

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include "grid.h"
#include "magneticfield.h"
#include "random.h"
#include "simulation.h"
#include "templates.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

#define ForAll(i,j,k) for (i=0; i<nx; i++) for (j=0; j<ny; j++) for (k=0; k<nz; k++)
#define ForInterior(i,j,k) for (i=1; i<nx-1; i++) for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)

const char *Tgrid::celldata_names[Tgrid::NCELLDATA] = {"u","ue","j","B"};
int Tgrid::cell_running_index = 0;
Tgrid::TPtrHash *Tgrid::hp = 0;
FieldCounter Tgrid::fieldCounter;

static bool hcFileAsciiFormat = false;

//! Constructor
Tgrid::Tgrid()
{
    MSGFUNCTIONCALL("Tgrid::Tgrid");
    MSGFUNCTIONEND("Tgrid::Tgrid");
}

//! Destructor
Tgrid::~Tgrid()
{
    MSGFUNCTIONCALL("Tgrid::~Tgrid");
    MSGFUNCTIONEND("Tgrid::~Tgrid");
}

/** \brief Grid initialization
 *
 * Basegrid size NOT including ghost cells
 * Basegrid (NOT including ghosts) lowest corner point
 * Grid spacing (isotropic)
 */
void Tgrid::init(int nx1, int ny1, int nz1, gridreal x1, gridreal y1, gridreal z1, gridreal bgdx1)
{
    MSGFUNCTIONCALL("Tgrid::init");
    static bool onlyOneTgridObject = false;
    if (onlyOneTgridObject == true) {
        ERRORMSG("only one Tgrid object can be defined");
        doabort();
    } else {
        onlyOneTgridObject = true;
    }
    nx = nx1 + 2;
    ny = ny1 + 2;
    nz = nz1 + 2;
    bgdx = bgdx1;
    saved_cellptr = 0;
    n_particles = 0;
    n_pdftables = 0;
    ave_ntimes = 0;
    previous_found_cell = 0;
    invbgdx = 1.0/bgdx;
    x_1 = x1 - bgdx;
    y_1 = y1 - bgdx;
    z_1 = z1 - bgdx;
    // "unit" is the smallest safely representable distance using gridreal
    inv_unit = invbgdx*65536.0;
    const int N = nx*ny*nz;
    mainlog << "|------------------------- GRID DETAILS -------------------------|\n";
    mainlog << "| Grid [nx,ny,nz]  = [" << nx << "," << ny << "," << nz << "] (ghosts included)\n"
            << "| Base cells       = " << N << " (with ghosts)\n"
            << "| Base cells       = " << nx1*ny1*nz1 << " (no ghosts)\n"
            << "| Base cell size   = " << bgdx/1e3 << " km = " << bgdx/Params::R_P << " R_P\n"
            << "| Box size [x,y,z] = [" << Params::box_X/1e3 << ", " << Params::box_Y/1e3 << ", " << Params::box_Z/1e3 << "] km  = [" << Params::box_X/Params::R_P << "," << Params::box_Y/Params::R_P << "," << Params::box_Z/Params::R_P << "] R_P\n"
            << "| Box volume       = " << Params::box_V << " m^3 = " << Params::box_V/cube(Params::R_P) << " R_P^3\n"
            << "| Gridreal unit    = " << 1.0/inv_unit << " m\n"
            << "| Bytes/cell       = " << approx_bytes_per_cell() << " (approx)\n";
    mainlog << "|----------------------------------------------------------------|\n\n";
    mainlog << "|--------------- COORDINATE LIMITS ---------------|\n";
    mainlog << "| " << Params::box_xmin/1e3 << " km (" << Params::box_xmin/Params::R_P << " R_P) < x < " << Params::box_xmax/1e3 << " km (" << Params::box_xmax/Params::R_P << " R_P)\n"
            << "| " << Params::box_ymin/1e3 << " km (" << Params::box_ymin/Params::R_P << " R_P) < y < " << Params::box_ymax/1e3 << " km (" << Params::box_ymax/Params::R_P << " R_P)\n"
            << "| " << Params::box_zmin/1e3 << " km (" << Params::box_zmin/Params::R_P << " R_P) < z < " << Params::box_zmax/1e3 << " km (" << Params::box_zmax/Params::R_P << " R_P)\n";
    mainlog << "|-------------------------------------------------|\n";
    cells = new TCellPtr [N];
    int i,j,k,c;
    // Create cells and faces, set up cell fields but not yet face fields
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c] = new Tcell;
        cells[c]->plist.init();
        cells[c]->haschildren = false;
        cells[c]->refine_it = cells[c]->recoarsen_it = cells[c]->forbid_psplit = false;
        cells[c]->refstatus = 0;
        cells[c]->level = 0;
        cells[c]->parent = 0;
        cells[c]->running_index = -123456;
        cells[c]->face[0][0] = (i > 0) ? cells[flatindex(i-1,j,k)]->face[0][1] : TFacePtr(0);
        cells[c]->face[0][1] = (i < nx-1) ? new Tface : 0;
        cells[c]->face[1][0] = (j > 0) ? cells[flatindex(i,j-1,k)]->face[1][1] : TFacePtr(0);
        cells[c]->face[1][1] = (j < ny-1) ? new Tface : 0;
        cells[c]->face[2][0] = (k > 0) ? cells[flatindex(i,j,k-1)]->face[2][1] : TFacePtr(0);
        cells[c]->face[2][1] = (k < nz-1) ? new Tface : 0;
        cells[c]->nc = 0.0;
        cells[c]->rho_q = 0.0;
        cells[c]->rho_q_bg = 0.0;
        cells[c]->ave_nc = 0.0;
#ifdef SAVE_PARTICLES_ALONG_ORBIT
        cells[c]->save_particles = false;
#endif
#ifdef SAVE_POPULATION_AVERAGES
        for (int l = 0; l < Params::POPULATIONS; ++l) {
            cells[c]->pop_ave_n.push_back(0.0);
            cells[c]->pop_ave_vx.push_back(0.0);
            cells[c]->pop_ave_vy.push_back(0.0);
            cells[c]->pop_ave_vz.push_back(0.0);
        }
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
        for (int l = 0; l < Params::POPULATIONS; ++l) {
            vector<datareal> temp;
            cells[c]->spectra.push_back(temp);
            for (int m = 0; m < Params::spectraNbins; ++m) {
                cells[c]->spectra[l].push_back(0.0);
            }
        }
#endif
        memset(&cells[c]->celldata[0][0],0,sizeof(datareal)*NCELLDATA*3);
        cells[c]->centroid[0] = x_1 + (i+0.5)*bgdx;
        cells[c]->centroid[1] = y_1 + (j+0.5)*bgdx;
        cells[c]->centroid[2] = z_1 + (k+0.5)*bgdx;
        cells[c]->size = bgdx;
        cells[c]->r2 = vecsqr(cells[c]->centroid);
        cells[c]->invsize = 1.0/bgdx;
        cells[c]->flatind = c;
    }
    // Set up neighbour fields
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->neighbour[0][0] = (i > 0)    ? cells[flatindex(i-1,j,k)] : TCellPtr(0);
        cells[c]->neighbour[0][1] = (i < nx-1) ? cells[flatindex(i+1,j,k)] : TCellPtr(0);
        cells[c]->neighbour[1][0] = (j > 0)    ? cells[flatindex(i,j-1,k)] : TCellPtr(0);
        cells[c]->neighbour[1][1] = (j < ny-1) ? cells[flatindex(i,j+1,k)] : TCellPtr(0);
        cells[c]->neighbour[2][0] = (k > 0)    ? cells[flatindex(i,j,k-1)] : TCellPtr(0);
        cells[c]->neighbour[2][1] = (k < nz-1) ? cells[flatindex(i,j,k+1)] : TCellPtr(0);
    }
    // Create nodes, set up node fields.
    // Must store the nodes somewhere, use temporary vector 'nodes' for it.
    TNodePtr *nodes = new TNodePtr [N];
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                const TNodePtr n = new Tnode;
                n->cell[0][0][0] = cells[c];
                n->cell[1][0][0] = cells[flatindex(i+1,j,  k  )];
                n->cell[0][1][0] = cells[flatindex(i,  j+1,k  )];
                n->cell[1][1][0] = cells[flatindex(i+1,j+1,k  )];
                n->cell[0][0][1] = cells[flatindex(i,  j,  k+1)];
                n->cell[1][0][1] = cells[flatindex(i+1,j,  k+1)];
                n->cell[0][1][1] = cells[flatindex(i,  j+1,k+1)];
                n->cell[1][1][1] = cells[flatindex(i+1,j+1,k+1)];
                n->centroid[0] = x_1 + (i+1)*bgdx;
                n->centroid[1] = y_1 + (j+1)*bgdx;
                n->centroid[2] = z_1 + (k+1)*bgdx;
                n->r2 = vecsqr(n->centroid);
                nodes[c] = n;
            }
    // Set up pointers from faces to nodes. Use right-pointing faces (dir=1)
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                // X-directed faces:
                if (j > 0 && k > 0) cells[c]->face[0][1]->node[0] = nodes[flatindex(i,j-1,k-1)];
                if (k > 0         ) cells[c]->face[0][1]->node[1] = nodes[flatindex(i,j,  k-1)];
                if (true          ) cells[c]->face[0][1]->node[2] = nodes[flatindex(i,j,  k)];
                if (j > 0         ) cells[c]->face[0][1]->node[3] = nodes[flatindex(i,j-1,k)];

                // Y-directed faces:
                if (i > 0 && k > 0) cells[c]->face[1][1]->node[0] = nodes[flatindex(i-1,j,k-1)];
                if (i > 0         ) cells[c]->face[1][1]->node[1] = nodes[flatindex(i-1,j,k  )];
                if (true          ) cells[c]->face[1][1]->node[2] = nodes[flatindex(i,  j,k)];
                if (k > 0         ) cells[c]->face[1][1]->node[3] = nodes[flatindex(i  ,j,k-1)];

                // Z-directed faces:
                if (i > 0 && j > 0) cells[c]->face[2][1]->node[0] = nodes[flatindex(i-1,j-1,k)];
                if (j > 0         ) cells[c]->face[2][1]->node[1] = nodes[flatindex(i,  j-1,k)];
                if (true          ) cells[c]->face[2][1]->node[2] = nodes[flatindex(i,  j,  k)];
                if (i > 0         ) cells[c]->face[2][1]->node[3] = nodes[flatindex(i-1,j,  k)];
            }
    delete [] nodes;
    MSGFUNCTIONEND("Tgrid::init");
}

#ifndef USE_SPHERICAL_COORDINATE_SYSTEM

//! Find a cell at r and return a pointer to the cell or 0 (NULL) if no cell found
Tgrid::TCellPtr Tgrid::findcell(const shortreal r[3], gridreal* lowercorner)
{
    if (previous_found_cell) {
        bool isinside = true;
        int d;
        const gridreal halfsize = 0.5*previous_found_cell->size;
        for (d=0; d<3; d++)
            if (r[d] < previous_found_cell->centroid[d] - halfsize ||
                r[d] > previous_found_cell->centroid[d] + halfsize) {
                isinside = false;
                break;
            }
        if (isinside) {
            if (lowercorner) {
                lowercorner[0] = previous_found_cell->centroid[0] - halfsize;
                lowercorner[1] = previous_found_cell->centroid[1] - halfsize;
                lowercorner[2] = previous_found_cell->centroid[2] - halfsize;
            }
            return previous_found_cell;
        }
    }
    const int i = int((r[0]-x_1)*invbgdx);
    const int j = int((r[1]-y_1)*invbgdx);
    const int k = int((r[2]-z_1)*invbgdx);
    if (i < 1 || i > nx-2 || j < 1 || j > ny-2 || k < 1 || k > nz-2) {
        return 0;
    }
    const int ijk = flatindex(i,j,k);
    if (lowercorner) {
        lowercorner[0] = x_1 + i*bgdx;
        lowercorner[1] = y_1 + j*bgdx;
        lowercorner[2] = z_1 + k*bgdx;
    }
    TCellPtr c = cells[ijk];
    while (c->haschildren) {
        const int chx = (r[0] > c->centroid[0]);
        const int chy = (r[1] > c->centroid[1]);
        const int chz = (r[2] > c->centroid[2]);
        if (lowercorner) {
            register const gridreal halfsize = 0.5*c->size;
            if (chx) lowercorner[0]+= halfsize;
            if (chy) lowercorner[1]+= halfsize;
            if (chz) lowercorner[2]+= halfsize;
        }
        c = c->child[chx][chy][chz];
    }
    previous_found_cell = c;
    return c;
}

#else

//! (SPHERICAL) Spherical version findcell
Tgrid::TCellPtr Tgrid::findcell(const shortreal r[3], gridreal* lowercorner)
{
    if (previous_found_cell) {
        bool isinside = true;
        int d;
        const gridreal halfsize[3] = {0.5*previous_found_cell->size, 0.5*previous_found_cell->sph_sizey, 0.5*previous_found_cell->sph_sizez};

        for (d=0; d<3; d++)
            if (r[d] < previous_found_cell->centroid[d] - halfsize[d] || r[d] > previous_found_cell->centroid[d] + halfsize[d]) {
                isinside = false;
                break;
            }
        if (isinside) {
            if (lowercorner) {
                lowercorner[0] = previous_found_cell->centroid[0] - halfsize[0];
                lowercorner[1] = previous_found_cell->centroid[1] - halfsize[1];
                lowercorner[2] = previous_found_cell->centroid[2] - halfsize[2];
            }
            return previous_found_cell;
        }
    }
    const int i = int((r[0]-x_1)*invbgdx);
    const int j = int((r[1]-y_1)*sph_invbgdy);
    const int k = int((r[2]-z_1)*sph_invbgdz);
    if (i < 1 || i > nx-2 || j < 1 || j > ny-2 || k < 1 || k > nz-2) {
        return 0;
    }
    const int ijk = flatindex(i,j,k);
    if (lowercorner) {
        lowercorner[0] = x_1 + i*bgdx;
        lowercorner[1] = y_1 + j*sph_bgdy;
        lowercorner[2] = z_1 + k*sph_bgdz;
    }
    TCellPtr c = cells[ijk];
    /*while (c->haschildren) {
    	const int chx = (r[0] > c->centroid[0]);
    	const int chy = (r[1] > c->centroid[1]);
    	const int chz = (r[2] > c->centroid[2]);
    	if (lowercorner) {
    		register const gridreal halfsize[3] = {0.5*c->size, 0.5*c->sizey, 0.5*c->sizez};
    		if (chx) lowercorner[0]+= halfsize[0];
    		if (chy) lowercorner[1]+= halfsize[1];
    		if (chz) lowercorner[2]+= halfsize[2];
    	}
    c = c->child[chx][chy][chz];
    }*/
    previous_found_cell = c;
    return c;
}

#endif

//! Find cell to a maximum level of maxlevel
Tgrid::TCellPtr Tgrid::findcell_to_maxlevel(const shortreal r[3], int maxlevel) const
{
    const int i = int((r[0]-x_1)*invbgdx);
    const int j = int((r[1]-y_1)*invbgdx);
    const int k = int((r[2]-z_1)*invbgdx);
    if (i < 0 || i > nx-1 || j < 0 || j > ny-1 || k < 0 || k > nz-1) return 0;
    const int ijk = flatindex(i,j,k);
    TCellPtr c = cells[ijk];
    while (c->haschildren && c->level < maxlevel) {
        const int chx = (r[0] > c->centroid[0]);
        const int chy = (r[1] > c->centroid[1]);
        const int chz = (r[2] > c->centroid[2]);
        c = c->child[chx][chy][chz];
    }
    return c;
}

//! Face average
real Tgrid::Tcell::faceave(int dim, int dir, TFaceDataSelect s) const
{
    if (neighbour[dim][dir] == 0) return 0;
    if (isrefined_face(dim,dir)) {
        const Trefintf *const r = refintf[dim][dir];
        return 0.25*(
                   r->face[0]->facedata[s] +
                   r->face[1]->facedata[s] +
                   r->face[2]->facedata[s] +
                   r->face[3]->facedata[s]);
    } else {
        return face[dim][dir]->facedata[s];
    }
}

//! Cell average (recursive)
void Tgrid::Tcell::childave(TCellDataSelect cs, real result[3]) const
{
    int d,a;
    if (haschildren) {
        result[0] = result[1] = result[2] = 0;
        real result1[3];
        for (a=0; a<8; a++) {
            child[0][0][a]->childave(cs,result1);
            for (d=0; d<3; d++) result[d]+= result1[d];
        }
        for (d=0; d<3; d++) result[d]*= 0.125;
    } else {
        for (d=0; d<3; d++) result[d] = celldata[cs][d];
    }
}

//! Cell average of rho_q (recursive)
real Tgrid::Tcell::childave_rhoq() const
{
    int a;
    real result = 0;
    if (haschildren) {
        for (a=0; a<8; a++) result+= child[0][0][a]->childave_rhoq();
        result*= 0.125;
    } else {
        result = rho_q;
    }
    return result;
}

//! Cell average of density (recursive)
real Tgrid::Tcell::childave_nc() const
{
    int a;
    real result = 0;
    if (haschildren) {
        for (a=0; a<8; a++) result+= child[0][0][a]->childave_nc();
        result*= 0.125;
    } else {
        result = nc;
    }
    return result;
}

//! face2cell interpolation (recursive)
void Tgrid::Tcell::FC_recursive(TFaceDataSelect fs, TCellDataSelect cs)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->FC_recursive(fs,cs);
    } else {
        int d;
        for (d=0; d<3; d++)
            celldata[cs][d] = 0.5*(faceave(d,0,fs) + faceave(d,1,fs));
    }
}

//! node2cell interpolation (recursive)
void Tgrid::Tcell::NC_recursive(TNodeDataSelect ns,TCellDataSelect cs)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->NC_recursive(ns,cs);
    } else {
        int dir,d,f,f2;
        datareal tempx,tempy,tempz;
        celldata[cs][0]=0.;
        celldata[cs][1]=0.;
        celldata[cs][2]=0.;
        for(dir=0; dir<3; dir++)for(d=0; d<2; d++) {
                tempx=tempy=tempz=0.;
                if (isrefined_face(dir,d)) {
                    for (f=0; f<4; f++) for (f2=0; f2<4; f2++) {
                            tempx+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][0];
                            tempy+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][1];
                            tempz+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][2];
                        }
                    tempx*=0.25;
                    tempy*=0.25;
                    tempz*=0.25;
                } else {
                    for(f=0; f<4; f++) {
                        tempx+=face[dir][d]->node[f]->nodedata[ns][0];
                        tempy+=face[dir][d]->node[f]->nodedata[ns][1];
                        tempz+=face[dir][d]->node[f]->nodedata[ns][2];
                    }
                }
                celldata[cs][0]+=tempx/24.;
                celldata[cs][1]+=tempy/24.;
                celldata[cs][2]+=tempz/24.;
            }
    }
}

//! node2cell interpolation for smoothing (recursive)
void Tgrid::Tcell::NC_smoothing_recursive()
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->NC_smoothing_recursive();
    } else {
        int dir,d,f,f2;
        datareal tempnc,temprhoq,tempvx,tempvy,tempvz;
        nc=0.;
        rho_q=0.;
        celldata[CELLDATA_Ji][0]=0.;
        celldata[CELLDATA_Ji][1]=0.;
        celldata[CELLDATA_Ji][2]=0.;
        for(dir=0; dir<3; dir++)for(d=0; d<2; d++) {
                tempnc=temprhoq=tempvx=tempvy=tempvz=0.;
                if (isrefined_face(dir,d)) {
                    for (f=0; f<4; f++) for (f2=0; f2<4; f2++) {
                            tempnc+=refintf[dir][d]->face[f]->node[f2]->nn;
                            temprhoq+=refintf[dir][d]->face[f]->node[f2]->nodedata[NODEDATA_UE][0];
                            tempvx+=refintf[dir][d]->face[f]->node[f2]->nodedata[NODEDATA_J][0];
                            tempvy+=refintf[dir][d]->face[f]->node[f2]->nodedata[NODEDATA_J][1];
                            tempvz+=refintf[dir][d]->face[f]->node[f2]->nodedata[NODEDATA_J][2];
                        }
                    tempnc*=0.25;
                    temprhoq*=0.25;
                    tempvx*=0.25;
                    tempvy*=0.25;
                    tempvz*=0.25;
                } else {
                    for(f=0; f<4; f++) {
                        tempnc+=face[dir][d]->node[f]->nn;
                        temprhoq+=face[dir][d]->node[f]->nodedata[NODEDATA_UE][0];
                        tempvx+=face[dir][d]->node[f]->nodedata[NODEDATA_J][0];
                        tempvy+=face[dir][d]->node[f]->nodedata[NODEDATA_J][1];
                        tempvz+=face[dir][d]->node[f]->nodedata[NODEDATA_J][2];
                    }
                }
                nc+=tempnc/24.;
                rho_q+=temprhoq/24.;
                celldata[CELLDATA_Ji][0]+=tempvx/24.;
                celldata[CELLDATA_Ji][1]+=tempvy/24.;
                celldata[CELLDATA_Ji][2]+=tempvz/24.;
            }
    }
}

//! cell2node interpolation
void Tgrid::Tnode::CN1(TCellDataSelect cs, TNodeDataSelect ns)
{
    int a;
    real result[3] = {0,0,0};
    gridreal weightsum = 0.0;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        register const gridreal weight = c->invsize;
        result[0]+= weight*c->celldata[cs][0];
        result[1]+= weight*c->celldata[cs][1];
        result[2]+= weight*c->celldata[cs][2];
        weightsum+= weight;
    }
    register const real invweightsum = 1.0/weightsum;
    nodedata[ns][0] = invweightsum*result[0];
    nodedata[ns][1] = invweightsum*result[1];
    nodedata[ns][2] = invweightsum*result[2];
}

//! cell2node interpolation of density
void Tgrid::Tnode::CN1_density()
{
    int a;
    real result = 0.0;
    gridreal weightsum = 0.0;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        register const gridreal weight = c->invsize;
        result+= weight*c->nc;
        weightsum+= weight;
    }
    nn = result/weightsum;
}

//! cell2node interpolation for smoothing
void Tgrid::Tnode::CN1_smoothing()
{
    int a;
    real result = 0.0;
    gridreal weightsum = 0.0;
    gridreal invweightsum;
    register Tcell *c;
    real chargedensity,currentdensity0,currentdensity1,currentdensity2;
    chargedensity=currentdensity0=currentdensity1=currentdensity2=0.0;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;

        register const gridreal weight = c->invsize;
        result+= weight*c->nc;//number density
        chargedensity+=weight*c->rho_q;//charge density
        currentdensity0+=weight*c->celldata[CELLDATA_Ji][0];//ion current density: x component
        currentdensity1+=weight*c->celldata[CELLDATA_Ji][1];//ion current density: y component
        currentdensity2+=weight*c->celldata[CELLDATA_Ji][2];//ion current density: z component
        weightsum+= weight;
    }
    if(weightsum==0.)invweightsum=0;
    else invweightsum=1./weightsum;
    nn = result*invweightsum;
    nodedata[NODEDATA_UE][0]=chargedensity*invweightsum;//NODEDATA_UE is used as a temporary storage place for charge density on node
    nodedata[NODEDATA_J][0]=currentdensity0*invweightsum;//NODEDATA_J is used as a temporary storage place for ion current density on node
    nodedata[NODEDATA_J][1]=currentdensity1*invweightsum;//NODEDATA_J is used as a temporary storage place for ion current density on node
    nodedata[NODEDATA_J][2]=currentdensity2*invweightsum;//NODEDATA_J is used as a temporary storage place for ion current density on node
}

//! cell2node interpolation of rho_q
void Tgrid::Tnode::CN1_rhoq()
{
    int a;
    gridreal weightsum = 0.0;
    gridreal invweightsum;
    register Tcell *c;
    real chargedensity=0.0;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;

        register const gridreal weight = c->invsize;
        chargedensity+=weight*c->rho_q;
        weightsum+= weight;
    }
    if(weightsum==0.)invweightsum=0;
    else invweightsum=1./weightsum;
    nodedata[NODEDATA_UE][0]=chargedensity*invweightsum;//NODEDATA_UE is used as a temporary storage place for charge density on node
}

//! Upwind nodedata by using cell2node interpolation
void Tgrid::Tnode::CN_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    int d,ch;
    real result[3] = {0,0,0};
    register real weightsum = 0.0;
    gridreal displacement[3];
    for(d=0; d<3; d++) {
        displacement[d] = nodedata[uns][d];
    }
    Tcell *c;
    // Find the smallest cell size among the neighbours
    gridreal basic_dx = -1;
    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!ch) continue;
        if (c->size < basic_dx || basic_dx < 0) basic_dx = c->size;
    }
    basic_dx*= 0.5;
    // let basic_dx be one half of the smallest touching cell size
    real norm = sqrt(sqr(displacement[0]) + sqr(displacement[1]) + sqr(displacement[2]));
    if (norm == 0) {    // nodedata velocity can be exactly zero (e.g., in the initial state), thus must protect ourselves
        displacement[0] = basic_dx;     // choose displacement pointing in x-direction (arbitrarily)
        displacement[1] = displacement[2] = 0;
    } else {    // the usual case: displacement is directed along velocity, but length is basic_dx
        norm = basic_dx/norm;
        for (d=0; d<3; d++) displacement[d]*= norm;
    }
    const gridreal upwind_x = centroid[0] - displacement[0];
    const gridreal upwind_y = centroid[1] - displacement[1];
    const gridreal upwind_z = centroid[2] - displacement[2];
    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!c) {
            errorlog << "Tnode::CN_donor1: nonexistent cell at " << Tr3v(centroid).toString() << ", upwind_x=" << upwind_x/3e6 << "\n";
            continue;
        }
        // notice: this is not strictly correct in a T-junction node
        const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
        //      const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
        const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
        result[0]+= weight*c->celldata[cs][0];
        result[1]+= weight*c->celldata[cs][1];
        result[2]+= weight*c->celldata[cs][2];
        weightsum+= weight;
    }
    const real invweightsum = 1.0/weightsum;
    for (d=0; d<3; d++) nodedata[ns][d] = invweightsum*result[d];
}

//! Set resistivity at a node
void Tgrid::Tnode::set_resistivity(ResistivityProfile res)
{
    eta = res.getValue(centroid);
    if (eta<0) {
        errorlog << "ERROR: Negative resistivity at " << Tr3v(centroid).toString()
                 << "\n Quitting!";
        doabort();
    }
}

//! Calculate the electric field at a node
void Tgrid::Tnode::calc_E1(void)
{
    // E = -U_e x (B + B_0) + eta*J
    real B0[3] = {0.0, 0.0, 0.0};
    // Add constant magnetic field B0 (centroid = exact node coordinates)
    addConstantMagneticField(centroid, B0);
    const real Bx = nodedata[NODEDATA_B][0] + B0[0];
    const real By = nodedata[NODEDATA_B][1] + B0[1];
    const real Bz = nodedata[NODEDATA_B][2] + B0[2];
    const real uex = nodedata[NODEDATA_UE][0];
    const real uey = nodedata[NODEDATA_UE][1];
    const real uez = nodedata[NODEDATA_UE][2];
    nodedata[NODEDATA_E][0] = By*uez - Bz*uey;
    nodedata[NODEDATA_E][1] = Bz*uex - Bx*uez;
    nodedata[NODEDATA_E][2] = Bx*uey - By*uex;
    const real jx = nodedata[NODEDATA_J][0];
    const real jy = nodedata[NODEDATA_J][1];
    const real jz = nodedata[NODEDATA_J][2];
    nodedata[NODEDATA_E][0] += eta*jx;
    nodedata[NODEDATA_E][1] += eta*jy;
    nodedata[NODEDATA_E][2] += eta*jz;
    // Check electric field cut value (CONSTRAINT)
    if(Params::Ecut > 0) {
        real E = sqrt( sqr(nodedata[NODEDATA_E][0]) + sqr(nodedata[NODEDATA_E][1]) + sqr(nodedata[NODEDATA_E][2]) );
        if(E > Params::Ecut) {
            real scaling = Params::Ecut/E;
            nodedata[NODEDATA_E][0] *= scaling;
            nodedata[NODEDATA_E][1] *= scaling;
            nodedata[NODEDATA_E][2] *= scaling;
#ifndef NO_DIAGNOSTICS
            // Increase counter
            Tgrid::fieldCounter.cutRateE += 1.0;
#endif
        }
    }
}

//! cell2node interpolation (recursive)
void Tgrid::Tcell::CN_recursive(TCellDataSelect cs, TNodeDataSelect ns)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->CN_recursive(cs,ns);
    } else {
        // NODE LOOP
        // Idea: To pass through all nodes, pass through all right-pointing faces,
        // and all upper-right corner points thereof
        // IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
        int d = 0; // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1(cs,ns);
        } else {
            face[d][1]->node[2]->CN1(cs,ns);    // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1(cs,ns);
            }
    }
}

//! cell2node interpolation for smoothing (recursive)
void Tgrid::Tcell::CN_smoothing_recursive()
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->CN_smoothing_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_smoothing();
        } else {
            face[d][1]->node[2]->CN1_smoothing(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_smoothing();
            }
    }
}

//! cell2node interpolation of rho_q (recursive)
void Tgrid::Tcell::CN_rhoq_recursive()
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->CN_rhoq_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_rhoq();
        } else {
            face[d][1]->node[2]->CN1_rhoq(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_rhoq();
            }
    }
}

//! Upwind nodedata by using cell2node interpolation (recursive)
void Tgrid::Tcell::CN_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->CN_donor_recursive(cs,ns,uns,dt);
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN_donor1(cs,ns,uns,dt);
        } else {
            face[d][1]->node[2]->CN_donor1(cs,ns,uns,dt);
            //assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN_donor1(cs,ns,uns,dt);
            }
    }
}

//! Calculate the electric field at a node (recursive)
void Tgrid::Tcell::calc_node_E_recursive(void)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->calc_node_E_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_E1();
        } else {
            face[d][1]->node[2]->calc_E1();
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_E1();
            }
    }
}

//! Calculate the polarization electric field in a cell and store in CELLDATA_TEMP2 (recursive)
void Tgrid::Tcell::calc_cell_E_recursive(void)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->calc_cell_E_recursive();
    } else {
        celldata[CELLDATA_TEMP2][0] = 0.;
        celldata[CELLDATA_TEMP2][1] = 0.;
        celldata[CELLDATA_TEMP2][2] = 0.;
        //electron pressure contribution
        if(r2>sqr(Params::R_zeroPolarizationField)) { //if the cell is outside of the background ionosphere density peak, include the polarization electric field.
            const real pressure_coef=Params::k_B*Params::Te/Params::e;//KTe/e for isothermal plasma.
            celldata[CELLDATA_TEMP2][0] -= pressure_coef*celldata[CELLDATA_TEMP1][0]/rho_q;
            celldata[CELLDATA_TEMP2][1] -= pressure_coef*celldata[CELLDATA_TEMP1][1]/rho_q;
            celldata[CELLDATA_TEMP2][2] -= pressure_coef*celldata[CELLDATA_TEMP1][2]/rho_q;
        }
    }
}

//! Reset nc, rho_q and CELLDATA_Ji in a cell (recursive)
void Tgrid::Tcell::zero_rhoq_nc_Vq_recursive()
{
    nc = rho_q = 0.0;
    celldata[CELLDATA_Ji][0] = 0.0;
    celldata[CELLDATA_Ji][1] = 0.0;
    celldata[CELLDATA_Ji][2] = 0.0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->zero_rhoq_nc_Vq_recursive();
    }
}

//! Calculate Ue in a cell (recursive)
void Tgrid::Tcell::calc_ue_recursive(void)
{
    if (haschildren) {
        for (int ch=0; ch<8; ch++) child[0][0][ch]->calc_ue_recursive();
    } else {
        // Check zero field radius
        if (r2 < Params::R_zeroFields2) {
            for (int d=0; d<3; d++) {
                celldata[CELLDATA_UE][d] = 0;
            }
            return;
        }
        // Charge density in thel cell
        real chargedensity = rho_q;
        const real invrho_q = 1.0/chargedensity;
        real ue[3];
        real ue2 = 0.0;
        for (int d=0; d<3; d++) {
#ifndef IGNORE_ELECTRIC_FIELD_HALL_TERM
            ue[d] = (celldata[CELLDATA_Ji][d] - celldata[CELLDATA_J][d])*invrho_q;
#else
            ue[d] = (celldata[CELLDATA_Ji][d])*invrho_q;
#endif
            ue2+= sqr(ue[d]);
        }
        // Check maximum electron fluid velocity (CONSTRAINT)
        if (Params::Ue_max > 0 && ue2 > Params::Ue_max2) {
            const real norm = Params::Ue_max/sqrt(ue2);
            for (int d=0; d<3; d++) ue[d]*= norm;
#ifndef NO_DIAGNOSTICS
            // Increase counter
            Tgrid::fieldCounter.cutRateUe += 1.0;
#endif
        }
        for (int d=0; d<3; d++) celldata[CELLDATA_UE][d] = ue[d];
    }
}

//! Calculate curl at cell faces (recursive)
void Tgrid::Tcell::FaceCurl_recursive(TNodeDataSelect nsB, TFaceDataSelect fsj, int d, real factor)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->FaceCurl_recursive(nsB,fsj,d,factor);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->Curl1(nsB,fsj,d,0.5*size,factor);
        } else {
            face[d][dir]->Curl1(nsB,fsj,d,size,factor);
        }
    }
}

//! node2face interpolation (recursive)
void Tgrid::Tcell::NF_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->NF_recursive(ns,fs,d);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->NF1(ns,fs,d);
        } else {
            face[d][dir]->NF1(ns,fs,d);
        }
    }
}

//! node2face interpolation of rho_q (recursive)
void Tgrid::Tcell::NF_rhoq_recursive(int d)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->NF_rhoq_recursive(d);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->NF_rhoq1();
        } else {
            face[d][dir]->NF_rhoq1();
        }
    }
}

//! Face propagate (recursive)
void Tgrid::Tcell::FacePropagate_recursive(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->FacePropagate_recursive(Bold,Bnew,dt);
    } else {
        int d;
        for (d=0; d<3; d++)
            if (isrefined_face(d,1)) {
                int f;
                for (f=0; f<4; f++) refintf[d][1]->face[f]->Propagate1(Bold,Bnew,dt);
            } else {
                face[d][1]->Propagate1(Bold,Bnew,dt);
            }
    }
}

//! Set resistivity at a node (recursive)
void Tgrid::Tcell::set_resistivity_recursive(ResistivityProfile res)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->set_resistivity_recursive(res);
    } else {
        // NODE LOOP
        int d = 0; // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->set_resistivity(res);
        } else {
            face[d][1]->node[2]->set_resistivity(res);
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->set_resistivity(res);
            }
    }
}

//! Set magnetic field (B1) on a face (recursive)
void Tgrid::Tcell::set_B_recursive(void (*func)(const gridreal[3], datareal[3]))
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->set_B_recursive(func);
    } else {
        int d;
        gridreal r[3];
        for (d=0; d<3; d++) {
            r[0]=centroid[0];
            r[1]=centroid[1];
            r[2]=centroid[2];
            r[d]+= 0.5*size;
            if (isrefined_face(d,1)) {
                int f;
                for (f=0; f<4; f++) {
                    // ripped from g.recoursen
                    const int f_plusminus_y[4] = {-1,+1,+1,-1};  // cyclic order of nodes in face
                    const int f_plusminus_z[4] = {-1,-1,+1,+1};
                    for (f=0; f<4; f++) {
                        r[(d+1)%3]+= (0.25*size)*f_plusminus_y[f];
                        r[(d+2)%3]+= (0.25*size)*f_plusminus_z[f];
                        refintf[d][1]->face[f]->set_B1(r,func,d);
                    }
                }
            } else {
                face[d][1]->set_B1(r,func,d);
            }
        }
    }
}

//! Set background charge density in a cell (recursive)
void Tgrid::Tcell::set_bgRhoQ_recursive(BackgroundChargeDensityProfile func)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->set_bgRhoQ_recursive(func);
    } else {
        this->rho_q_bg = func.getValue(centroid);
    }
}

//! Constructor for the magnetic log
MagneticLog::MagneticLog()
{
    avgBx = 0.0;
    avgBy = 0.0;
    avgBz = 0.0;
    avgB = 0.0;
    avgDivB = 0.0;
    energyB = 0.0;
    maxDivB = 0.0;
    maxB = 0.0;
    posMaxB[0] = 0.0;
    posMaxB[1] = 0.0;
    posMaxB[2] = 0.0;
    maxDxDivBperB = 0.0;
}

//! Calculate magnetic log values
void Tgrid::calc_facediv(TFaceDataSelect fs, MagneticLog& result) const
{
    MagneticLog BLog;
    // Loop thru the grid cells
    int i,j,k;
    ForInterior(i,j,k) {
        int c = flatindex(i,j,k);
        cells[c]->calc_facediv_recursive(fs,BLog);
        // Comparate current max values to new ones
        if(BLog.maxDivB > result.maxDivB) {
            result.maxDivB = BLog.maxDivB;
        }
        if(BLog.maxB > result.maxB) {
            result.maxB = BLog.maxB;
            result.posMaxB[0] = BLog.posMaxB[0];
            result.posMaxB[1] = BLog.posMaxB[1];
            result.posMaxB[2] = BLog.posMaxB[2];
        }
        if(BLog.maxDxDivBperB > result.maxDxDivBperB) {
            result.maxDxDivBperB = BLog.maxDxDivBperB;
        }
        // Sum average values and energy
        result.avgDivB += BLog.avgDivB;
        result.avgBx += BLog.avgBx;
        result.avgBy += BLog.avgBy;
        result.avgBz += BLog.avgBz;
        result.avgB += BLog.avgB;
        result.energyB += BLog.energyB;
    }
    // Normalize average values
    real normCoeff = Params::nx * Params::ny * Params::nz;
    result.avgBx /= normCoeff;
    result.avgBy /= normCoeff;
    result.avgBz /= normCoeff;
    result.avgB /= normCoeff;
    result.avgDivB /= normCoeff;
    // Energy dimension
    result.energyB *= cube(bgdx)/(2 * Params::mu_0);
}

//! Calculate magnetic log values (recursive)
void Tgrid::Tcell::calc_facediv_recursive(TFaceDataSelect fs, MagneticLog& result) const
{
    if (haschildren) {
        MagneticLog BLog;
        // Loop thru the child cells
        for (int ch=0; ch<8; ch++) {
            child[0][0][ch]->calc_facediv_recursive(fs,BLog);
            // Comparate current max values to new ones
            if(BLog.maxDivB > result.maxDivB) {
                result.maxDivB = BLog.maxDivB;
            }
            if(BLog.maxB > result.maxB) {
                result.maxB = BLog.maxB;
                result.posMaxB[0] = BLog.posMaxB[0];
                result.posMaxB[1] = BLog.posMaxB[1];
                result.posMaxB[2] = BLog.posMaxB[2];
            }
            if(BLog.maxDxDivBperB > result.maxDxDivBperB) {
                result.maxDxDivBperB = BLog.maxDxDivBperB;
            }
            // Sum average values and energy
            result.avgBx += BLog.avgBx;
            result.avgBy += BLog.avgBy;
            result.avgBz += BLog.avgBz;
            result.avgB += BLog.avgB;
            result.avgDivB += BLog.avgDivB;
            result.energyB += BLog.energyB;
        }
        // Child cell normalization
        result.avgDivB *= 0.125;
        result.avgBx *= 0.125;
        result.avgBy *= 0.125;
        result.avgBz *= 0.125;
        result.avgB *= 0.125;
        result.energyB *= 0.125;
    } else {
        // Flux out of cell, aveBxyz and posMaxB
        real flux = 0.0;
        real tempAvgB[3];
        for (int d=0; d<3; d++) {
            flux += faceave(d,1,fs) - faceave(d,0,fs);
            tempAvgB[d] = 0.5 *(faceave(d,1,fs) + faceave(d,0,fs));
            result.posMaxB[d] = centroid[d];
        }
        result.avgBx = tempAvgB[0];
        result.avgBy = tempAvgB[1];
        result.avgBz = tempAvgB[2];
        // B^2 in this cell
        result.energyB = sqr(result.avgBx) + sqr(result.avgBy) + sqr(result.avgBz);
        // B in this cell
        result.maxB = result.avgB = sqrt(result.energyB);
        // divB in this cell
        result.avgDivB = result.maxDivB = fabs(flux)*invsize;
        // DxDivBperB in this cell
        if (result.avgB > 1e-20) {
            result.maxDxDivBperB = size*result.avgDivB/result.avgB;
        } else {
            result.maxDxDivBperB = 0.0;
        }
    }
}

//! Calculate gradient of rho_q (recursive)
void Tgrid::Tcell::CalcGradient_rhoq_recursive(void)
{
    if (haschildren) {
        for (int ch=0; ch<8; ch++) child[0][0][ch]->CalcGradient_rhoq_recursive();
    } else {
        for (int d=0; d<3; d++)celldata[CELLDATA_TEMP1][d]= (faceave(d,1,Tgrid::FACEDATA_MINUSDB) - faceave(d,0,Tgrid::FACEDATA_MINUSDB))/size;
    }
}

//! Return the number of refined faces in a cell
int Tgrid::Tcell::Nfaces() const
{
    if (haschildren) return 0;
    int d,dir,result=0;
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++)
            if (isrefined_face(d,dir)) result+= 4;
            else result++;
    return result;
}

//! Return the number of refined cells in a cell (recursive)
int Tgrid::Tcell::Ncells_recursive() const
{
    int result = 1;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) result+= child[0][0][ch]->Ncells_recursive();
    }
    return result;
}

//! Calculate curl on a face
void Tgrid::Tface::Curl1(TNodeDataSelect nsB, TFaceDataSelect fsj, int d, gridreal dx, real factor)
{
    // Use nodedata nsB to calculate fsj = curl(nsB)/mu0 on a face
    // <j> = (1/S)*int(dS.j) = (1/(S*mu0))int(dl.B) = (1/(dx*mu0))*int(dl.B)|.dl=1
    const int diry = (d+1)%3;
    const int dirz = (d+2)%3;
    const real twice_integral_01_y =
        sidenode[0]
        ? 0.5*(node[0]->nodedata[nsB][diry] + 2*sidenode[0]->nodedata[nsB][diry] + node[1]->nodedata[nsB][diry])
        : node[0]->nodedata[nsB][diry] + node[1]->nodedata[nsB][diry];
    const real twice_integral_12_z =
        sidenode[1]
        ? 0.5*(node[1]->nodedata[nsB][dirz] + 2*sidenode[1]->nodedata[nsB][dirz] + node[2]->nodedata[nsB][dirz])
        : node[1]->nodedata[nsB][dirz] + node[2]->nodedata[nsB][dirz];
    const real twice_integral_23_y =
        sidenode[2]
        ? 0.5*(node[2]->nodedata[nsB][diry] + 2*sidenode[2]->nodedata[nsB][diry] + node[3]->nodedata[nsB][diry])
        : node[2]->nodedata[nsB][diry] + node[3]->nodedata[nsB][diry];
    const real twice_integral_30_z =
        sidenode[3]
        ? 0.5*(node[3]->nodedata[nsB][dirz] + 2*sidenode[3]->nodedata[nsB][dirz] + node[0]->nodedata[nsB][dirz])
        : node[3]->nodedata[nsB][dirz] + node[0]->nodedata[nsB][dirz];
    const real twice_integral = twice_integral_01_y + twice_integral_12_z - twice_integral_23_y - twice_integral_30_z;
    facedata[fsj] = factor*(0.5/dx)*twice_integral; // notice factor - (1 for Faraday, 1/Params::mu_0 for Ampere)
}

//! node2face interpolation
void Tgrid::Tface::NF1(TNodeDataSelect ns, TFaceDataSelect fs, int d)
{
    // Face number (0,1,2)
    //         (0)  F01_dL01  (1)
    //          -------x-------
    //          |             |
    //          |             |
    // F30_dL30 x             x F12_dL12
    //          |             |
    //          |             |
    //          -------x-------
    //         (3)  F23_dL23  (2)
    const real node_0_nodedata = node[0]->nodedata[ns][d];
    const real node_1_nodedata = node[1]->nodedata[ns][d];
    const real node_2_nodedata = node[2]->nodedata[ns][d];
    const real node_3_nodedata = node[3]->nodedata[ns][d];
    const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata);
    const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata);
    const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata);
    const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata);
    const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/4;
    facedata[fs] = integral;
}

//! node2face interpolation or rho_q
void Tgrid::Tface::NF_rhoq1(void)
{
    const real twice_integral_01 =
        sidenode[0]
        ? 0.5*(node[0]->nodedata[NODEDATA_UE][0] + 2*sidenode[0]->nodedata[NODEDATA_UE][0] + node[1]->nodedata[NODEDATA_UE][0])
        : node[0]->nodedata[NODEDATA_UE][0] + node[1]->nodedata[NODEDATA_UE][0];
    const real twice_integral_12 =
        sidenode[1]
        ? 0.5*(node[1]->nodedata[NODEDATA_UE][0] + 2*sidenode[1]->nodedata[NODEDATA_UE][0] + node[2]->nodedata[NODEDATA_UE][0])
        : node[1]->nodedata[NODEDATA_UE][0] + node[2]->nodedata[NODEDATA_UE][0];
    const real twice_integral_23 =
        sidenode[2]
        ? 0.5*(node[2]->nodedata[NODEDATA_UE][0] + 2*sidenode[2]->nodedata[NODEDATA_UE][0] + node[3]->nodedata[NODEDATA_UE][0])
        : node[2]->nodedata[NODEDATA_UE][0] + node[3]->nodedata[NODEDATA_UE][0];
    const real twice_integral_30 =
        sidenode[3]
        ? 0.5*(node[3]->nodedata[NODEDATA_UE][0] + 2*sidenode[3]->nodedata[NODEDATA_UE][0] + node[0]->nodedata[NODEDATA_UE][0])
        : node[3]->nodedata[NODEDATA_UE][0] + node[0]->nodedata[NODEDATA_UE][0];
    const real twice_integral = twice_integral_01 + twice_integral_12 + twice_integral_23 + twice_integral_30;
    facedata[FACEDATA_MINUSDB] = 0.125*twice_integral;
}

//! face2r interpolation
void Tgrid::faceintpol(const shortreal r[3], TFaceDataSelect s, real result[3])
{
    int d;
    gridreal t,lowercorner[3];
    Tcell *const c = findcell(r,lowercorner);
    if (!c) {
        errorlog << "ERROR [Tgrid::faceintpol]: findcell returned null for r=" << Tr3v(r).toString() << "\n";
        doabort();
    }
    if (c->anyrefined_face()) {
        for (d=0; d<3; d++) {
            t = (r[d] - lowercorner[d])*c->invsize;
            result[d] = (1-t)*c->faceave(d,0,s) + t*c->faceave(d,1,s);
        }
    } else {
        // faster branch for fully regular case
        t = (r[0] - lowercorner[0])*c->invsize;
        result[0] = (1-t)*c->face[0][0]->facedata[s] + t*c->face[0][1]->facedata[s];
        t = (r[1] - lowercorner[1])*c->invsize;
        result[1] = (1-t)*c->face[1][0]->facedata[s] + t*c->face[1][1]->facedata[s];
        t = (r[2] - lowercorner[2])*c->invsize;
        result[2] = (1-t)*c->face[2][0]->facedata[s] + t*c->face[2][1]->facedata[s];
    }
    saved_cellptr = c;
}

//! face2cell interpolation
void Tgrid::FC(TFaceDataSelect fs, TCellDataSelect cs)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->FC_recursive(fs,cs);
    }
}

//! node2cell interpolation for smoothing
void Tgrid::NC_smoothing()
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->NC_smoothing_recursive();
    }
}

//! node2cell interpolation
void Tgrid::NC(TNodeDataSelect ns,TCellDataSelect cs)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->NC_recursive(ns,cs);
    }
}

//! cell2node interpolation
void Tgrid::CN(TCellDataSelect cs, TNodeDataSelect ns)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->CN_recursive(cs,ns);
            }
}

//! cell2node interpolation for smoothing
void Tgrid::CN_smoothing()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->CN_smoothing_recursive();
            }
}

//! cell2node interpolation of rho_q
void Tgrid::CN_rhoq()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->CN_rhoq_recursive();
            }
}

//! Set resistivity at a node
void Tgrid::set_resistivity(ResistivityProfile res)
{
    MSGFUNCTIONCALL("Tgrid::set_resistivity");
    int i,j,k;
    for (i=0; i<nx-1; i++)
        for (j=0; j<ny-1; j++)
            for (k=0; k<nz-1; k++) {
                cells[flatindex(i,j,k)]->set_resistivity_recursive(res);
            }
    MSGFUNCTIONEND("Tgrid::set_resistivity");
}

//! Upwind nodedata by using cell2node interpolation
void Tgrid::CN_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->CN_donor_recursive(cs,ns,uns,dt);
            }
}

//! Reset nc, rho_q and CELLDATA_Ji in a cell
void Tgrid::zero_rhoq_nc_Vq()
{
    int i,j,k,c;
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->zero_rhoq_nc_Vq_recursive();
    }
}

//! Accumulate Particle-In-Cell quantities in the grid (recursive)
inline gridreal Tgrid::accumulate_PIC_recursive(const TBoxDef& cloudbox, Tcell *c, gridreal invvol, const shortreal v[3], real w, int popid)
{
    // This function is recursive, but a good compiler will inline the first-level calls in accumulate_PIC anyway.
    gridreal accum = 0.0;
    if (c->haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) accum+= accumulate_PIC_recursive(cloudbox,c->child[0][0][ch],invvol,v,w,popid);
    } else {
        TBoxDef cellbox;
        const gridreal halfdx = 0.5*size(c);
        cellbox.lowx = c->centroid[0] - halfdx;
        cellbox.lowy = c->centroid[1] - halfdx;
        cellbox.lowz = c->centroid[2] - halfdx;
        cellbox.size = size(c);
        accum = intersection_volume(cloudbox,cellbox)*invvol;
        // Particle number contribution from a macroparticle to this cell
        register const real accum_w = w*accum;
        // Charge contribution from a macroparticle to this cell
        datareal charge = accum_w*Params::pops[popid]->q;
        if(Params::pops[popid]->getAccumulate() == true) {
            // Add particle number contribution to the cell
            c->nc += accum_w;

            // Add charge contribution to the cell
            c->rho_q += charge;

            // Add ion current contribution to the cell
            for (int d=0; d<3; d++) {
                c->celldata[CELLDATA_Ji][d] += charge*v[d];
            }
        } else {
            accum = 0.0;
        }

        if (Params::averaging == true) {
#ifdef SAVE_POPULATION_AVERAGES
            c->pop_ave_n[popid] += accum_w;
            c->pop_ave_vx[popid] += accum_w*v[0];
            c->pop_ave_vy[popid] += accum_w*v[1];
            c->pop_ave_vz[popid] += accum_w*v[2];
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
            const real v2 = vecsqr(v);
            const real spectraAccum = sqrt(v2)*accum_w;
            const unsigned int Nbins = Params::spectraV2BinsPerPop[popid].size() - 1;
            if(v2 < Params::spectraV2BinsPerPop[popid][0] && Params::spectraEminAll == true) {
                c->spectra[popid][0] += spectraAccum;
            } else if(v2 > Params::spectraV2BinsPerPop[popid][Nbins] && Params::spectraEmaxAll == true) {
                c->spectra[popid][Nbins] += spectraAccum;
            } else {
                for(unsigned int i=0; i<Nbins; ++i) {
                    if(v2 >= Params::spectraV2BinsPerPop[popid][i] && v2 < Params::spectraV2BinsPerPop[popid][i+1]) {
                        c->spectra[popid][i] += spectraAccum;
                        break;
                    }
                }
            }
#endif
        }
    }
    return accum;
}

//! Accumulate Particle-In-Cell quantities in the grid
void Tgrid::accumulate_PIC(const shortreal r[3], const shortreal v[3], real w, int popid)
{
    gridreal centroid[3];
    Tcell *const c = findcell(r,centroid);
    if (!c) {
        // Do not abort if no cell is found. Particle removed in pass_with_relocate afterwards.
        errorlog << "ERROR [Tgrid::accumulate_PIC]: findcell returned null for r="
                 << Tr3v(r).toString() << "\n"
                 << "   v=" << Tr3v(v).toString()
                 << " popIdStr=" << Params::pops[popid]->getIdStr() << endl;
        return; //doabort();
    }
    const gridreal halfdx = 0.5*size(c);
    bool movetoright[3];        // false if move to left, true if move to right, along dimension d
    for (int d=0; d<3; d++) {
        centroid[d]+= halfdx;
        //if (r[d]==centroid[d]) errorlog << "Problem:::"<<endl;
        movetoright[d] = (r[d] > centroid[d]);
    }
    Tcell *C[8];        // eight cells that define the stencil, some of them can be refined
    C[0] = c;
    C[1] = moveto(c,0,movetoright[0]);      // move along x, right or left
    C[2] = moveto(c,1,movetoright[1]);      // move along y, right or left
    C[3] = moveto(C[1],1,movetoright[1]);
    C[4] = moveto(c,2,movetoright[2]);      // move along z, right or left
    C[5] = moveto(C[4],0,movetoright[0]);
    C[6] = moveto(C[4],1,movetoright[1]);
    C[7] = moveto(C[5],1,movetoright[1]);
    // Check if all cells are of the same level. This is a necessary requirement for "regularity" (see below)
    int a;
    bool all_same_level = true;
    const int clevel = c->level;
    for (a=1; a<8; a++) if (C[a]->level != clevel) {
            all_same_level = false;
            break;
        }
    // Replace duplicate cells by null pointers. Duplicates can occur if the neighbour is coarser.
    // Cells 0,1,2,4 cannot be duplicates anyway so exclude them from the search.
    int b;
    bool found;
    // We want 'a' to range 3,5,6,7. It is most efficient to write it all out explicitly.
    a = 3;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 5;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 6;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 7;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    const gridreal cloudsize0 = size(c);  // initially, the cloud is of the same size as cell c
    // C1,C2,C4 are the direct neighbours, check if any of them is refined
    // Now checks also the diagonal neighbours and does not scale linearly
    // as before but with a step directly to halfsize that's the size in
    // the refinement. Is a correct way of handling this!
    const bool C1haschildren = C[1]->haschildren;
    const bool C2haschildren = C[2]->haschildren;
    const bool C3haschildren = C[3] && C[3]->haschildren;
    const bool C4haschildren = C[4]->haschildren;
    const bool C5haschildren = C[5] && C[5]->haschildren;
    const bool C6haschildren = C[6] && C[6]->haschildren;
    const bool C7haschildren = C[7] && C[7]->haschildren;
    const bool any_refined = (C1haschildren || C2haschildren || C3haschildren
                              || C4haschildren || C5haschildren || C6haschildren || C7haschildren);
    gridreal cloudsize;
    if (any_refined) {
        cloudsize = halfdx;
    } else {
        cloudsize = cloudsize0;
    }
    const gridreal halfcloudsize = 0.5*cloudsize;
    const gridreal invvol = 1.0/(cloudsize*cloudsize*cloudsize);
    gridreal accum = 0.0;
    // Pass through the eight cells. Exclude nulls (those duplicates we removed above).
    if (!any_refined && all_same_level) {
        // Omit the halfdx shift in intersection box definitions, can be done since intersection volume is translation invariant
        // the essential thing is that cloudsize and cellboxsize are the same in the regular case
        TCellPtr c1;
        for (a=0; a<8; a++) {
            c1 = C[a];
            gridreal accum1 = intersection_volume_samesize_nochecks(r,c1->centroid,cloudsize)*invvol;
            // Handle correctly boundary cases by "reflecting" ghost cell density
            // into corresponding interior cell. Algorithm:
            // (1) find the cell's i,j,k basegrid indices from its flatindex,
            // (2) if i=0 set i=1, if i=nx-1 set i=nx-2, etc., (3) reassign it from basegrid
            // Notice that this is in the 'regular grid' branch only, thus it works correctly only if
            // grid refinement does not touch the box boundary (doing so would cause other problems
            // in any case so this is not a new restriction).
            if (c1->level == 0) {
                int i,j,k;
                decompose(c1->flatind,i,j,k);		// find (i,j,k) indices of the basegrid cell
                bool anymod = false;
                if (i==0) {
                    i = 1;    // modify it if needed
                    anymod=true;
                }
                if (i==nx-1) {
                    i = nx-2;
                    anymod=true;
                }
                if (j==0) {
                    j = 1;
                    anymod=true;
                }
                if (j==ny-1) {
                    j = ny-2;
                    anymod=true;
                }
                if (k==0) {
                    k = 1;
                    anymod=true;
                }
                if (k==nz-1) {
                    k = nz-2;
                    anymod=true;
                }
                if (anymod) {
                    c1 = cells[flatindex(i,j,k)];	// compute back the cell from the basegrid
                }
            }
            // Particle number contribution from a macroparticle to this cell
            register const real accum_w = w*accum1;
            // Charge contribution from a macroparticle to this cell
            datareal charge = accum_w*Params::pops[popid]->q;
            if(Params::pops[popid]->getAccumulate() == true) {
                // Add particle number contribution to the cell
                c1->nc += accum_w;
                // Add charge contribution to the cell
                c1->rho_q += charge;
                // Add charge times velocity contribution to the cell
                for (int d=0; d<3; d++) {
                    c1->celldata[CELLDATA_Ji][d] += charge*v[d];
                }
            } else {
                accum1 = 0.0;
            }
            if(Params::averaging == true) {
#ifdef SAVE_POPULATION_AVERAGES
                c1->pop_ave_n[popid] += accum_w;
                c1->pop_ave_vx[popid] += accum_w*v[0];
                c1->pop_ave_vy[popid] += accum_w*v[1];
                c1->pop_ave_vz[popid] += accum_w*v[2];
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
                const real v2 = vecsqr(v);
                const real spectraAccum = sqrt(v2)*accum_w;
                const unsigned int Nbins = Params::spectraV2BinsPerPop[popid].size()-1;
                if(v2 < Params::spectraV2BinsPerPop[popid][0]) {
                    c1->spectra[popid][0] += spectraAccum;
                } else if(v2 > Params::spectraV2BinsPerPop[popid][Nbins]) {
                    c1->spectra[popid][Nbins] += spectraAccum;
                } else {
                    for(unsigned int i=0; i<Nbins; ++i) {
                        if(v2 >= Params::spectraV2BinsPerPop[popid][i] && v2 < Params::spectraV2BinsPerPop[popid][i+1]) {
                            c1->spectra[popid][i] += spectraAccum;
                            break;
                        }
                    }
                }
#endif
            }
            accum += accum1;
        }
    } else {
        TBoxDef cloudbox;
        cloudbox.lowx = r[0] - halfcloudsize;
        cloudbox.lowy = r[1] - halfcloudsize;
        cloudbox.lowz = r[2] - halfcloudsize;
        cloudbox.size = cloudsize;
        for (a=0; a<8; a++) if (C[a]) {
                accum+= accumulate_PIC_recursive(cloudbox,C[a],invvol,v,w,popid);
            }
    }
}

//! Finalize accumulate Particle-In-Cell quantities in the grid (recursive)
void Tgrid::finalize_accum_recursive(Tcell *c)
{
    if (c->haschildren) {
        for (int ch=0; ch<8; ch++) finalize_accum_recursive(c->child[0][0][ch]);
    } else {
        // Cell volume
        gridreal invvol = 1./volume(c);
        // Normalize particle number in the cell by the cell volume -> particle number density
        c->nc *= invvol;
        // Normalize charge in the cell by the cell volume -> charge density
        c->rho_q *= invvol;
        // Normalize ion current in the cell by the cell volume -> ion current density
        for (int d=0; d<3; d++) {
            c->celldata[CELLDATA_Ji][d] *= invvol;
        }
        // Add background charge density
        c->rho_q += c->rho_q_bg;
        // Check charge density min value (CONSTRAINT)
        if ( c->rho_q < Params::rho_q_min ) {
            c->rho_q = Params::rho_q_min;
#ifndef NO_DIAGNOSTICS
            // Increase counter
            Tgrid::fieldCounter.cutRateRhoQ += 1.0;
#endif
        }
        // Averaging
        if (Params::averaging == true) {
            c->ave_nc+= c->nc;
            for (int d=0; d<3; d++) {
                if (c->neighbour[d][1] == 0) {
                    continue;
                }
                if (c->isrefined_face(d,1)) {
                    for (int f=0; f<4; f++) {
                        c->refintf[d][1]->face[f]->facedata[FACEDATA_AVEB] += c->refintf[d][1]->face[f]->facedata[FACEDATA_B];
                    }
                } else {
                    c->face[d][1]->facedata[FACEDATA_AVEB] += c->face[d][1]->facedata[FACEDATA_B];
                }
            }
        }
    }
}

//! Finalize accumulate Particle-In-Cell quantities in the grid
void Tgrid::finalize_accum()
{
    Neumann_rhoq();
    int i,j,k,c;
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        finalize_accum_recursive(cells[c]);
    }
    if (Params::averaging == true) {
        ave_ntimes++;
    }
}

//! Copy cell data
void Tgrid::copy_celldata(int cTo, int cFrom, TCellDataSelect cs)
{
    int d;
    if (cells[cFrom]->haschildren) {
        real avevec[3];
        cells[cFrom]->childave(cs, avevec);
        for (d=0; d<3; d++) cells[cTo]->celldata[cs][d] = avevec[d];
    } else {
        for (d=0; d<3; d++)
            cells[cTo]->celldata[cs][d] = cells[cFrom]->celldata[cs][d];
    }
}

//! Copy cell rho_q
void Tgrid::copy_rhoq(int cTo, int cFrom)
{
    if (cells[cFrom]->haschildren) {
        cells[cTo]->rho_q = cells[cFrom]->childave_rhoq();
    } else {
        cells[cTo]->rho_q = cells[cFrom]->rho_q;
    }
}

//! Copy cell data for smoothing
void Tgrid::copy_smoothing(int cTo, int cFrom)
{
    int d;
    if (cells[cFrom]->haschildren) {
        real avevec[3];
        cells[cFrom]->childave(CELLDATA_Ji, avevec);
        for (d=0; d<3; d++) cells[cTo]->celldata[CELLDATA_Ji][d] = avevec[d];//ion current density
        cells[cTo]->rho_q = cells[cFrom]->childave_rhoq();//charge density
        cells[cTo]->nc = cells[cFrom]->childave_nc();//number density
    } else {
        for (d=0; d<3; d++)
            cells[cTo]->celldata[CELLDATA_Ji][d] = cells[cFrom]->celldata[CELLDATA_Ji][d];//ion current density
        cells[cTo]->rho_q = cells[cFrom]->rho_q;//charge density
        cells[cTo]->nc = cells[cFrom]->nc;//number density
    }
}

//! Calculate Ue
void Tgrid::calc_ue(void)
{
    //! Ue = (Ji - j)/rho_q
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->calc_ue_recursive();
    }
}

//! Calculate electric field at nodes
void Tgrid::calc_node_E(void)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->calc_node_E_recursive();
            }
}

//! Calculate electric field in cells
void Tgrid::calc_cell_E(void)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->calc_cell_E_recursive();
    }
}

//! Calculate curl on cell faces
void Tgrid::FaceCurl(TNodeDataSelect nsB, TFaceDataSelect fsj, real factor)
{
    //1. Ampere's law j=curl(B)/mu0  2. Faraday's induction dB/dt=-curl(E)
    // Factor is for case 1: 1/Params::mu_0  and for case 2: 1
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->FaceCurl_recursive(nsB,fsj,0, factor);
                if (i > 0 && k > 0) cells[c]->FaceCurl_recursive(nsB,fsj,1, factor);
                if (i > 0 && j > 0) cells[c]->FaceCurl_recursive(nsB,fsj,2, factor);
            }
}

//! node2face interpolation
void Tgrid::NF(TNodeDataSelect ns, TFaceDataSelect fs)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->NF_recursive(ns,fs,0);
                if (i > 0 && k > 0) cells[c]->NF_recursive(ns,fs,1);
                if (i > 0 && j > 0) cells[c]->NF_recursive(ns,fs,2);
            }
}

// node2face interpolation of rho_q
void Tgrid::NF_rhoq()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->NF_rhoq_recursive(0);
                if (i > 0 && k > 0) cells[c]->NF_rhoq_recursive(1);
                if (i > 0 && j > 0) cells[c]->NF_rhoq_recursive(2);
            }
}

//! Set magnetic field in cells and faces
void Tgrid::set_B(void (*f)(const gridreal r[3], datareal[3]))
{
    MSGFUNCTIONCALL("Tgrid::set_B");
    // Set B-field at all the faces:
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->set_B_recursive(f);
            }
    // Interpolate it also to cells:
    FC(FACEDATA_B,CELLDATA_B);
    MSGFUNCTIONEND("Tgrid::set_B");
}

//! Set background charge density in cells
void Tgrid::set_bgRhoQ(BackgroundChargeDensityProfile func)
{
    MSGFUNCTIONCALL("Tgrid::set_bgRhoQ");
    // Set rho_q_bg-field in all cells:
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->set_bgRhoQ_recursive(func);
            }
    MSGFUNCTIONEND("Tgrid::set_bgRhoQ");
}

//! Face propagate
void Tgrid::FacePropagate(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->FacePropagate_recursive(Bold,Bnew,dt);
            }
}
//! Calculate the gradient of rho_q
void Tgrid::CalcGradient_rhoq(void)
{
    int i,j,k;
    ForInterior(i,j,k) {
        cells[flatindex(i,j,k)]->CalcGradient_rhoq_recursive();
    }
}

//! Return the number of cells with ghosts included
int Tgrid::Ncells_with_ghosts() const
{
    const int Nbase = nx*ny*nz;
    int result = 0;
    int i;
    for (i=0; i<Nbase; i++) result+= cells[i]->Ncells_recursive();
    return result;
}

//! Return the number of cells without ghosts
int Tgrid::Ncells_without_ghosts() const
{
    int result = 0;
    int i,j,k;
    ForInterior(i,j,k) result+= cells[flatindex(i,j,k)]->Ncells_recursive();
    return result;
}

//! Neumann boundary conditions
void Tgrid::Neumann(TCellDataSelect cs)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i+1,j,k),cs);
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i-1,j,k),cs);
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
#ifdef PERIODIC_FIELDS_Y
            copy_celldata(flatindex(i,j,k),flatindex(i,ny-2,k),cs);
#else
            copy_celldata(flatindex(i,j,k),flatindex(i,j+1,k),cs);
#endif
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
#ifdef PERIODIC_FIELDS_Y
            copy_celldata(flatindex(i,j,k),flatindex(i,1,k),cs);
#else
            copy_celldata(flatindex(i,j,k),flatindex(i,j-1,k),cs);
#endif
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata(flatindex(i,j,k),flatindex(i,j,k+1),cs);
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata(flatindex(i,j,k),flatindex(i,j,k-1),cs);
}

//! Neumann boundary conditions of rho_q
void Tgrid::Neumann_rhoq()
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_rhoq(flatindex(i,j,k),flatindex(i+1,j,k));
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_rhoq(flatindex(i,j,k),flatindex(i-1,j,k));
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
#ifdef PERIODIC_FIELDS_Y
            copy_rhoq(flatindex(i,j,k),flatindex(i,ny-2,k));
#else
            copy_rhoq(flatindex(i,j,k),flatindex(i,j+1,k));
#endif
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
#ifdef PERIODIC_FIELDS_Y
            copy_rhoq(flatindex(i,j,k),flatindex(i,1,k));
#else
            copy_rhoq(flatindex(i,j,k),flatindex(i,j-1,k));
#endif
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_rhoq(flatindex(i,j,k),flatindex(i,j,k+1));
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_rhoq(flatindex(i,j,k),flatindex(i,j,k-1));
}

//! Neumann boundary conditions for smoothing
void Tgrid::Neumann_smoothing()
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_smoothing(flatindex(i,j,k),flatindex(i+1,j,k));
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_smoothing(flatindex(i,j,k),flatindex(i-1,j,k));
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_smoothing(flatindex(i,j,k),flatindex(i,j+1,k));
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_smoothing(flatindex(i,j,k),flatindex(i,j-1,k));
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_smoothing(flatindex(i,j,k),flatindex(i,j,k+1));
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_smoothing(flatindex(i,j,k),flatindex(i,j,k-1));
}

//! Smoothing of nc, rho_q and CELLDATA_Ji
void Tgrid::smoothing()
{
    for(int n=0; n<Params::densitySmoothingNumber; n++) {
        Neumann_smoothing();//set up Neumann boundary(can also be other boundary condition) for the particle related quantities
        CN_smoothing();//Cell to Node interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
        NC_smoothing();//Node to Cell interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
    }
    if(Params::densitySmoothingNumber > 0) {
        Neumann_smoothing();//set up Neumann boundary(can also be other boundary condition) for the particle related quantities
    }
    return;
}

//! Smoothing of electric field
void Tgrid::smoothing_E()
{
    for(int n=0; n<Params::electricFieldSmoothingNumber; n++) {
        NC(NODEDATA_E,CELLDATA_TEMP1);//Node to Cell interpolation. NODEDATA_E is interpolated to CELLDATA_TEMP1. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
        Neumann(CELLDATA_TEMP1);//Set up Neumann boundary(can also be other boundary condition) for CELLDATA_TEMP1 (actually saves CELLDATA_E).
        CN(CELLDATA_TEMP1,NODEDATA_E);//Cell to Node interpolation for E. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
    }
    return;
}

//! Field boundary conditions
void Tgrid::boundarypass(int dim, bool toRight, void (*funcB)(datareal cdata[NCELLDATA][3], int))
{
    int i,j,k;
    switch (dim) {
    case 0: // x-boundaries
        i = toRight ? nx-1 : 0;
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {
#ifdef RECONNECTION_GEOMETRY
                if(toRight == true) {
                    (*funcB)(cells[flatindex(i,j,k)]->celldata, 0);
                } else                {
                    (*funcB)(cells[flatindex(i,j,k)]->celldata, 3);
                }
#else
                (*funcB)(cells[flatindex(i,j,k)]->celldata, 0);
#endif
            }
        }
        break;
    case 1: // y-boundaries
        j = toRight ? ny-1 : 0;
        for (i=0; i<nx; i++) {
            for (k=0; k<nz; k++) {
                (*funcB)(cells[flatindex(i,j,k)]->celldata, 1);
            }
        }
        break;
    case 2: // z-boundaries
        k = toRight ? nz-1 : 0;
        for (i=0; i<nx; i++) {
            for (j=0; j<ny; j++) {
                (*funcB)(cells[flatindex(i,j,k)]->celldata, 2);
            }
        }
        break;
    }
}

// =================================================================================
// ================================ GRID REFINEMENT ================================
// =================================================================================

//! (GRID REFINEMENT)
int Tgrid::TPtrHash::HashFunction(const gridreal r[3], int iq[3]) const
{
    gptr->intcoords(r,iq);
    unsigned long x = (iq[0] << 7);
    x+= iq[1];
    x+= (iq[2] >> 7);
    return int(x % (unsigned long)(PTR_HASHTABLE_SIZE));
}

//! (GRID REFINEMENT) Return old entry or null, and add
void *Tgrid::TPtrHash::add(const gridreal r[3], void *value)
{
    int iq[3];
    const int key = HashFunction(r,iq);
    Thashnode *p;
    for (p=pool[key]; p; p=p->next)
        if (iq[0] == p->iq[0] && iq[1] == p->iq[1] && iq[2] == p->iq[2]) {
            void *const result = p->data;
            p->data = value;
            return result;
        }
    p = new Thashnode;
    p->next = pool[key];
    p->iq[0] = iq[0];
    p->iq[1] = iq[1];
    p->iq[2] = iq[2];
    p->data = value;
    pool[key] = p;
    return 0;
}

//! (GRID REFINEMENT) Same as add but check that the old value, if exists, is the same
void *Tgrid::TPtrHash::add_unique(const gridreal r[3], void *value)
{
    int iq[3];
    const int key = HashFunction(r,iq);
    Thashnode *p;
    for (p=pool[key]; p; p=p->next)
        if (iq[0] == p->iq[0] && iq[1] == p->iq[1] && iq[2] == p->iq[2]) {
            if (p->data != value) {
                errorlog << "ERROR [Tgrid::TPtrHash::add_unique]: Data is not unique\n";
                doabort();
            }
            return p->data;
        }
    p = new Thashnode;
    p->next = pool[key];
    p->iq[0] = iq[0];
    p->iq[1] = iq[1];
    p->iq[2] = iq[2];
    p->data = value;
    pool[key] = p;
    return 0;
}

//! (GRID REFINEMENT) Return entry, or null if does not exist
void *Tgrid::TPtrHash::read(const gridreal r[3]) const
{
    int iq[3];
    const int key = HashFunction(r,iq);
    Thashnode *p;
    for (p=pool[key]; p; p=p->next)
        if (iq[0] == p->iq[0] && iq[1] == p->iq[1] && iq[2] == p->iq[2])
            return p->data;
    return 0;
}

//! (GRID REFINEMENT)
void Tgrid::TPtrHash::remove(const gridreal r[3])
{
    int iq[3];
    const int key = HashFunction(r,iq);
    Thashnode *p,*q;
    bool found = false;
    for (p=pool[key]; p; p=p->next)
        if (iq[0] == p->iq[0] && iq[1] == p->iq[1] && iq[2] == p->iq[2]) {
            found = true;
            if (p == pool[key]) {
                // first entry in synonym chain
                pool[key] = p->next;
            } else {
                // not first entry
                q = pool[key];
                while (q->next != p) q = q->next;
                q->next = p->next;
            }
            delete p;
            break;
        }
    if (!found) {
        errorlog << "*** Tgrid::TPtrHash::remove" << Tr3v(r).toString() << " failed\n";
    }
}

//! (GRID REFINEMENT)
void Tgrid::TPtrHash::delete_all_as_nodeptr()
{
    int i;
    Thashnode *p;
    for (i=0; i<PTR_HASHTABLE_SIZE; i++)
        for (p=pool[i]; p; p=p->next)
            delete (Tgrid::Tnode *)p->data;
}

//! (GRID REFINEMENT)
Tgrid::TPtrHash::~TPtrHash()
{
    int i;
    Thashnode *p;
    for (i=0; i<PTR_HASHTABLE_SIZE; i++)
        while (pool[i]) {
            p = pool[i];
            pool[i] = p->next;
            delete p;
        }
}

//! (GRID REFINEMENT)
void Tgrid::update_neighbours(Tcell *c)
{
    int d,dir;
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            gridreal r[3];
            r[0] = c->centroid[0];
            r[1] = c->centroid[1];
            r[2] = c->centroid[2];
            if (dir) r[d]+= c->size;
            else r[d]-= c->size;
            c->neighbour[d][dir] = findcell_to_maxlevel(r,c->level);
        }
}

//! (GRID REFINEMENT)
void Tgrid::Tcell::take27neighbours(TCellPtr celltab[3][3][3])
{
/*
  Child numbering:
  ch=0  000  chdir[0]=0, chdir[1]=0, chdir[2]=0
  1  001  chdir[0]=0, chdir[1]=0, chdir[2]=1
  2  010  chdir[0]=0, chdir[1]=1, chdir[2]=0
  3  011  chdir[0]=0, chdir[1]=1, chdir[2]=1
  4  100  chdir[0]=1, chdir[1]=0, chdir[2]=0
  5  101  chdir[0]=1, chdir[1]=0, chdir[2]=1
  6  110  chdir[0]=1, chdir[1]=1, chdir[2]=0
  7  111  chdir[0]=1, chdir[1]=1, chdir[2]=1
  GetBit(ch,0) = chdir[2], GetBit(ch,1) = chdir[1], GetBit(ch,2) = chdir[0]
  GetBit(ch,2-d) = chdir[d]
  child[0][0][ch] = child[chdir[0]][chdir[1]][chdir[2]]
*/
    // Gather all 27 touching cells in celltab
    celltab[1][1][1] = this;
    celltab[0][1][1] = neighbour[0][0];     // -x
    celltab[2][1][1] = neighbour[0][1];     // +x
    celltab[1][0][1] = neighbour[1][0];     // -y
    celltab[1][2][1] = neighbour[1][1];     // +y
    celltab[1][1][0] = neighbour[2][0];     // -z
    celltab[1][1][2] = neighbour[2][1];     // +z
    celltab[0][0][1] = celltab[1][0][1]->neighbour[0][0];
    celltab[2][0][1] = celltab[1][0][1]->neighbour[0][1];
    celltab[0][2][1] = celltab[1][2][1]->neighbour[0][0];
    celltab[2][2][1] = celltab[2][1][1]->neighbour[1][1];
    celltab[1][0][0] = celltab[1][1][0]->neighbour[1][0];
    celltab[1][2][0] = celltab[1][1][0]->neighbour[1][1];
    celltab[1][0][2] = celltab[1][1][2]->neighbour[1][0];
    celltab[1][2][2] = celltab[1][1][2]->neighbour[1][1];
    celltab[0][1][0] = celltab[0][1][1]->neighbour[2][0];
    celltab[2][1][0] = celltab[2][1][1]->neighbour[2][0];
    celltab[0][1][2] = celltab[0][1][1]->neighbour[2][1];
    celltab[2][1][2] = celltab[2][1][1]->neighbour[2][1];
    // thus far we have 1+6+3*4 = 19 entries, 8 missing (the corners)
    celltab[0][0][0] = celltab[0][0][1]->neighbour[2][0];
    celltab[2][0][0] = celltab[1][0][0]->neighbour[0][1];
    celltab[0][2][0] = celltab[0][1][0]->neighbour[1][1];
    celltab[2][2][0] = celltab[1][2][0]->neighbour[0][1];
    celltab[0][0][2] = celltab[0][0][1]->neighbour[2][1];
    celltab[2][0][2] = celltab[2][0][1]->neighbour[2][1];
    celltab[0][2][2] = celltab[1][2][2]->neighbour[0][0];
    celltab[2][2][2] = celltab[2][2][1]->neighbour[2][1];
}

//! (GRID REFINEMENT)
void Tgrid::Tnode::update_cell_pointers(gridreal size, Tgrid& g)
{
    int dirx1,diry1,dirz1;
    gridreal r[3];
    const gridreal epsilon = size/16;   // corner node could be 3 (==dim) times refined; and half of that makes 2^16
    // epsilon is the displacement we move from the node towards each of the (dirx1,diry1,dirz1) directions
    // epsilon must be just small enough, but making it extremely small risks roundoff error.
    // there should be plenty of room however, epsilon/16 is small enough but not any smaller.
    // in case of doubt you can make it even smaller (size/32 maybe)
    for (dirx1=0; dirx1<2; dirx1++) for (diry1=0; diry1<2; diry1++) for (dirz1=0; dirz1<2; dirz1++) {
                r[0] = centroid[0] + (2*dirx1-1)*epsilon;
                r[1] = centroid[1] + (2*diry1-1)*epsilon;
                r[2] = centroid[2] + (2*dirz1-1)*epsilon;
                cell[dirx1][diry1][dirz1] = g.findcell(r);
            }
}

//! (GRID REFINEMENT) Refine a cell. Fails if any neighbour is larger.
bool Tgrid::Tcell::refine(Tgrid& g)
{
    int d,dir,a,s,ch,f1,f2,f,dirx,diry,dirz,chdir[3];
    TCellPtr c,celltab[3][3][3];        // [x][y][z]
    Tgrid::TNodePtr n,nodetab[3][3][3];
    Tgrid::TFacePtr facetab[3][2][2][3];    // [dim][+-y][+-z][left,middle,right]
    Tgrid::TPtrHash nodehash(g);
    // Check if any neighbour is larger, if yes, do nothing and return false
    const int ourlevel = level;     // avoid using bit operations but once
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) if (neighbour[d][dir]->level < ourlevel) return false;
    take27neighbours(celltab);
    // Pass through all 6 directions and all faces.
    // Add every corner node found in nodehash.
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            if (isrefined_face(d,dir)) {
                for (f1=0; f1<4; f1++) for (f2=0; f2<4; f2++)
                        nodehash.add_unique(refintf[d][dir]->face[f1]->node[f2]);
            } else {
                for (f=0; f<4; f++) nodehash.add_unique(face[d][dir]->node[f]);
            }
        }
    // Pass through all cells in celltab which may possible have a "side" node.
    // Excluded from this search are: (a) the central cell 111, (b) its direct neighbours (211,101,etc),
    // (c) corner cells (002,222,etc.). There are 12 cells in the list.
    // We define the list explicitly. Only those cells which have children OR which have refined faces are considered.
    static const int cell_searchlist[12][3] = {
        {0,0,1}, {2,0,1}, {0,2,1}, {2,2,1},
        {1,0,0}, {0,1,0}, {2,1,0}, {1,2,0},
        {1,0,2}, {0,1,2}, {2,1,2}, {1,2,2}
    };
    for (a=0; a<12; a++) {
        c = celltab[cell_searchlist[a][0]][cell_searchlist[a][1]][cell_searchlist[a][2]];
        if (!c->haschildren && !c->anyrefined_face()) continue;
        if (c->haschildren) {
            for (ch=0; ch<8; ch++) for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
                        // Refined children are not considered; they may exist but can not possibly
                        // point towards our cell.
                        // Notice that c->child[ch] might itself have children, but it doesn't harm here.
                        if (!c->child[0][0][ch]->isrefined_face(d,dir)) {
                            for (f=0; f<4; f++) {
                                nodehash.add_unique(c->child[0][0][ch]->face[d][dir]->node[f]);
                            }
                        }
                    }
        } else {
            // c has no children, but some refined face(s)
            for (d=0; d<3; d++) for (dir=0; dir<2; dir++) if (c->isrefined_face(d,dir)) {
                        for (f1=0; f1<4; f1++) for (f2=0; f2<4; f2++) {
                                nodehash.add_unique(c->refintf[d][dir]->face[f1]->node[f2]);
                            }
                    }
        }
    }
    // nodehash ready, except those that we have to create using new
    // Transfer entries from nodehash to nodetab:
    for (dirx=0; dirx<3; dirx++) for (diry=0; diry<3; diry++) for (dirz=0; dirz<3; dirz++) {
                gridreal r[3];
                r[0] = centroid[0] + 0.5*size*(dirx-1);
                r[1] = centroid[1] + 0.5*size*(diry-1);
                r[2] = centroid[2] + 0.5*size*(dirz-1);
                n = (Tgrid::Tnode *)nodehash.read(r);   // null, or the node
                if (n == 0) {
                    n = new Tnode;
                    for (d=0; d<3; d++) {
                        n->centroid[d] = r[d];
                    }
                    n->r2 = vecsqr(n->centroid);
                    nodehash.add_unique(n);
                }
                nodetab[dirx][diry][dirz] = n;
            }
    // nodehash now completely ready, it may contain spurious elements, however, which would be fatal
    // (see below) unless treated.
    // nodetab ready. Every entry is non-null, and centroid fields are initialized, but not yet cell[] pointers.
    // Initialize the face table to null
    for (a=0; a<36; a++) facetab[0][0][0][a] = 0;
    // Start filling the face table.
    // Pass through the 6 directions. If refined, take the 4 faces from existing ones,
    // if not, create them (maximum number of faces created in this phase = 6*2 = 24)
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            if (neighbour[d][dir]->haschildren) {
                // take 4 faces from existing ones
                for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) {
                        chdir[d] = !dir;
                        chdir[(d+1)%3] = diry;
                        chdir[(d+2)%3] = dirz;
                        facetab[d][diry][dirz][2*dir] = neighbour[d][dir]->child[chdir[0]][chdir[1]][chdir[2]]->face[d][!dir];
                    }
            } else {
                // create 4 new faces
                for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) {
                        facetab[d][diry][dirz][2*dir] = new Tface;
                    }
            }
        }
    // Create the 12 intra-faces. These are facetab[:][:][:][1] (3*2*2=12)
    for (d=0; d<3; d++) for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) {
                facetab[d][diry][dirz][1] = new Tface;
            }
    // facetab ready
    // Now all tables (celltab, nodetab, facetab) are ready.
    // Delete any Trefintf pieces (they would become inaccessible)
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) if (isrefined_face(d,dir)) delete refintf[d][dir];
    // Create the eight children cells and set up pointers from our cell to them.
    // Also set pointers from them to us (parent), and their cell level, and refinements status.
    for (a=0; a<8; a++) {
        chdir[0] = a/4;             // x
        chdir[1] = (a/2) % 2;       // y
        chdir[2] = a % 2;           // z
        c = new Tcell;
        child[0][0][a] = c;
        c->plist.init();
        c->haschildren = c->refine_it = c->recoarsen_it = c->forbid_psplit = false;
        c->refstatus = 0;       // no child can have a refined neighbour at this stage, this we know for sure
        c->level = ourlevel + 1;
        c->parent = this;
        c->flatind = a;
        c->running_index = -1234;       // arbitrary illegal value to ease debugging (not important though)
        c->rho_q_bg = 0.0;
        for (d=0; d<3; d++) c->centroid[d] = centroid[d] + size*(chdir[d] ? +0.25 : -0.25);
        c->r2 = vecsqr(c->centroid);
        c->size = 0.5*size;
        c->invsize = 1.0/c->size;
        c->nc = nc;
        // Averaging
        if (Params::averaging == true) {
            c->ave_nc = ave_nc;
#ifdef SAVE_POPULATION_AVERAGES
            for (int i = 0; i < Params::POPULATIONS; ++i) {
                c->pop_ave_n.push_back(pop_ave_n[i]);
                c->pop_ave_vx.push_back(pop_ave_vx[i]);
                c->pop_ave_vy.push_back(pop_ave_vy[i]);
                c->pop_ave_vz.push_back(pop_ave_vz[i]);
            }
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
            for (int i = 0; i < Params::POPULATIONS; ++i) {
                c->spectra.push_back(spectra[i]);
            }
#endif
        }
        for (s=0; s<NCELLDATA; s++) for (d=0; d<3; d++) c->celldata[s][d] = celldata[s][d];
        // set c->face pointers
        for (d=0; d<3; d++) {
            diry = chdir[(d+1)%3];
            dirz = chdir[(d+2)%3];
            for (dir=0; dir<2; dir++)
                c->face[d][dir] = facetab[d][diry][dirz][chdir[d]+dir];
        }
        // c->neighbour must be filled later
    }
    // Mark it that we now have children
    haschildren = true;
    // Set pointers from faces to nodes
    for (d=0; d<3; d++) for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) for (a=0; a<3; a++) {
                    gridreal rfacec[3],rnode[3];    // face center, node position
                    rfacec[0] = centroid[0];
                    rfacec[1] = centroid[1];
                    rfacec[2] = centroid[2];
                    rfacec[d]+= (a-1)*(0.5*size);
                    rfacec[(d+1)%3]+= (diry-0.5)*(0.5*size);
                    rfacec[(d+2)%3]+= (dirz-0.5)*(0.5*size);
                    // rfacec ready
                    static const int f_plusminus_y[4] = {-1,+1,+1,-1};  // cyclic order of nodes in face
                    static const int f_plusminus_z[4] = {-1,-1,+1,+1};
                    for (f=0; f<4; f++) {
                        rnode[0] = rfacec[0];
                        rnode[1] = rfacec[1];
                        rnode[2] = rfacec[2];
                        rnode[(d+1)%3]+= (0.25*size)*f_plusminus_y[f];
                        rnode[(d+2)%3]+= (0.25*size)*f_plusminus_z[f];
                        facetab[d][diry][dirz][a]->node[f] = (Tgrid::Tnode*)nodehash.read(rnode);
                    }
                }
    // Now everything is ready except pointers from nodes to cells, and neighbour pointers of cells.
    // Also the face/refintf fields of neighbouring cells must be updated.
    // The tree/centroid structure is ready so we can use the cell lookup functions already.
    // Set up node->cell[] pointers.
    for (dirx=0; dirx<3; dirx++) for (diry=0; diry<3; diry++) for (dirz=0; dirz<3; dirz++) {
                n = nodetab[dirx][diry][dirz];
                n->update_cell_pointers(size,g);
            }
    // Set up face/refintf fields in neighbouring cells
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            chdir[0] = chdir[1] = chdir[2] = 1;
            chdir[d] = 2*dir;
            c = celltab[chdir[0]][chdir[1]][chdir[2]];
            if (c->haschildren) {
                // neighbour c has children
                for (a=0; a<8; a++) {
                    if (GetBit(a,2-d) == dir) continue;
                    c->child[0][0][a]->face[d][!dir] = facetab[d][GetBit(a,2-(d+1)%3)][GetBit(a,2-(d+2)%3)][2*dir];
                }
            } else {
                // neighbour c is leaf cell, must create refintf piece
                // can delete the old face; it would remain an inaccessible face between two non-leaf cells
                // notice that it may have to be recreated if the grid is recoarsened
                delete c->face[d][!dir];
                c->setrefined_face(d,!dir,true);
                Tgrid::Trefintf *ref;
                c->refintf[d][!dir] = ref = new Trefintf;
                ref->face[0] = facetab[d][0][0][2*dir];
                ref->face[1] = facetab[d][1][0][2*dir];
                ref->face[2] = facetab[d][0][1][2*dir];
                ref->face[3] = facetab[d][1][1][2*dir];
            }
        }
    // Transfer everything from nodetab to nodehash2,
    // thus nodehash2 is free of the spurious elements possibly contained in nodehash.
    Tgrid::TPtrHash nodehash2(g);
    for (dirx=0; dirx<3; dirx++) for (diry=0; diry<3; diry++) for (dirz=0; dirz<3; dirz++)
                nodehash2.add_unique(nodetab[dirx][diry][dirz]);
    // Set up neighbour pointers
    for (a=0; a<8; a++) g.update_neighbours(child[0][0][a]);
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            chdir[0] = chdir[1] = chdir[2] = 1;
            chdir[d] = 2*dir;
            c = celltab[chdir[0]][chdir[1]][chdir[2]];
            if (c->haschildren) {
                // neighbour c has children, update the children's neighbours (but not c's neighbours need to be updated)
                for (a=0; a<8; a++) if (!c->child[0][0][a]->haschildren) g.update_neighbours(c->child[0][0][a]);
            }
            if (!c->haschildren) {
                int d1,dir1;
                for (d1=0; d1<3; d1++) if (d1!=d) for (dir1=0; dir1<2; dir1++) {
                            if (c->isrefined_face(d1,dir1)) continue;
                            TFacePtr modface = c->face[d1][dir1];
                            int nfound = 0;
                            for (f=0; f<4; f++) {
                                const TNodePtr n1 = modface->node[f];
                                const TNodePtr n2 = modface->node[(f+1)%4];
                                gridreal rsidenode[3];
                                int d2;
                                for (d2=0; d2<3; d2++) rsidenode[d2] = 0.5*(n1->centroid[d2] + n2->centroid[d2]);
                                TNodePtr sidenod = (Tgrid::Tnode *)nodehash2.read(rsidenode);
                                if (sidenod) {
                                    modface->sidenode[f] = sidenod;
                                    nfound++;
                                }
                            }
                        }
            }
        }
    return true;
}

//! (GRID REFINEMENT) Recoarsen a cell. Fails if the cell has no children, has grandchildren, or any neighbour of the cell's children has children.
bool Tgrid::Tcell::recoarsen(Tgrid& g)
{
    if (!haschildren) return false;
    int ch,d,dir,n,diry,dirz,chdir[3],s,f;
    for (ch=0; ch<8; ch++) {
        if (child[0][0][ch]->haschildren) return false;
        for (d=0; d<3; d++) for (dir=0; dir<2; dir++)
                if (child[0][0][ch]->neighbour[d][dir]->haschildren) return false;
    }
    // now we know that recoarsening will succeed
    TCellPtr c,celltab[3][3][3];        // [x][y][z]
    Tgrid::TPtrHash nodehash_retain(g), nodehash_remove(g);
    take27neighbours(celltab);
    // Take all 27 nodes in nodehash_remove
    // Use the 8 children and their faces. (Check that no refintf pieces are found.)
    for (ch=0; ch<8; ch++) {
        c = child[0][0][ch];
        for (d=0; d<3; d++) for (dir=0; dir<2; dir++) for (n=0; n<4; n++)
                    nodehash_remove.add_unique(c->face[d][dir]->node[n]);
    }
    // Pass through the 6 neighbours using celltab.
    // - If neighbour is not refined, check that it has a refintf towards us and delete it and its faces.
    // Create in its place a single face.
    // - If neighbour is refined, create a refintf from us to it using the existing faces
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            chdir[0] = chdir[1] = chdir[2] = 1;
            chdir[d] = 2*dir;
            c = celltab[chdir[0]][chdir[1]][chdir[2]];
            if (c->haschildren) {
                // Neighbour c is refined, create a refintf from us to it using existing faces
                setrefined_face(d,dir,true);    // mark that we are now refined wrt. that (d,dir)
                refintf[d][dir] = new Trefintf;
                int childsel[3];
                for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) {
                        childsel[d] = !dir;
                        childsel[(d+1)%3] = diry;
                        childsel[(d+2)%3] = dirz;
                        refintf[d][dir]->face[2*dirz + diry] = c->child[childsel[0]][childsel[1]][childsel[2]]->face[d][!dir];
                    }
            } else {
                // Neighbour c is not refined, check that it has a refintf towards us and delete it and its faces.
                // Create in its place a single face. Average the face quantities.
                Tgrid::TFacePtr newface = new Tface;
                for (s=0; s<NFACEDATA; s++) newface->facedata[s] = 0;
                for (f=0; f<4; f++) for (s=0; s<NFACEDATA; s++) newface->facedata[s]+= c->refintf[d][!dir]->face[f]->facedata[s];
                for (s=0; s<NFACEDATA; s++) newface->facedata[s]*= 0.25;
                for (f=0; f<4; f++) delete c->refintf[d][!dir]->face[f];
                delete c->refintf[d][!dir];
                c->setrefined_face(d,!dir,false);
                c->face[d][!dir] = newface;
                face[d][dir] = newface;
                setrefined_face(d,dir,false);
                // Set up the new face's nodes
                gridreal facecenter[3], r[3];
                facecenter[0] = centroid[0];
                facecenter[1] = centroid[1];
                facecenter[2] = centroid[2];
                facecenter[d]+= (dir-0.5)*size;
                static const int f_plusminus_y[4] = {-1,+1,+1,-1};  // cyclic order of nodes in face
                static const int f_plusminus_z[4] = {-1,-1,+1,+1};
                for (f=0; f<4; f++) {
                    r[0] = facecenter[0];
                    r[1] = facecenter[1];
                    r[2] = facecenter[2];
                    r[(d+1)%3]+= (0.5*size)*f_plusminus_y[f];
                    r[(d+2)%3]+= (0.5*size)*f_plusminus_z[f];
                    newface->node[f] = (Tgrid::Tnode *)nodehash_remove.read(r);
                }
            }
        }
    // Take nc and celldata as averages from children
    nc = 0;
    for (ch=0; ch<8; ch++) nc+= child[0][0][ch]->nc;
    nc*= 0.125;
    // Averaging
    if (Params::averaging == true) {
        ave_nc = 0;
        for (ch=0; ch<8; ch++) {
            ave_nc+= child[0][0][ch]->ave_nc;
        }
        ave_nc*= 0.125;
#ifdef SAVE_POPULATION_AVERAGES
        for (int i = 0; i < Params::POPULATIONS; ++i) {
            pop_ave_n[i] = 0.0;
            pop_ave_vx[i] = 0.0;
            pop_ave_vy[i] = 0.0;
            pop_ave_vz[i] = 0.0;
            for (ch=0; ch<8; ch++) {
                pop_ave_n[i] += child[0][0][ch]->pop_ave_n[i];
                pop_ave_vx[i] += child[0][0][ch]->pop_ave_vx[i];
                pop_ave_vy[i] += child[0][0][ch]->pop_ave_vy[i];
                pop_ave_vz[i] += child[0][0][ch]->pop_ave_vz[i];
            }
            pop_ave_n[i] *= 0.125;
            pop_ave_vx[i] *= 0.125;
            pop_ave_vy[i] *= 0.125;
            pop_ave_vz[i] *= 0.125;
        }
#endif
    }
    for (s=0; s<NCELLDATA; s++) for (d=0; d<3; d++) celldata[s][d] = 0;
    for (ch=0; ch<8; ch++) for (s=0; s<NCELLDATA; s++) for (d=0; d<3; d++) celldata[s][d]+= child[0][0][ch]->celldata[s][d];
    for (s=0; s<NCELLDATA; s++) for (d=0; d<3; d++) celldata[s][d]*= 0.125;
    // Delete intra-cell faces (12)
    for (d=0; d<3; d++) for (ch=0; ch<8; ch++) {
            if (GetBit(ch,2-d)) continue;
            delete child[0][0][ch]->face[d][1];
        }
    // Delete children cells
    for (ch=0; ch<8; ch++) {
        delete child[0][0][ch];
    }
    haschildren = false;
    // Go through the 6 neighbours, excluding their outer faces.
    // Also exclude leaf cells. Remove all accessible nodes from nodehash_remove,
    // but before it update their cell pointers (update_cell_pointers()).
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            chdir[0] = chdir[1] = chdir[2] = 1;
            chdir[d] = 2*dir;
            c = celltab[chdir[0]][chdir[1]][chdir[2]];
            if (!c->haschildren) continue;  // exclude leaf cells
            for (ch=0; ch<8; ch++) {
                if (GetBit(ch,2-d) == dir) continue;    // exclude outer row of cells
                int d1,dir1;
                for (d1=0; d1<3; d1++) for (dir1=0; dir1<2; dir1++) {
                        if (d1==d && dir1==dir) continue;   // exclude outer faces
                        Tgrid::TFacePtr face1 = c->child[0][0][ch]->face[d1][dir1];
                        for (f=0; f<4; f++) {
                            const Tgrid::TNodePtr nod = face1->node[f];
                            if (nodehash_remove.read(nod->centroid)) {
                                nod->update_cell_pointers(size,g);
                                nodehash_remove.remove(nod->centroid);
                                nodehash_retain.add_unique(nod);
                            }
                        }
                    }
            }
        }
    // Similarly, pass through the faces of *this (only non-refined faces)
    // to get hand to the 8 corner nodes (for refined faces, the nodes have already
    // been processed by the previous loop).
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            if (isrefined_face(d,dir)) continue;
            Tgrid::TFacePtr face1 = face[d][dir];
            for (f=0; f<4; f++) {
                const Tgrid::TNodePtr nod = face1->node[f];
                if (nodehash_remove.read(nod->centroid)) {
                    nod->update_cell_pointers(size,g);
                    nodehash_remove.remove(nod->centroid);
                    nodehash_retain.add_unique(nod);
                }
            }
        }
    // Update neighbour pointers of neighbour's children
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) {
            chdir[0] = chdir[1] = chdir[2] = 1;
            chdir[d] = 2*dir;
            c = celltab[chdir[0]][chdir[1]][chdir[2]];
            if (c->haschildren) {
                for (ch=0; ch<8; ch++) g.update_neighbours(c->child[0][0][ch]);
            }
            // Remove those sidenodes which are in nodehash_remove (because those nodes will be soon deleted)
            if (!c->haschildren) {
                int d1,dir1;
                for (d1=0; d1<3; d1++) for (dir1=0; dir1<2; dir1++)
                        if (!c->isrefined_face(d1,dir1)) {
                            for (f=0; f<4; f++) {
                                Tgrid::TNodePtr& sidenod = c->face[d1][dir1]->sidenode[f];
                                if (sidenod == 0) continue;
                                if (nodehash_remove.read(sidenod->centroid)) sidenod = 0;   // note that sidenod is ref. variable
                            }
                        }
            }
        }
    // Pass through potential side nodes of the faces of *this (4*6 = 24).
    // Add those sidenodes which are in nodehash_retain.
    for (d=0; d<3; d++) for (dir=0; dir<2; dir++) if (!isrefined_face(d,dir)) {
                for (f=0; f<4; f++) {
                    const Tgrid::TNodePtr n1 = face[d][dir]->node[f];
                    const Tgrid::TNodePtr n2 = face[d][dir]->node[(f+1)%4];
                    gridreal r[3];
                    int d1;
                    for (d1=0; d1<3; d1++) r[d1] = 0.5*(n1->centroid[d1] + n2->centroid[d1]);
                    Tgrid::TNodePtr sidenod = (Tgrid::Tnode *)nodehash_retain.read(r);
                    if (sidenod) face[d][dir]->sidenode[f] = sidenod;
                }
            }
    // Delete those nodes that remain in nodehash_remove
    nodehash_remove.delete_all_as_nodeptr();
    return true;
}

//! (GRID REFINEMENT) Returns the number of cells marked
int Tgrid::Tcell::mark_refinement_recursive(GridRefinementProfile refFunc)
{
    int result = 0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) result+= child[0][0][ch]->mark_refinement_recursive(refFunc);
    } else {
        if (size > refFunc.getValue(centroid)) {
            refine_it = true;
            result++;
        }
    }
    return result;
}

//! (GRID REFINEMENT) Returns the number of cells marked (recursive)
int Tgrid::Tcell::mark_recoarsening_recursive(gridreal (*mindx)(const gridreal[]))
{
    /*
     * Set recoarsen_it bit for ALL children cells when needed
     * Return value is the number of PARENT cells to be recoarsened
     * (that is, the total number of recoarsen_it bits set is
     * 8 times higher than the return value)
     */
    int result = 0;
    if (haschildren) {
        int ch;
        bool has_grandchildren = false;
        for (ch=0; ch<8; ch++) if (child[0][0][ch]->haschildren) {
                has_grandchildren = true;
                break;
            }
        if (has_grandchildren) {
            for (ch=0; ch<8; ch++) result+= child[0][0][ch]->mark_recoarsening_recursive(mindx);
        } else {
            // if all children must be recoarsened, do so, otherwise nothing
            bool do_all = true;
            for (ch=0; ch<8; ch++) if (child[0][0][ch]->size >= (*mindx)(child[0][0][ch]->centroid)) {
                    do_all = false;
                    break;
                }
            if (do_all) {
                for (ch=0; ch<8; ch++) child[0][0][ch]->recoarsen_it = true;
                result = 1;
            }
        }
    }
    return result;
}

//! (GRID REFINEMENT) Call refine for nonleaf cells (recursive)
void Tgrid::Tcell::refine_recursive(Tgrid& g)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->refine_recursive(g);
    } else {
        if (refine_it) {
            if (!refine(g)) {
                errorlog << "warning: refine() failed at " << Tr3v(centroid).toString() << "\n";
            }
            refine_it = false;
        }
    }
}

//! (GRID REFINEMENT) Call recoarsen for nonleaf cells (recursive)
void Tgrid::Tcell::recoarsen_recursive(Tgrid& g)
{
    if (haschildren) {
        int ch;
        bool has_grandchildren = false;
        for (ch=0; ch<8; ch++) if (child[0][0][ch]->haschildren) {
                has_grandchildren = true;
                break;
            }
        if (!has_grandchildren && child[0][0][0]->recoarsen_it) {
            if (!recoarsen(g)) {
                errorlog << "warning: recoarsen() failed at " << Tr3v(centroid).toString() << "\n";
                for (ch=0; ch<8; ch++) child[0][0][ch]->recoarsen_it = false;
            }
        } else {
            for (ch=0; ch<8; ch++) child[0][0][ch]->recoarsen_recursive(g);
        }
    }
}

//! (GRID REFINEMENT) Refine grid
void Tgrid::Refine(GridRefinementProfile refFunc)
{
    MSGFUNCTIONCALL("Tgrid::Refine");
    mainlog << "Refining the grid:\n";
    Params::currentGridRefinementLevel = 0;
    // Do refinements
    vector<int> levelNRefinedCells;
    for (int level=0; level<Params::maxGridRefinementLevel; level++) {
        mainlog << " level " << level << ": ";
        int nmarked = 0;
        int i,j,k;
        ForInterior(i,j,k) {
            nmarked+= cells[flatindex(i,j,k)]->mark_refinement_recursive(refFunc);
        }
        mainlog << "refining " << nmarked << " cells\n";
        if (nmarked == 0) {
            break;
        }
        ForInterior(i,j,k) {
            cells[flatindex(i,j,k)]->refine_recursive(*this);
        }
        levelNRefinedCells.push_back(nmarked);
        ++Params::currentGridRefinementLevel;
    }
    levelNRefinedCells.push_back(0);
    vector<int> levelNCells;
    // Zero (base) level cells
    levelNCells.push_back(nx*ny*nz - levelNRefinedCells[0]);
    // Nth level cells
    for(unsigned int ii=1; ii < levelNRefinedCells.size(); ++ii) {
        levelNCells.push_back(levelNRefinedCells[ii-1]*8 - levelNRefinedCells[ii]);
    }
    // Report cells
    mainlog << "|--------------------- GRID CELL REFINEMENTS ---------------------|\n";
    int totalCells = 0;
    for(unsigned int ii=0; ii < levelNCells.size(); ++ii) {
        mainlog << "| Cells in level " << static_cast<int>(ii) << " : " << levelNCells[ii];
        if(ii == 0) {
            mainlog << " (base cells including ghosts)";
        }
        mainlog<< "\n";
        totalCells += levelNCells[ii];
    }
    int totalCellsWithoutGhosts = totalCells -nx*ny*nz + (nx-2)*(ny-2)*(nz-2);
    mainlog << "| Total cells in all levels                     : " << Ncells_with_ghosts() << "\n";
    mainlog << "| Total cells in all levels (no ghosts)         : " << Ncells_without_ghosts() << "\n";
    mainlog << "| Total cells in all levels (no parents)        : " << totalCells << "\n";
    mainlog << "| Total cells in all levels (no parents/ghosts) : " << totalCellsWithoutGhosts << "\n";
    mainlog << "| Total macroparticles in the box (average)     : " << totalCellsWithoutGhosts*Params::macroParticlesPerCell << "\n";
    mainlog << "|-----------------------------------------------------------------|\n";
    // reset the cached cell pointer since it may have been invalidated
    previous_found_cell = 0;
    MSGFUNCTIONEND("Tgrid::Refine");
}

//! (GRID REFINEMENT) Recoarsen grid
void Tgrid::recoarsen(gridreal (*mindx)(const gridreal[]))
{
    int i,j,k;
    int level;
    for (level=0; level<10; level++) {
        int nmarked = 0;
        ForInterior(i,j,k)
        nmarked+= cells[flatindex(i,j,k)]->mark_recoarsening_recursive(mindx);
        if (nmarked == 0) break;
        mainlog << "Tgrid::recoarsen: recoarsening = " << nmarked << " cells\n";
        ForInterior(i,j,k)
        cells[flatindex(i,j,k)]->recoarsen_recursive(*this);
    }
    previous_found_cell = 0;        // reset the cached cell pointer since it may have been invalidated
}

// =================================================================================
// ================================ HC-FILE WRITING ================================
// =================================================================================

static int n_int_bytes; //!< must be set before calling WriteInt (is set by Tgrid::hcwrite_MHD)

//! Write integer
inline void WriteInt(ostream& o, int x)
{
    int i;
    for (i=0; i<n_int_bytes; i++)
        o.put((unsigned char)((x >> 8*i) & 0xFF));
}

//! Enunerate child cells
void Tgrid::Tcell::enum_children_recursive()
{
    if (haschildren) {
        int dirx,diry,dirz;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->running_index = Tgrid::cell_running_index++;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->enum_children_recursive();
    }
}

//! Calculate fluid parameters in a cell
void Tgrid::Tcell::cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, vector<int> popId)
{
    n = plist.calc_mass(popId)/(Params::pops[popId[0]]->m*size*size*size);
    vx = vy = vz = 0;
    plist.calc_avev(vx, vy, vz, popId);
    const real T = 0.5*plist.calc_avemv2(vx, vy, vz, popId);
    P = n*T;
}

//! Write population and plasma quantities (recursive)
void Tgrid::Tcell::writeMHD_children_recursive(ostream& o,const int filetype,vector<int> popId) const
{
    if (haschildren) {
        int dirx,diry,dirz;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeMHD(o,filetype, popId);
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeMHD_children_recursive(o,filetype,popId);
    }
}

//! Write population and plasma quantities
void Tgrid::Tcell::writeMHD(ostream& o,const int filetype,vector<int> popId) const
{
    if (filetype == 1 || filetype == 2) {
        popId.clear();
    }
    if (haschildren) {
        if (hcFileAsciiFormat == false) {
            o.put('N');
            WriteInt(o,parent ? parent->running_index : -1);
            WriteInt(o,child[0][0][0]->running_index);
        } else {
            o << "N " << (parent ? parent->running_index : -1) << ' ' << child[0][0][0]->running_index << '\n';
        }
    } else {
        if (hcFileAsciiFormat == false) {
            o.put('L');
            WriteInt(o,parent ? parent->running_index : -1);
        } else {
            o << "L " << (parent ? parent->running_index : -1) << ' ';
        }
        real rho = 0;
        real rhovx = 0;
        real rhovy = 0;
        real rhovz = 0;
        real U1 = 0;
        real B1x = 0;
        real B1y = 0;
        real B1z = 0;
        real B0x = 0;
        real B0y = 0;
        real B0z = 0;
        real n=0,vx=0,vy=0,vz=0;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        gridreal dV = size*size*size;
#else
        gridreal dV = sph_dV;
        if (dV<0.0) dV = -dV;  // This because r = x and x can be (and usually is!) negative
#endif
        if(filetype == 0) { // Population(s)
            // This gives correct mass density, but hcvis may assume that
            // m = mp and then, e.g., for O+ populations n comes out wrong
            // in hcvis. The volume weighting (accumulation) is not used!
            rho = plist.calc_mass(popId)/dV;
            // This gives correct number density if all the populations
            // to be saved in this hc-file have the same particle mass
            // (as should be usually the case)!
            n = rho/Params::pops[popId[0]]->m;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
            // This is B after the propagations in this timestep
            B1x = 0.5*(faceave(0,0,FACEDATA_B) + faceave(0,1,FACEDATA_B));
            B1y = 0.5*(faceave(1,0,FACEDATA_B) + faceave(1,1,FACEDATA_B));
            B1z = 0.5*(faceave(2,0,FACEDATA_B) + faceave(2,1,FACEDATA_B));
#else
            B1x = celldata[CELLDATA_B][0];
            B1y = celldata[CELLDATA_B][1];
            B1z = celldata[CELLDATA_B][2];
#endif
            // Average velocity among the selected populations in the cell.
            // The volume weighting (accumulation) is not used!
            plist.calc_avev(vx,vy,vz,popId);
        } else if(filetype == 1) { // Average
            // Use proton mass here => rho is generally incorrect, but
            // hcvis can display correct number density for temporal
            // average files. Accumulated value!
            if (Params::bg_in_avehcfile==false) {
                rho = Params::m_p * ave_nc;
            } else {
                rho = Params::m_p * (ave_nc + rho_q_bg/Params::e);
            }
            // This is correct temporal average number density. It has little
            // use, since temporal average files cannot have correct U. To get
            // correct U, something like CELLDATA_AVE_VQ, CELLDATA_AVE_T etc.
            // would be needed. Accumulated value!
            n = ave_nc;
            // B is correct in temporal average hc-file
            B1x = 0.5*(faceave(0,0,FACEDATA_AVEB) + faceave(0,1,FACEDATA_AVEB));
            B1y = 0.5*(faceave(1,0,FACEDATA_AVEB) + faceave(1,1,FACEDATA_AVEB));
            B1z = 0.5*(faceave(2,0,FACEDATA_AVEB) + faceave(2,1,FACEDATA_AVEB));
            // At least n (remember to use proton mass for mp in hcvis)
            // and B come out correct in hcvis when displaying temporal
            // average hc-files.
            vx = celldata[CELLDATA_Ji][0]/rho_q;
            vy = celldata[CELLDATA_Ji][1]/rho_q;
            vz = celldata[CELLDATA_Ji][2]/rho_q;
        } else if(filetype == 2) { // Plasma
            // Total mass density in the cell.
            // The volume weighting (accumulation) is not used!
            rho = plist.calc_mass()/dV;
            // Total particle number density in the cell.
            // The volume weighting (accumulation) is not used!
            // Accumulated nc could be used but then rho would
            // not consistent with this value.
            n = plist.calc_weight()/dV;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
            B1x = 0.5*(faceave(0,0,FACEDATA_B) + faceave(0,1,FACEDATA_B));
            B1y = 0.5*(faceave(1,0,FACEDATA_B) + faceave(1,1,FACEDATA_B));
            B1z = 0.5*(faceave(2,0,FACEDATA_B) + faceave(2,1,FACEDATA_B));
#else
            gridreal dS0[3] = {sph_dS[0], sph_dS[1], sph_dS[2]};
            gridreal dS1[3] = {sph_dS_next[0], sph_dS_next[1], sph_dS_next[2]};
            gridreal dSc[3] = {sph_dS_centroid[0], sph_dS_centroid[1], sph_dS_centroid[2]};
            B1x = 0.5*(faceave(0,0,FACEDATA_B)*dS0[0] + faceave(0,1,FACEDATA_B)*dS1[0])/dSc[0];
            B1y = 0.5*(faceave(1,0,FACEDATA_B)*dS0[1] + faceave(1,1,FACEDATA_B)*dS1[1])/dSc[1];
            B1z = 0.5*(faceave(2,0,FACEDATA_B)*dS0[2] + faceave(2,1,FACEDATA_B)*dS1[2])/dSc[2];
#endif
            // Average (MHD single fluid) velocity in the cell.
            // The volume weighting (accumulation) is not used!
            plist.calc_U(vx,vy,vz);
        }
#ifdef SAVE_POPULATION_AVERAGES
        else if(filetype == 3) { // Average population(s)
            for(unsigned int i=0; i< popId.size(); i++) {
                int j = popId[i];
                n += pop_ave_n[j];
                vx += pop_ave_n[j]*pop_ave_vx[j];
                vy += pop_ave_n[j]*pop_ave_vy[j];
                vz += pop_ave_n[j]*pop_ave_vz[j];
            }
            real inv_n=0;
            if(n>0) {
                inv_n = 1.0/n;
            }
            vx *= inv_n;
            vy *= inv_n;
            vz *= inv_n;
            rho = n*Params::pops[popId[0]]->m;
            B1x = 0.5*(faceave(0,0,FACEDATA_AVEB) + faceave(0,1,FACEDATA_AVEB));
            B1y = 0.5*(faceave(1,0,FACEDATA_AVEB) + faceave(1,1,FACEDATA_AVEB));
            B1z = 0.5*(faceave(2,0,FACEDATA_AVEB) + faceave(2,1,FACEDATA_AVEB));
        }
#endif
        else {
            ERRORMSG("internal error");
        }
        // Particle momentum density
        rhovx = rho*vx;
        rhovy = rho*vy;
        rhovz = rho*vz;
        const real T = 0.5*plist.calc_avemv2(vx,vy,vz,popId);
        const real P = n*T;
        // Total energy without constant magnetic field
        U1 = P/(Params::gamma-1) + 0.5*rho*( sqr(vx) + sqr(vy) + sqr(vz) ) + ( sqr(B1x) + sqr(B1y) + sqr(B1z) )/(2*Params::mu_0);
        // Constant magnetic field
        real B0[3] = {0.0, 0.0, 0.0};
        fastreal r[3] = {centroid[0], centroid[1], centroid[2]};
        addConstantMagneticField(r, B0);
        B0x = B0[0];
        B0y = B0[1];
        B0z = B0[2];
        if (hcFileAsciiFormat == false) {
            const int n = 11;
            float xf[n];
            xf[0] = rho;
            xf[1] = rhovx;
            xf[2] = rhovy;
            xf[3] = rhovz;
            xf[4] = U1;
            xf[5] = B1x;
            xf[6] = B1y;
            xf[7] = B1z;
            xf[8] = B0x;
            xf[9] = B0y;
            xf[10] = B0z;
            ByteConversion(sizeof(float),(unsigned char*)xf,n);
            WriteFloatsToFile(o,xf,n);
        } else {
            o << rho << ' '
              << rho*vx << ' '
              << rho*vy << ' '
              << rho*vz << ' '
              << U1 << ' '
              << B1x << ' '
              << B1y << ' '
              << B1z << ' '
              << B0x << ' '
              << B0y << ' '
              << B0z << ' ';
        }
        if (hcFileAsciiFormat == true) {
            o << '\n';
        }
    }
}

//! Write debug quantities (recursive)
void Tgrid::Tcell::writeDBUG_children_recursive(ostream& o) const
{
    if (haschildren) {
        int dirx,diry,dirz;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeDBUG(o);
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeDBUG_children_recursive(o);
    }
}

//! Write debug quantities
void Tgrid::Tcell::writeDBUG(ostream& o) const
{
    if (haschildren) {
        if (hcFileAsciiFormat == false) {
            o.put('N');
            WriteInt(o,parent ? parent->running_index : -1);
            WriteInt(o,child[0][0][0]->running_index);
        } else {
            o << "N " << (parent ? parent->running_index : -1) << ' ' << child[0][0][0]->running_index << '\n';
        }
    } else {
        if (hcFileAsciiFormat == false) {
            o.put('L');
            WriteInt(o,parent ? parent->running_index : -1);
        } else {
            o << "L " << (parent ? parent->running_index : -1) << ' ';
        }
        real rho = plist.Nparticles();
        real rhovx = celldata[CELLDATA_J][0];
        real rhovy = celldata[CELLDATA_J][1];
        real rhovz = celldata[CELLDATA_J][2];
        real U1 = rho_q;
        real Btotx,Btoty,Btotz,Ex,Ey,Ez;
        // Magnetic field
        real B0[3] = {0.0,0.0,0.0};
        addConstantMagneticField(centroid,B0);
        Btotx = celldata[CELLDATA_B][0] + B0[0];
        Btoty = celldata[CELLDATA_B][1] + B0[1];
        Btotz = celldata[CELLDATA_B][2] + B0[2];
        // Electric field Ue x B (+grad pe)
        Ex = Btoty*celldata[CELLDATA_UE][2] - Btotz*celldata[CELLDATA_UE][1];
        Ey = Btotz*celldata[CELLDATA_UE][0] - Btotx*celldata[CELLDATA_UE][2];
        Ez = Btotx*celldata[CELLDATA_UE][1] - Btoty*celldata[CELLDATA_UE][0];
        if(Params::electronPressure==true) {
            Ex += celldata[CELLDATA_TEMP2][0];
            Ey += celldata[CELLDATA_TEMP2][1];
            Ez += celldata[CELLDATA_TEMP2][2];
        }
        if (hcFileAsciiFormat == false) {
            const int n = 11;
            float xf[n];
            xf[0] = rho;
            xf[1] = rhovx;
            xf[2] = rhovy;
            xf[3] = rhovz;
            xf[4] = U1;
            xf[5] = Btotx;
            xf[6] = Btoty;
            xf[7] = Btotz;
            xf[8] = Ex;
            xf[9] = Ey;
            xf[10] = Ez;
            ByteConversion(sizeof(float),(unsigned char*)xf,n);
            WriteFloatsToFile(o,xf,n);
        } else {
            o << rho << ' '
              << rhovx << ' '
              << rhovy << ' '
              << rhovz << ' '
              << U1 << ' '
              << Btotx << ' '
              << Btoty << ' '
              << Btotz << ' '
              << Ex << ' '
              << Ey << ' '
              << Ez;
        }
        if (hcFileAsciiFormat == true) {
            o << '\n';
        }
    }
}

//! Write extra quantities (recursive)
void Tgrid::Tcell::writeEXTRA_children_recursive(ostream& o, ScalarField* s) const
{
    if (haschildren) {
        int dirx,diry,dirz;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeEXTRA(o,s);
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeEXTRA_children_recursive(o,s);
    }
}

//! Write extra quantities
void Tgrid::Tcell::writeEXTRA(ostream& o, ScalarField* s) const
{
    if (haschildren) {
        if (hcFileAsciiFormat == false) {
            o.put('N');
            WriteInt(o,parent ? parent->running_index : -1);
            WriteInt(o,child[0][0][0]->running_index);
        } else {
            o << "N " << (parent ? parent->running_index : -1) << ' ' << child[0][0][0]->running_index << '\n';
        }
    } else {
        if (hcFileAsciiFormat == false) {
            o.put('L');
            WriteInt(o,parent ? parent->running_index : -1);
        } else {
            o << "L " << (parent ? parent->running_index : -1) << ' ';
        }
        gridreal r[3] = {centroid[0], centroid[1], centroid[2]};
        //real rho = (*func)(r,id);
        real rho = s->getValue(r);
        if (hcFileAsciiFormat == false) {
            const int n = 1;
            float xf[n];
            xf[0] = rho;
            ByteConversion(sizeof(float),(unsigned char*)xf,n);
            WriteFloatsToFile(o,xf,n);
        } else {
            o << rho;
        }
        if (hcFileAsciiFormat == true) {
            o << '\n';
        }
    }
}

#ifdef SAVE_PARTICLE_CELL_SPECTRA

//! Write spectra (recursive)
void Tgrid::Tcell::writeSPECTRA_children_recursive(ostream& o,vector<int> popId) const
{
    if (haschildren) {
        int dirx,diry,dirz;
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeSPECTRA(o,popId);
        for (dirz=0; dirz<2; dirz++) for (diry=0; diry<2; diry++) for (dirx=0; dirx<2; dirx++)
                    child[dirx][diry][dirz]->writeSPECTRA_children_recursive(o,popId);
    }
}

//! Write spectra
void Tgrid::Tcell::writeSPECTRA(ostream& o,vector<int> popId) const
{
    if (haschildren) {
        if (hcFileAsciiFormat == false) {
            o.put('N');
            WriteInt(o,parent ? parent->running_index : -1);
            WriteInt(o,child[0][0][0]->running_index);
        } else {
            o << "N " << (parent ? parent->running_index : -1) << ' ' << child[0][0][0]->running_index << '\n';
        }
    } else {
        if (hcFileAsciiFormat == false) {
            o.put('L');
            WriteInt(o,parent ? parent->running_index : -1);
        } else {
            o << "L " << (parent ? parent->running_index : -1) << ' ';
        }
        if (hcFileAsciiFormat == false) {
            const int n = Params::spectraNbins;
            float xf[n];
            for(int i=0; i<n; i++) {
                xf[i] = 0.0;
                for(unsigned int j=0; j<popId.size(); ++j) {
                    xf[i] += spectra[popId[j]][i];
                }
            }
            ByteConversion(sizeof(float),(unsigned char*)xf,n);
            WriteFloatsToFile(o,xf,n);
        } else {
            const int n = Params::spectraNbins;
            for(int i=0; i<n; i++) {
                real sum = 0.0;
                for(unsigned int j=0; j<popId.size(); ++j) {
                    sum += spectra[popId[j]][i];
                }
                o << sum << ' ';
            }
        }
        if (hcFileAsciiFormat == true) {
            o << '\n';
        }
    }
}

#endif

//! Encode cell info
inline unsigned char EncodeCellInfo(int celltype)
{
    unsigned char ch;
    ch = (unsigned char)celltype;
    if (true /*may_subdivide(i)*/) ch|= 0x80;
    if (true /*may_recoarsen(i)*/) ch|= 0x40;
    return ch;
}

/** \brief Number of bytes
 *
 * Take care that e.g. NumberOfBytes(50000) is 3, not 2, because the sign bit must
 * also be included. That's why we compute (nbits+8)/8, otherwise it would be
 * (nbits+7)/8.
 */
static int NumberOfBytes(int x)
{
    if (x < 0) x = -x;
    const int maxnbits = 8*sizeof(x);
    int i,nbits=0;
    for (i=0; i<maxnbits; i++)
        if (((x >> i) & 0x1) != 0) nbits = i+1;
    return ((nbits+8)/8);
}

#ifdef SAVE_PARTICLES_ALONG_ORBIT

static ofstream ofple; //!< File where to save particle along orbit
static int ple_flatint_temp; //!< flatindex of a cell particle is located in

//! Flag cells along orbit given in file fn in which particles are saved
void Tgrid::set_save_particles_orbit(const char *fn)
{
    vector< vector<real> > orbit;
    orbit = readRealsFromFile(fn);
    mainlog << "SAVING PARTICLES ALONG SPACECRAFT ORBIT: " << fn << "\n";
    if(orbit.size() <= 0) {
        ERRORMSG2("cannot read orbit file",string(fn));
        doabort();
    }
    ofstream of;
    of.open("particles_along_orbit_cellindices.dat");
    of << scientific;
    of.precision(5);
    for(unsigned int i = 0; i < orbit.size(); ++i) {
        if(orbit[i].size() != 3) {
            ERRORMSG("Bad orbit file line structure, skipping...");
            continue;
        }
        shortreal r[3] = {orbit[i][0],orbit[i][1],orbit[i][2]};
        gridreal lowercorner[3];
        Tcell *const c = findcell(r,lowercorner);
        if(!c) {
            continue;
        }
        c->save_particles = true;
        of << c->flatind << " " << c->centroid[0] << " " << c->centroid[1] << " " << c->centroid[2] << endl;
    }
    of.close();
}

//! Write a particle in a file
bool write_ple(TLinkedParticle& p)
{
    ofple
            << Params::t << " " << ple_flatint_temp << " " << p.popid << " " << p.w << " "
            << p.x << " " << p.y << " " << p.z << " "
            << p.vx << " " << p.vy << " " << p.vz << endl;
    ofple << flush;
    return true;
}

//! Write particles along orbit
void Tgrid::particles_write()
{
    static bool init_done =false;
    if(init_done == false) {
        ofple.open("particles_along_orbit.dat");
        ofple << scientific;
        ofple.precision(3);
        init_done = true;
    }
    for (int i=0; i<nx-1; i++) for (int j=0; j<ny-1; j++) for (int k=0; k<nz-1; k++) {
                int c = flatindex(i,j,k);
                cells[c]->particles_write_recursive();
            }
}

//! Write particles along orbit (recursive)
void Tgrid::Tcell::particles_write_recursive()
{
    if(haschildren) {
        for (int ch=0; ch<8; ch++) {
            child[0][0][ch]->particles_write_recursive();
        }
    } else {
        if(save_particles == false) {
            return;
        }
        ple_flatint_temp = flatind;
        plist.pass(write_ple);
    }
}

#endif

//! Write population and plasma hc-file
bool Tgrid::hcwrite_MHD(const char *fn,string ascbin,string hctype,vector<int> popId)
{
    if(ascbin.compare("binary") == 0) {
        hcFileAsciiFormat = false;
    } else if(ascbin.compare("ascii") == 0) {
        hcFileAsciiFormat = true;
    }
    // File types: 0 = populations, 1 = average, 2 = plasma
    int filetype;
    if(hctype.compare("populations") == 0) {
        filetype = 0;
    } else if(hctype.compare("average") == 0) {
        filetype = 1;
    } else if(hctype.compare("plasma") == 0) {
        filetype = 2;
    }
#ifdef SAVE_POPULATION_AVERAGES
    else if(hctype.compare("populations_ave") == 0) {
        filetype = 3;
    }
#endif
    else {
        ERRORMSG("hc-file type not recognized");
        filetype = - 1;
        return false;
    }
    if(filetype == 0 && popId.size() <= 0) {
        ERRORMSG("no population ids");
        return false;
    }
    if(filetype == 1 && Params::averaging == false) {
        ERRORMSG("averaging turned off");
        return false;
    }
    if(filetype == 1 || filetype == 2) {
        popId.clear();
    }
#ifdef SAVE_POPULATION_AVERAGES
    if(filetype == 3 && Params::averaging == 0) {
        ERRORMSG("no averaging used");
        return false;
    }
#endif
    const int ncells = Ncells_with_ghosts();
    ofstream o(fn);
    if (!o.good()) {
        return false;
    }
    const int oldprec = o.precision();
    o.precision(16);
    int c;
    const int Nbase = nx*ny*nz;
    // Enumerate cells
    for (c=0; c<Nbase; c++) {
        cells[c]->running_index = c;
    }
    cell_running_index = Nbase;
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->enum_children_recursive();
    }
    n_int_bytes = NumberOfBytes(ncells);
    o << "# " << Params::codeVersion << "\n";
    o << "# filename: " << fn << "\n";
    o << "# t = " << Params::t << "\n";
    o << "# gamma = " << Params::gamma << "\n";
    o << "# mu0 = " << Params::mu_0 << "\n";
    o << "# R = " << Params::R_P << "\n";
    if (filetype == 0) {
        o << "# POPULATION FILE CONTENTS\n";
        o << "#  rho = mass density (NGP)\n";
        o << "#  rhovx = momentum density x-comp (NGP)\n";
        o << "#  rhovy = momentum density y-comp (NGP)\n";
        o << "#  rhovz = momentum density z-comp (NGP)\n";
        o << "#  U1 = total energy density without constant magnetic fields (NGP)\n";
        o << "#  B1x = interaction magnetic field x-comp\n";
        o << "#  B1y = interaction magnetic field y-comp\n";
        o << "#  B1z = interaction magnetic field z-comp\n";
        o << "#  B0x = constant magnetic field x-comp\n";
        o << "#  B0y = constant magnetic field y-comp\n";
        o << "#  B0z = constant magnetic field z-comp\n";
        o << "# POPULATION FILE CONTENTS\n";
        o << "# populations = ";
        // Check particle masses and charges
        real mm = Params::pops[popId[0]]->m;
        real qq = Params::pops[popId[0]]->q;
        for (unsigned int i = 0; i < popId.size(); i++) {
            o << Params::pops[popId[i]]->getIdStr() << " ";
            if (Params::pops[popId[i]]->m != mm) {
                o << "(WARNING: particle mass) ";
                errorlog << "WARNING [Tgrid::hcwrite_MHD]: particle mass not same for all populations in a hc-file (" << fn << ")\n";
            }
            if (Params::pops[popId[i]]->q != qq) {
                o << "(WARNING: particle charge) ";
                errorlog << "WARNING [Tgrid::hcwrite_MHD]: particle charge not same for all populations in a hc-file (" << fn << ")\n";
            }
        }
        o << "\n";
        o << "# m = " << Params::pops[popId[0]]->m << "\n";
        o << "# q = " << Params::pops[popId[0]]->q << "\n";
    } else if(filetype == 1) {
        o << "# AVERAGE FILE CONTENTS\n";
        if (Params::bg_in_avehcfile == 1) {
            o << "#  rho = (average_number_density+bgdensity)*mp (temporal average, PIC) \n";
        } else {
            o << "#  rho = average_number_density*mp (temporal average, PIC) \n";
        }
        o << "#  rhovx = momentum density x-comp (mixed temp.ave./instantaneous, PIC)\n";
        o << "#  rhovy = momentum density y-comp (mixed temp.ave./instantaneous, PIC)\n";
        o << "#  rhovz = momentum density z-comp (mixed temp.ave./instantaneous, PIC)\n";
        o << "#  U1 = total energy density without constant magnetic fields (mixed, PIC/NGP)\n";
        o << "#  B1x = interaction magnetic field x-comp (temporal average)\n";
        o << "#  B1y = interaction magnetic field y-comp (temporal average)\n";
        o << "#  B1z = interaction magnetic field z-comp (temporal average)\n";
        o << "#  B0x = constant magnetic field x-comp\n";
        o << "#  B0y = constant magnetic field y-comp\n";
        o << "#  B0z = constant magnetic field z-comp\n";
        o << "# AVERAGE FILE CONTENTS\n";
        o << "# populations = all (average file)\n";
        o << "# m = average\n";
        o << "# q = average\n";
    } else if(filetype == 2) {
        o << "# PLASMA FILE CONTENTS\n";
        o << "#  rho = mass density (NGP)\n";
        o << "#  rhovx = momentum density x-comp (NGP)\n";
        o << "#  rhovy = momentum density y-comp (NGP)\n";
        o << "#  rhovz = momentum density z-comp (NGP)\n";
        o << "#  U1 = total energy density without constant magnetic fields (NGP)\n";
        o << "#  B1x = interaction magnetic field x-comp\n";
        o << "#  B1y = interaction magnetic field y-comp\n";
        o << "#  B1z = interaction magnetic field z-comp\n";
        o << "#  B0x = constant magnetic field x-comp\n";
        o << "#  B0y = constant magnetic field y-comp\n";
        o << "#  B0z = constant magnetic field z-comp\n";
        o << "# PLASMA FILE CONTENTS\n";
        o << "# populations = all (plasma file)\n";
        o << "# m = plasma\n";
        o << "# q = plasma\n";
    }
#ifdef SAVE_POPULATION_AVERAGES
    else if (filetype == 3) {
        o << "# AVERAGE POPULATION FILE CONTENTS\n";
        o << "#  rho = mass density (NGP)\n";
        o << "#  rhovx = momentum density x-comp (NGP)\n";
        o << "#  rhovy = momentum density y-comp (NGP)\n";
        o << "#  rhovz = momentum density z-comp (NGP)\n";
        o << "#  U1 = total energy density without constant magnetic fields (NGP)\n";
        o << "#  B1x = interaction magnetic field x-comp\n";
        o << "#  B1y = interaction magnetic field y-comp\n";
        o << "#  B1z = interaction magnetic field z-comp\n";
        o << "#  B0x = constant magnetic field x-comp\n";
        o << "#  B0y = constant magnetic field y-comp\n";
        o << "#  B0z = constant magnetic field z-comp\n";
        o << "# AVERAGE POPULATION FILE CONTENTS\n";
        o << "# populations = ";
        // Check particle masses and charges
        real mm = Params::pops[popId[0]]->m;
        real qq = Params::pops[popId[0]]->q;
        for (unsigned int i = 0; i < popId.size(); i++) {
            o << Params::pops[popId[i]]->getIdStr() << " ";
            if (Params::pops[popId[i]]->m != mm) {
                o << "(WARNING: particle mass) ";
                errorlog << "WARNING [Tgrid::hcwrite_MHD]: particle mass not same for all populations in a hc-file (" << fn << ")\n";
            }
            if (Params::pops[popId[i]]->q != qq) {
                o << "(WARNING: particle charge) ";
                errorlog << "WARNING [Tgrid::hcwrite_MHD]: particle charge not same for all populations in a hc-file (" << fn << ")\n";
            }
        }
        o << "\n";
        o << "# m = " << Params::pops[popId[0]]->m << "\n";
        o << "# q = " << Params::pops[popId[0]]->q << "\n";
    }
#endif
    else {
        ERRORMSG("internal error");
        return false;
    }
    o << "dim = 3\n";
    o << "ncd = 11\n";
    o << "nsd = 0\n";
    o << "maxnc = " << ncells << "\n";
    o << "realformat = " << (hcFileAsciiFormat ? "ascii" : "float") << "\n";
    o << "n1 = " << nx << "\n";
    o << "n2 = " << ny << "\n";
    o << "n3 = " << nz << "\n";
    o << "xmin0 = " << x_1+bgdx << "\n";  //real grid limits - not ghost cells
    o << "xmin1 = " << y_1+bgdx << "\n";
    o << "xmin2 = " << z_1+bgdx << "\n";
    o << "xmax0 = " << x_1-bgdx + nx*bgdx << "\n";  // ditto
    o << "xmax1 = " << y_1-bgdx + ny*bgdx << "\n";
    o << "xmax2 = " << z_1-bgdx + nz*bgdx << "\n";
    o << "dx = " << bgdx << "\n";
    o << "type = hc\n";
    o << "ncells = " << ncells << "\n";
    o << "nablocks = 0\n";
    o << "bytes_per_int = " << n_int_bytes << "\n";     // no meaning if ascii
    o << "freelist1 = -1234567\n";
    o << "freelist2 = -1234567\n";
    o << "eoh\n";
    o.precision(oldprec);
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) {
            if (hcFileAsciiFormat == false) {
                o.put('N');
                WriteInt(o,-1);
                WriteInt(o,cells[c]->child[0][0][0]->running_index);
            } else {
                o << "N -1 " << cells[c]->child[0][0][0]->running_index << "\n";
            }
        } else {
            cells[c]->writeMHD(o,filetype,popId);
        }
    }
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->writeMHD_children_recursive(o,filetype,popId);
    }
    int i,j,k;
    ForAll(i,j,k) {
        // celltype=0: interior, 1: ghost, 2:dead
        int celltype = (i==0 || i==nx-1) + (j==0 || j==ny-1) + (k==0 || k==nz-1);
        if (celltype > 2) celltype = 2;
        c = flatindex(i,j,k);
        if (hcFileAsciiFormat == false) {
            o.put(EncodeCellInfo(celltype));
            o.put((unsigned char)'\0');
        } else {
            o << celltype << " 1 1 0\n";
        }
    }
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        if (cells[c]->haschildren) {
            const int n = cells[c]->Ncells_recursive() - 1;
            int cnt;
            if (hcFileAsciiFormat == false) {
                const unsigned char ch = EncodeCellInfo(0);
                for (cnt=0; cnt<n; cnt++) {
                    o.put(ch);
                    o.put((unsigned char)'\0');
                }
            } else {
                for (cnt=0; cnt<n; cnt++) o << "0 1 1 0\n";
            }
        }
    }
    if (o.good()) {
        mainlog << "Tgrid::hcwrite_MHD: Wrote \"" << fn << "\"\n";
    } else {
        errorlog << "*** Tgrid::hcwrite_MHD: Could not write \"" << fn << "\" completely - disk full?\n";
    }
    return o.good();
}

//! Write debug hc-file
bool Tgrid::hcwrite_DBUG(const char *fn)
{
    const int ncells = Ncells_with_ghosts();
    ofstream o(fn);
    if (!o.good()) {
        return false;
    }
    const int oldprec = o.precision();
    o.precision(16);
    int c;
    const int Nbase = nx*ny*nz;
    // Enumerate cells
    for (c=0; c<Nbase; c++) {
        cells[c]->running_index = c;
    }
    cell_running_index = Nbase;
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->enum_children_recursive();
    }
    // Magnetic and electric field to cells
    g.FC(Tgrid::FACEDATA_B,Tgrid::CELLDATA_B);
    g.Neumann(Tgrid::CELLDATA_B);
    n_int_bytes = NumberOfBytes(ncells);
    o << "# " << Params::codeVersion << "\n";
    o << "# DBUG FILE CONTENTS\n";
    o << "#  rho = macroparticle occupations\n";
    o << "#  rhovx = CELLDATA_Jx\n";
    o << "#  rhovy = CELLDATA_Jy\n";
    o << "#  rhovz = CELLDATA_Jz\n";
    o << "#  U1 = rho_q\n";
    o << "#  B1x = CELLDATA_Bx\n";
    o << "#  B1y = CELLDATA_By\n";
    o << "#  B1z = CELLDATA_Bz\n";
    if(Params::electronPressure==true) {
        o << "#  B0x = Ex, Electric field Ue x B + grad(pe)\n";
        o << "#  B0y = Ey, Electric field Ue x B + grad(pe)\n";
        o << "#  B0z = Ez, Electric field Ue x B + grad(pe)\n";
    } else {
        o << "#  B0x = Ex, Electric field Ue x B\n";
        o << "#  B0y = Ey, Electric field Ue x B\n";
        o << "#  B0z = Ez, Electric field Ue x B\n";
    }
    o << "# DBUG FILE CONTENTS\n";
    o << "# filename: " << fn << "\n";
    o << "# t = " << Params::t << "\n";
    o << "# gamma = " << Params::gamma << "\n";
    o << "# mu0 = " << Params::mu_0 << "\n";
    o << "# R = " << Params::R_P << "\n";
    o << "dim = 3\n";
    o << "ncd = 11\n";
    o << "nsd = 0\n";
    o << "maxnc = " << ncells << "\n";
    o << "realformat = " << (hcFileAsciiFormat ? "ascii" : "float") << "\n";
    o << "n1 = " << nx << "\n";
    o << "n2 = " << ny << "\n";
    o << "n3 = " << nz << "\n";
    o << "xmin0 = " << x_1+bgdx << "\n";  //real grid limits - not ghost cells
    o << "xmin1 = " << y_1+bgdx << "\n";
    o << "xmin2 = " << z_1+bgdx << "\n";
    o << "xmax0 = " << x_1-bgdx + nx*bgdx << "\n";  // ditto
    o << "xmax1 = " << y_1-bgdx + ny*bgdx << "\n";
    o << "xmax2 = " << z_1-bgdx + nz*bgdx << "\n";
    o << "dx = " << bgdx << "\n";
    o << "type = hc\n";
    o << "ncells = " << ncells << "\n";
    o << "nablocks = 0\n";
    o << "bytes_per_int = " << n_int_bytes << "\n";     // no meaning if ascii
    o << "freelist1 = -1234567\n";
    o << "freelist2 = -1234567\n";
    o << "eoh\n";
    o.precision(oldprec);
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) {
            if (hcFileAsciiFormat == false) {
                o.put('N');
                WriteInt(o,-1);
                WriteInt(o,cells[c]->child[0][0][0]->running_index);
            } else {
                o << "N -1 " << cells[c]->child[0][0][0]->running_index << "\n";
            }
        } else {
            cells[c]->writeDBUG(o);
        }
    }
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->writeDBUG_children_recursive(o);
    }
    int i,j,k;
    ForAll(i,j,k) {
        // celltype=0: interior, 1: ghost, 2:dead
        int celltype = (i==0 || i==nx-1) + (j==0 || j==ny-1) + (k==0 || k==nz-1);
        if (celltype > 2) celltype = 2;
        c = flatindex(i,j,k);
        if (hcFileAsciiFormat == false) {
            o.put(EncodeCellInfo(celltype));
            o.put((unsigned char)'\0');
        } else {
            o << celltype << " 1 1 0\n";
        }
    }
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        if (cells[c]->haschildren) {
            const int n = cells[c]->Ncells_recursive() - 1;
            int cnt;
            if (hcFileAsciiFormat == false) {
                const unsigned char ch = EncodeCellInfo(0);
                for (cnt=0; cnt<n; cnt++) {
                    o.put(ch);
                    o.put((unsigned char)'\0');
                }
            } else {
                for (cnt=0; cnt<n; cnt++) o << "0 1 1 0\n";
            }
        }
    }
    if (o.good()) {
        mainlog << "Tgrid::hcwrite_DBUG: Wrote \"" << fn << "\"\n";
    } else {
        errorlog << "*** Tgrid::hcwrite_DBUG: Could not write \"" << fn << "\" completely - disk full?\n";
    }
    return o.good();
}

//! Write extra hc-file
bool Tgrid::hcwrite_EXTRA(string fileName, ScalarField* s, string hcHeaderStr)
{
    const char* fn = fileName.c_str();
    const int ncells = Ncells_with_ghosts();
    ofstream o(fn);
    if (!o.good()) {
        return false;
    }
    const int oldprec = o.precision();
    o.precision(16);
    int c;
    const int Nbase = nx*ny*nz;
    // Enumerate cells
    for (c=0; c<Nbase; c++) {
        cells[c]->running_index = c;
    }
    cell_running_index = Nbase;
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->enum_children_recursive();
    }
    n_int_bytes = NumberOfBytes(ncells);
    o << "# " << Params::codeVersion << "\n";
    o << "# EXTRA FILE CONTENTS\n";
    o << "#  rho = scalar field: " << hcHeaderStr << "\n";
    o << "# EXTRA FILE CONTENTS\n";
    o << "# filename: " << fileName << "\n";
    o << "# gamma = " << Params::gamma << "\n";
    o << "# mu0 = " << Params::mu_0 << "\n";
    o << "# R = " << Params::R_P << "\n";
    o << "dim = 3\n";
    o << "ncd = 1\n";
    o << "nsd = 0\n";
    o << "maxnc = " << ncells << "\n";
    o << "realformat = " << (hcFileAsciiFormat ? "ascii" : "float") << "\n";
    o << "n1 = " << nx << "\n";
    o << "n2 = " << ny << "\n";
    o << "n3 = " << nz << "\n";
    o << "xmin0 = " << x_1+bgdx << "\n";
    o << "xmin1 = " << y_1+bgdx << "\n";
    o << "xmin2 = " << z_1+bgdx << "\n";
    o << "xmax0 = " << x_1-bgdx + nx*bgdx << "\n";
    o << "xmax1 = " << y_1-bgdx + ny*bgdx << "\n";
    o << "xmax2 = " << z_1-bgdx + nz*bgdx << "\n";
    o << "dx = " << bgdx << "\n";
    o << "type = hc\n";
    o << "ncells = " << ncells << "\n";
    o << "nablocks = 0\n";
    o << "bytes_per_int = " << n_int_bytes << "\n";     // no meaning if ascii
    o << "freelist1 = -1234567\n";
    o << "freelist2 = -1234567\n";
    o << "eoh\n";
    o.precision(oldprec);
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) {
            if (hcFileAsciiFormat == false) {
                o.put('N');
                WriteInt(o,-1);
                WriteInt(o,cells[c]->child[0][0][0]->running_index);
            } else {
                o << "N -1 " << cells[c]->child[0][0][0]->running_index << "\n";
            }
        } else {
            cells[c]->writeEXTRA(o,s);
        }
    }
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->writeEXTRA_children_recursive(o,s);
    }
    int i,j,k;
    ForAll(i,j,k) {
        // celltype=0: interior, 1: ghost, 2:dead
        int celltype = (i==0 || i==nx-1) + (j==0 || j==ny-1) + (k==0 || k==nz-1);
        if (celltype > 2) celltype = 2;
        c = flatindex(i,j,k);
        if (hcFileAsciiFormat == false) {
            o.put(EncodeCellInfo(celltype));
            o.put((unsigned char)'\0');
        } else {
            o << celltype << " 1 1 0\n";
        }
    }
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        if (cells[c]->haschildren) {
            const int n = cells[c]->Ncells_recursive() - 1;
            int cnt;
            if (hcFileAsciiFormat == false) {
                const unsigned char ch = EncodeCellInfo(0);
                for (cnt=0; cnt<n; cnt++) {
                    o.put(ch);
                    o.put((unsigned char)'\0');
                }
            } else {
                for (cnt=0; cnt<n; cnt++) o << "0 1 1 0\n";
            }
        }
    }
    if (o.good()) {
        mainlog << "Tgrid::hcwrite_EXTRA: Wrote \"" << fn << "\"\n";
    } else {
        errorlog << "*** Tgrid::hcwrite_EXTRA: Could not write \"" << fn << "\" completely - disk full?\n";
    }
    return o.good();
}

#ifdef SAVE_PARTICLE_CELL_SPECTRA

//! Write spectra hc-file
bool Tgrid::hcwrite_SPECTRA(const char *fn,string ascbin,vector<int> popId)
{
    if(ascbin.compare("binary") == 0) {
        hcFileAsciiFormat = false;
    } else if(ascbin.compare("ascii") == 0) {
        hcFileAsciiFormat = true;
    }
    const int ncells = Ncells_with_ghosts();
    ofstream o(fn);
    if (!o.good()) {
        return false;
    }
    const int oldprec = o.precision();
    o.precision(16);
    int c;
    const int Nbase = nx*ny*nz;
    // Enumerate cells
    for (c=0; c<Nbase; c++) {
        cells[c]->running_index = c;
    }
    cell_running_index = Nbase;
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->enum_children_recursive();
    }
    n_int_bytes = NumberOfBytes(ncells);
    o << "# " << Params::codeVersion << "\n";
    o << "# filename: " << fn << "\n";
    o << "# t = " << Params::t << "\n";
    o << "# gamma = " << Params::gamma << "\n";
    o << "# mu0 = " << Params::mu_0 << "\n";
    o << "# R = " << Params::R_P << "\n";
    o << "# SPECTRA FILE CONTENTS\n";
    for(int i=0; i<Params::spectraNbins; ++i) {

        o << "# " << int2string(i,4) << ". bin: E (eV) = " << Params::spectraEnergyBins_eV[i] << " " << Params::spectraEnergyBins_eV[i+1] << "\n";
    }
    o << "# SPECTRA FILE CONTENTS\n";
    o << "# populations = ";
    // Check particle masses and charges
    real mm = Params::pops[popId[0]]->m;
    real qq = Params::pops[popId[0]]->q;
    for (unsigned int i = 0; i < popId.size(); i++) {
        o << Params::pops[popId[i]]->getIdStr() << " ";
        if (Params::pops[popId[i]]->m != mm) {
            o << "(WARNING: particle mass) ";
            errorlog << "WARNING [Tgrid::hcwrite_SPECTRA]: particle mass not same for all populations in a hc-file (" << fn << ")\n";
        }
        if (Params::pops[popId[i]]->q != qq) {
            o << "(WARNING: particle charge) ";
            errorlog << "WARNING [Tgrid::hcwrite_SPECTRA]: particle charge not same for all populations in a hc-file (" << fn << ")\n";
        }
    }
    o << "\n";
    o << "# m = " << Params::pops[popId[0]]->m << "\n";
    o << "# q = " << Params::pops[popId[0]]->q << "\n";
    o << "# spectra = 1\n";
    o << "# Emin_eV = " << Params::spectraEmin_eV << "\n";
    o << "# Emax_eV = " << Params::spectraEmax_eV << "\n";
    o << "# Nbins = " << Params::spectraNbins << "\n";
    o << "# logBins = " << Params::spectraLogBins << "\n";
    o << "# EminAllLowEnergies = " << Params::spectraEminAll << "\n";
    o << "# EmaxAllHighEnergies = " << Params::spectraEmaxAll << "\n";
    o << "# spectraMethod = " << Params::spectraMethod << "\n";
    o << "# spectraUnit = " << Params::spectraUnit << "\n";
    o << "dim = 3\n";
    o << "ncd = " << Params::spectraNbins << "\n";
    o << "nsd = 0\n";
    o << "maxnc = " << ncells << "\n";
    o << "realformat = " << (hcFileAsciiFormat ? "ascii" : "float") << "\n";
    o << "n1 = " << nx << "\n";
    o << "n2 = " << ny << "\n";
    o << "n3 = " << nz << "\n";
    o << "xmin0 = " << x_1+bgdx << "\n";  //real grid limits - not ghost cells
    o << "xmin1 = " << y_1+bgdx << "\n";
    o << "xmin2 = " << z_1+bgdx << "\n";
    o << "xmax0 = " << x_1-bgdx + nx*bgdx << "\n";  // ditto
    o << "xmax1 = " << y_1-bgdx + ny*bgdx << "\n";
    o << "xmax2 = " << z_1-bgdx + nz*bgdx << "\n";
    o << "dx = " << bgdx << "\n";
    o << "type = hc\n";
    o << "ncells = " << ncells << "\n";
    o << "nablocks = 0\n";
    o << "bytes_per_int = " << n_int_bytes << "\n";     // no meaning if ascii
    o << "freelist1 = -1234567\n";
    o << "freelist2 = -1234567\n";
    o << "eoh\n";
    o.precision(oldprec);
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) {
            if (hcFileAsciiFormat == false) {
                o.put('N');
                WriteInt(o,-1);
                WriteInt(o,cells[c]->child[0][0][0]->running_index);
            } else {
                o << "N -1 " << cells[c]->child[0][0][0]->running_index << "\n";
            }
        } else {
            cells[c]->writeSPECTRA(o,popId);
        }
    }
    for (c=0; c<Nbase; c++) {
        if (cells[c]->haschildren) cells[c]->writeSPECTRA_children_recursive(o,popId);
    }
    int i,j,k;
    ForAll(i,j,k) {
        // celltype=0: interior, 1: ghost, 2:dead
        int celltype = (i==0 || i==nx-1) + (j==0 || j==ny-1) + (k==0 || k==nz-1);
        if (celltype > 2) celltype = 2;
        c = flatindex(i,j,k);
        if (hcFileAsciiFormat == false) {
            o.put(EncodeCellInfo(celltype));
            o.put((unsigned char)'\0');
        } else {
            o << celltype << " 1 1 0\n";
        }
    }
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        if (cells[c]->haschildren) {
            const int n = cells[c]->Ncells_recursive() - 1;
            int cnt;
            if (hcFileAsciiFormat == false) {
                const unsigned char ch = EncodeCellInfo(0);
                for (cnt=0; cnt<n; cnt++) {
                    o.put(ch);
                    o.put((unsigned char)'\0');
                }
            } else {
                for (cnt=0; cnt<n; cnt++) o << "0 1 1 0\n";
            }
        }
    }
    if (o.good()) {
        mainlog << "Tgrid::hcwrite_SPECTRA: Wrote \"" << fn << "\"\n";
    } else {
        errorlog << "*** Tgrid::hcwrite_SPECTRA: Could not write \"" << fn << "\" completely - disk full?\n";
    }
    return o.good();
}

#endif

// =================================================================================
// ================================= BREAKPOINTING =================================
// =================================================================================


//! (BREAKPOINTING) Write data
template <class T>
ostream& writeData(ostream& os, T& data)
{
    os.write(reinterpret_cast<const char*>(&data), sizeof(data));
    return os;
}

//! (BREAKPOINTING) Read data
template <class T>
istream& readData(istream& is, T& data)
{
    is.read(reinterpret_cast<char*>(&data), sizeof(data));
    return is;
}

//! (BREAKPOINTING) Write particle
struct Tgrid::writeParticle {

    writeParticle(ostream& o) : os_(o) { }
    bool operator() (TLinkedParticle& p) const {
        writeData(os_, p.x);
        writeData(os_, p.y);
        writeData(os_, p.z);
        writeData(os_, p.vx);
        writeData(os_, p.vy);
        writeData(os_, p.vz);
        writeData(os_, p.w);
        writeData(os_, p.popid);
        return true;
    }
private:
    ostream& os_;
};

//! (BREAKPOINTING) Write magnetic field
struct Tgrid::writeMagneticField {
    writeMagneticField(ostream& o) : os_(o) { }
    void operator() (Tcell& cell) const {
        const int dims = 3;
        for (int dim = 0; dim < dims; ++dim)
            if (!cell.isrefined_face(dim, 1)) {
                TFacePtr face = cell.face[dim][1];
                if (face != 0)
                    writeData(os_, face->facedata[FACEDATA_B]);
            } else {
                const int faces = 4;
                for (int f = 0; f < faces; ++f) {
                    Trefintf* ref = cell.refintf[dim][1];
                    if (ref != 0)
                        writeData(os_, ref->face[f]->facedata[FACEDATA_B]);
                }
            }
    }
private:
    ostream& os_;
};

//! (BREAKPOINTING) Read magnetic field
struct Tgrid::readMagneticField {
    readMagneticField(istream& i) : is_(i) { }
    void operator() (Tcell& cell) {
        const int dims = 3;
        for (int dim = 0; dim < dims; ++dim)
            if (!cell.isrefined_face(dim, 1)) {
                TFacePtr face = cell.face[dim][1];
                if (face != 0)
                    readData(is_, face->facedata[FACEDATA_B]);
            } else {
                const int faces = 4;
                for (int f = 0; f < faces; ++f) {
                    Trefintf* ref = cell.refintf[dim][1];
                    if (ref != 0)
                        readData(is_, ref->face[f]->facedata[FACEDATA_B]);
                }
            }
    }
private:
    istream& is_;
};

//! (BREAKPOINTING) Write breakpoint
void Tgrid::dumpState(ostream& os)
{
    portrand.save(os); // write random number generator state
    writeData(os, Params::t);
    writeData(os, Params::cnt_dt);
    writeData(os, n_particles);
    // write particles
    particle_pass(writeParticle(os));
    // write magnetic field
    cellPass(writeMagneticField(os));
#ifndef NO_DIAGNOSTICS
    for (unsigned int i = 0; i < Params::pops.size(); ++i) {
        writeData(os, Params::diag.pCounter[i]);
    }
    writeData(os, fieldCounter);
#endif
}

//! (BREAKPOINTING) Read breakpoint
void Tgrid::readState(istream& is)
{
    portrand.load(is); // restore random number generator state
    readData(is, Params::t);
    readData(is, Params::cnt_dt);
    int nPart;
    readData(is, nPart);
    // read particles
    shortreal x, y, z, vx, vy, vz, w;
    int popid;
    for (int n = 0; n < nPart; ++n) {
        readData(is, x);
        readData(is, y);
        readData(is, z);
        readData(is, vx);
        readData(is, vy);
        readData(is, vz);
        readData(is, w);
        readData(is, popid);
        if(popid < 0 || popid > static_cast<int>(Params::pops.size())) {
            WARNINGMSG2("readState: particle with a bad popid, ignoring",popid);
            continue;
        }
        addparticle(x, y, z, vx, vy, vz, w, popid);
    }
    // read magnetic field
    readMagneticField mReader(is);
    cellPass(mReader);
#ifndef NO_DIAGNOSTICS
    // This may not work as planned
    /*for (unsigned int i = 0; i < Params::pops.size(); ++i) {
        readData(is, Params::diag.pCounter[i]);
    }
    readData(is, fieldCounter);*/
#endif
}

// =================================================================================
// =================================== PARTICLES ===================================
// =================================================================================

//! Add particle into in the grid with given parameters
void Tgrid::addparticle(shortreal x, shortreal y, shortreal z,
                        shortreal vx, shortreal vy, shortreal vz,
                        shortreal w, int popid, bool inject)
{
    const shortreal r[3] = {x,y,z};
    TCellPtr c = findcell(r);
    if (!c) {
        errorlog << "WARNING: Tgrid::addparticle" << Tr3v(r).toString()
                 << " idStr=" << Params::pops[popid]->getIdStr()
                 << " out of box (not created)\n";
        return;
    }
    c->plist.add(x,y,z,vx,vy,vz,w,popid);
#ifndef NO_DIAGNOSTICS
    if(inject==true) {
        Params::diag.pCounter[popid]->increaseInjectCounters(vx,vy,vz,w);
    }
#endif
    n_particles++;
}

//! Find the particle's list
TParticleList *Tgrid::find_plist(const TLinkedParticle& p)
{
    const shortreal r[3] = {p.x,p.y,p.z};
    TCellPtr c = findcell(r);
    if (!c) {
        errorlog << "WARNING: Tgrid::find_plist" << Tr3v(r).toString() << " out of box\n";
        return 0;
    }
    return &(c->plist);
}

//! Returns the number of particles deleted
int Tgrid::Tcell::particle_pass_recursive(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), bool relocate)
{
    int ndel = 0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) ndel+= child[0][0][ch]->particle_pass_recursive(op,relocate);
    } else {
        // Set arguments for the particle pass function
        ParticlePassArgs a;
        a.rho_q = this->rho_q;
        a.size = this->size;
        if (relocate) {
            ndel = plist.pass_with_relocate(op,a);
        } else {
            ndel = plist.pass(op,a);
        }
    }
    return ndel;
}

//! Pass all particles in the list to the function op
int Tgrid::particle_pass(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), bool relocate)
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

//! Return the number of macro particles (recursive)
int Tgrid::Tcell::Nparticles_recursive() const
{
    int result = 0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) result+= child[0][0][ch]->Nparticles_recursive();
    } else {
        result = plist.Nparticles();
    }
    return result;
}

//! Return the number of macro particles
int Tgrid::Nparticles() const
{
    return n_particles;
}

/** \brief Split&Join probability function
 *
 * If Params::splitJoinDeviation[1]==1 use new stepfunction probability.
 * If the number of particles inside this cell deviates too much from the 'optimal' value n_target,
 * try splitting or joining. The number of particles actually split and joined is returned
 * (counting the cells' children also). Probabilistic algorithm. Ensures that cell never has fewer
 * than 1 particle (splitting probability is 1 if there is only 1 particle) (except in the
 * extremely unlikely case that this splitting fails because the particle is too close (5%)
 * to the cell boundary).
 */
void Tgrid::Tcell::split_and_join_recursive(int& nsplit, int& njoined)
{
    nsplit = njoined = 0;
    if (forbid_psplit) return;
    if (haschildren) {
        int ch;
        int nsplit1,njoined1;
        for (ch=0; ch<8; ch++) {
            child[0][0][ch]->split_and_join_recursive(nsplit1,njoined1);
            nsplit+= nsplit1;
            njoined+= njoined1;
        }
    } else {
        const int npart = plist.Nparticles(), n_target = Params::macroParticlesPerCell;
        if (npart==n_target || npart==0) {
            return;
        }
        gridreal mincell[3],maxcell[3];
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        const real halfdx=0.5*size;
        for (int d=0; d<3; d++) {
            mincell[d] = centroid[d] - halfdx;
            maxcell[d] = centroid[d] + halfdx;
        }
#else
        for (int d=0; d<3; d++) {
            mincell[d] = -0.5*abs(sph_dl[d]);
            maxcell[d] =  0.5*abs(sph_dl[d]);
        }
#endif
        if (Params::splitJoinDeviation[1]==1) { // new method
            real n_split=(1-Params::splitJoinDeviation[0])*n_target;
            real n_join=(1+Params::splitJoinDeviation[0])*n_target;
            if (npart<n_split) {
                int splits=1;
                if (npart<n_split/3.0) {
                    splits=5;
                } else {
                    if (npart<n_split/2.0) {
                        splits=3;
                    }
                }
                // do splits
                if (Params::useMacroParticleSplitting) {
                    for (int i_split=0; i_split<splits; i_split++) {
                        nsplit += Params::splittingFunction.doSplitting(mincell,maxcell,plist);
                    }
                }
            } else {
                if (npart>n_join) {
                    int joins=1;
                    //do fast joins, regular joins scale as NxN, stalls if large number of particles are created in a cell,
                    // e.g., in base cell in an initialization with many refinements.
                    if (npart>3.0*n_join) {
                        joins = npart - int(3*n_join); //do them all at once -save time in progation
                    } else {
                        if (npart>2.0*n_join) {
                            joins=3;
                        }
                    }
                    // do joins
                    if (Params::useMacroParticleJoining) {
                        if (joins>8) { //do fast joins
                            for (int i_join=0; i_join<joins; i_join++) {
                                njoined -= Params::joiningFunction.doJoining(mincell,maxcell,plist,1);
                            }
                        } else {
                            for (int i_join=0; i_join<joins; i_join++) {
                                njoined -= Params::joiningFunction.doJoining(mincell,maxcell,plist,0);
                            }
                        }
                    }
                } else return; // n_split<=npart<+=n_join
            }
        } else { // original method
            real sjprob;
            if (npart<n_target) {
                sjprob = 1.0 - (Params::splitjoin_a*(npart-1.0))/(n_target-1.0);
            } else {
                sjprob = (npart-n_target)/real(n_target);
                sjprob -= (Params::splitjoin_a-1)*(n_target-1)/(n_target+Params::splitjoin_a-1.0);
            }
            if (sjprob <= 0) {
                return;
            }
            if (sjprob < 1) {
                if (uniformrnd() > sjprob) {
                    return;
                }
            }
            int n=0;
            if (npart < n_target && Params::useMacroParticleSplitting) {
                n = Params::splittingFunction.doSplitting(mincell,maxcell,plist);
            } else if (Params::useMacroParticleJoining) {
                n = Params::joiningFunction.doJoining(mincell,maxcell,plist,0);
            }
            if (n > 0) {
                nsplit+= n;
            } else if (n < 0) {
                njoined-= n;
            }
        }
    }
}

//! Macro particles split&join 
void Tgrid::split_and_join(int& nsplit, int& njoined)
{
    int i,j,k,nsplit1,njoined1;
    nsplit = njoined = 0;
    ForAll(i,j,k) {
        cells[flatindex(i,j,k)]->split_and_join_recursive(nsplit1,njoined1);
        nsplit+= nsplit1;
        njoined+= njoined1;
    }
    n_particles -= njoined;
}

//! Set cells where split&join is forbidden (recursive)
int Tgrid::Tcell::forbid_split_and_join_recursive(ForbidSplitAndJoinProfile forb)
{
    int result = 0;
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) result+= child[0][0][ch]->forbid_split_and_join_recursive(forb);
    } else {
        int dirx,diry,dirz,d;
        gridreal r[3];
        const gridreal halfdx = 0.5*size;
        bool do_forbid = false;
        for (dirx=0; dirx<2; dirx++) for (diry=0; diry<2; diry++) for (dirz=0; dirz<2; dirz++) {
                    for (d=0; d<3; d++) r[d] = centroid[d];
                    r[0]+= (dirx==0) ? -halfdx : +halfdx;
                    r[1]+= (diry==0) ? -halfdx : +halfdx;
                    r[2]+= (dirz==0) ? -halfdx : +halfdx;
                    if (forb.getValue(r)) {
                        do_forbid = true;
                        goto exitloops;
                    }
                }
exitloops:
        forbid_psplit = do_forbid;
        if (do_forbid) result++;
    }
    return result;
}

//! Set cells where split&join is forbidden
int Tgrid::forbid_split_and_join(ForbidSplitAndJoinProfile forb)
{
    int i,j,k;
    int result = 0;
    ForAll(i,j,k) {
        result+= cells[flatindex(i,j,k)]->forbid_split_and_join_recursive(forb);
    }
    return result;
}

//! Start temporal average (recursive)
void Tgrid::Tcell::begin_average_recursive()
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->begin_average_recursive();
    } else {
        ave_nc = 0.0;
        int d;
        for (d=0; d<3; d++) {
            if (neighbour[d][1] == 0) continue;
            if (isrefined_face(d,1)) {
                int f;
                for (f=0; f<4; f++) refintf[d][1]->face[f]->facedata[FACEDATA_AVEB] = 0;
            } else {
                face[d][1]->facedata[FACEDATA_AVEB] = 0;
            }
        }
#ifdef SAVE_POPULATION_AVERAGES
        for(int i = 0; i < Params::POPULATIONS; i++) {
            pop_ave_n[i] = 0.0;
            pop_ave_vx[i] = 0.0;
            pop_ave_vy[i] = 0.0;
            pop_ave_vz[i] = 0.0;
        }
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
        for(unsigned int i=0; i<spectra.size(); ++i) {
            for(unsigned int j=0; j<spectra[i].size(); ++j) {
                spectra[i][j] = 0.0;
            }
        }
#endif
    }
}

//! End temporal average (recursive)
void Tgrid::Tcell::end_average_recursive(real inv_ave_ntimes)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->end_average_recursive(inv_ave_ntimes);
    } else {
        ave_nc*= inv_ave_ntimes;
        int d;
        for (d=0; d<3; d++) {
            if (neighbour[d][1] == 0) continue;
            if (isrefined_face(d,1)) {
                int f;
                for (f=0; f<4; f++) refintf[d][1]->face[f]->facedata[FACEDATA_AVEB]*= inv_ave_ntimes;
            } else {
                face[d][1]->facedata[FACEDATA_AVEB]*= inv_ave_ntimes;
            }
        }
#if defined(SAVE_POPULATION_AVERAGES) || defined(SAVE_PARTICLE_CELL_SPECTRA)
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
        gridreal invvol = 1.0/(size*size*size);
#else
        gridreal invvol = 1.0/abs(sph_dV);
#endif
#endif
#ifdef SAVE_POPULATION_AVERAGES
        for (int i = 0; i < Params::POPULATIONS; ++i) {
            gridreal inv_sum_w = 0;
            if(pop_ave_n[i] > 0) {
                inv_sum_w = 1.0/pop_ave_n[i];
            }
            pop_ave_n[i] *= invvol;
            pop_ave_vx[i] *= inv_sum_w;
            pop_ave_vy[i] *= inv_sum_w;
            pop_ave_vz[i] *= inv_sum_w;
            pop_ave_n[i] *= inv_ave_ntimes;
        }
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
        for(unsigned int i=0; i<spectra.size(); ++i) {
            const real norm = invvol*inv_ave_ntimes/(Params::spectra_dE_eV[i]*4*pi);
            for(unsigned int j=0; j<spectra[i].size(); ++j) {
                spectra[i][j] *= norm;
            }
        }
#endif
    }
}

//! Start temporal average
void Tgrid::begin_average()
{
    int i,j,k;
    ForAll(i,j,k) cells[flatindex(i,j,k)]->begin_average_recursive();
    ave_ntimes = 0;
}

//! End temporal average
bool Tgrid::end_average()
{
    if (ave_ntimes <= 0) {
        errorlog << "*** Tgrid::end_average: ave_ntimes=" << ave_ntimes << "\n";
        return false;
    }
    const real inv_ave_ntimes = 1.0/ave_ntimes;
    int i,j,k;
    ForAll(i,j,k) cells[flatindex(i,j,k)]->end_average_recursive(inv_ave_ntimes);
    ave_ntimes = 0;
    return true;
}

//! Prepare Probability Density Functions in the grid (recursive)
void Tgrid::Tcell::prepare_PDF_recursive(ScalarField* pdffunc, Tgrid* g)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->prepare_PDF_recursive(pdffunc,g);
    } else {
        const real fval = pdffunc->getValue(centroid);
        g->pdftables[g->n_pdftables].cpdf[g->pdftable_iterator] = fval*(size*size*size);
        g->pdftables[0].cellptrs[g->pdftable_iterator++] = this;
    }
}

//! Prepare Probability Density Functions in the grid
void Tgrid::prepare_PDF(ScalarField* pdffunc, TPDF_ID& pdfid, real& cumsumvalue)
{
    if (n_pdftables >= MAX_PDFTABLES) {
        errorlog << "ERROR [Tgrid::prepare_PDF]: too many PDF tables already allocated (max is " << MAX_PDFTABLES << ")\n";
        doabort();
    }
    pdfid = n_pdftables;
    if (n_pdftables == 0) {
        const int n = Ncells_without_ghosts();
        pdftables[n_pdftables].n = n;
        pdftables[n_pdftables].cellptrs = new TCellPtr [n];
    } else {
        pdftables[n_pdftables].n = pdftables[0].n;
        pdftables[n_pdftables].cellptrs = 0;
    }
    pdftables[n_pdftables].cpdf = new shortreal [pdftables[n_pdftables].n];
    pdftable_iterator = 0;
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->prepare_PDF_recursive(pdffunc,this);
    }
    for (i=0; i<=n_pdftables; i++) pdftables[i].n = pdftable_iterator + 1;
    // plus 1 because cumsum table contains one more element (because first element is 0)
    if (pdftable_iterator+1 > pdftables[0].n) {
        errorlog << "ERROR [Tgrid::prepare_PDF]: internal error in prepare_PDF\n";
        doabort();
    }
    // now cpdf[0..pdftable_iterator-1] contains pdf, cumsum it
    cumsumvalue = 0;
    for (c=0; c<pdftable_iterator; c++) {
        cumsumvalue+= pdftables[n_pdftables].cpdf[c];
        pdftables[n_pdftables].cpdf[c] = cumsumvalue;
    }
    // shift one step to right, and set first to 0:
    for (c=pdftable_iterator; c>=1; c--) pdftables[n_pdftables].cpdf[c] = pdftables[n_pdftables].cpdf[c-1];
    pdftables[n_pdftables].cpdf[0] = 0;
    if (cumsumvalue <= 0 || !finite(cumsumvalue)) {
        errorlog << "ERROR [Tgrid::prepare_PDF]: bad cumsumvalue (" << cumsumvalue << ")\n";
        doabort();
    } else {
        const real coeff = 1/cumsumvalue;
        for (c=0; c<pdftable_iterator; c++) pdftables[n_pdftables].cpdf[c]*= coeff;
        // now first entry is 0 and last 1
    }
    n_pdftables++;
}

/** \brief Locate
 *
 * Given array xx[1..n] and given value x, returns value j such that x is between
 * xx[j] and xx[j+1]. xx must be monotonic, either increasing or decreasing.
 * j=0 is returned to indicate that x is out of range.
 */
static int locate(const shortreal xx[], int n, real x)
{
    int ju=n+1,jm,jl=0,result=0;
    const bool sscnd = (xx[n] >= xx[1]);
    while (ju-jl > 1) {
        jm = (ju+jl) >> 1; // compute a midpoint
        if ((x >= xx[jm]) == sscnd)
            jl = jm; // replace either the lower limit
        else
            ju = jm; // or the upper limit, as appropriate
    }
    if (x == xx[1])
        result = 1;
    else if (x == xx[n])
        result = n-1;
    else
        result = jl;
    if (result >= n) result = 0; // always return 0 for out-of-range condition
    return result;
}

//! Generate random point in a cell
void Tgrid::Tcell::generate_random_point(gridreal r[3])
{
    int d;
    for (d=0; d<3; d++) {
        r[d]= centroid[d] + (1.0-2*Params::box_eps)*size*(uniformrnd()-0.5);
    }
}

//! Generate random point according to given PDF
void Tgrid::generate_random_point(const TPDF_ID& pdfid, gridreal r[3])
{
    const real x = uniformrnd();
    if (pdfid < 0 || pdfid >= n_pdftables) {
        errorlog << "*** Tgrid::generate_random_point(pdfid=" << pdfid << ") is out of range 0.." << n_pdftables-1 << "\n";
        return;
    }
    const int j = locate(&pdftables[pdfid].cpdf[-1],pdftables[pdfid].n,x);
    if (j == 0) {
        errorlog << "*** Internal error in generate_random_point(), x=" << x << "\n";
        return;
    }
    TCellPtr cptr = pdftables[0].cellptrs[j-1];
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    cptr->generate_random_point(r);
#else
    cptr->generate_random_point(r);
#endif
}

//! NGP interpolation
void Tgrid::cellintpol(const shortreal r[3], TCellDataSelect s, real result[3])
{
    int d;
    Tcell *const c = findcell(r);
    if(!c) {
        ERRORMSG("NULL cell pointer");
        doabort();
    }
    for (d=0; d<3; d++) result[d] = c->celldata[s][d];
    saved_cellptr = c;
}

//! NGP interpolation (use saved cellptr)
void Tgrid::cellintpol(TCellDataSelect s, real result[3]) const
{
    int d;
    for (d=0; d<3; d++) result[d] = saved_cellptr->celldata[s][d];
}

//! Interpolation of fluid parameters in a cell (use saved cellptr)
void Tgrid::cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, vector<int> popId)
{
    saved_cellptr->cellintpol_fluid(n, vx, vy, vz, P, popId);
}

//! Interpolation of fluid parameters in a cell
void Tgrid::cellintpol_fluid(const shortreal r[3], real& n, real& vx, real& vy, real& vz, real& P, vector<int> popId)
{
    Tcell *const c = findcell(r);
    if(!c) {
        ERRORMSG("NULL cell pointer");
        doabort();
    }
    saved_cellptr = c;
    cellintpol_fluid(n, vx, vy, vz, P, popId);
}

//! Constructor
FieldCounter::FieldCounter()
{
    reset();
}

//! Reset field counters
void FieldCounter::reset()
{
    cutRateE = 0.0;
    cutRateRhoQ = 0.0;
    cutRateUe = 0.0;
    cutRateMaxVw = 0.0;
    resetTimestep = Params::cnt_dt;
}

//! Finalize results
void FieldCounter::finalize()
{
    real timeSteps = static_cast<real>(1+Params::cnt_dt - Tgrid::fieldCounter.resetTimestep);
    if(timeSteps > 0) {
        cutRateE /= timeSteps;
        cutRateRhoQ /= timeSteps;
        cutRateUe /= timeSteps;
        cutRateMaxVw /= timeSteps;
    } else {
        cutRateE = 0;
        cutRateRhoQ = 0;
        cutRateUe = 0;
        cutRateMaxVw = 0;
    }
}

// =================================================================================
// ============================ JSTAG AND NODEUE SCHEMES ===========================
// =================================================================================

/** \brief Face boundary conditions
 *
 * Copies J to ghost faces, so J can be calculated for nodes in contact with
 * ghosts in the same way as for inner nodes.
 */
void Tgrid::boundary_faces(TFaceDataSelect cs)
{
    int i,j,k;
    i = 0;
    for(j=0; j<ny; j++) for(k=0; k<nz; k++) { // -X-wall
            if(cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i+1,j,k)]->face[1][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i+1,j,k)]->face[2][1]->facedata[cs];
        }
    i = nx-1;
    for(j=0; j<ny; j++) for(k=0; k<nz; k++) { // +X-wall
            if(cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i-1,j,k)]->face[1][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i-1,j,k)]->face[2][1]->facedata[cs];
        }
    j = 0;
    for(i=0; i<nx; i++) for(k=0; k<nz; k++) { // -Y-wall
#ifdef PERIODIC_FIELDS_Y
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,ny-2,k)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,ny-2,k)]->face[2][1]->facedata[cs];
#else
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j+1,k)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,j+1,k)]->face[2][1]->facedata[cs];
#endif
        }
    j = ny-1;
    for(i=0; i<nx; i++) for(k=0; k<nz; k++) { // +Y-wall
#ifdef PERIODIC_FIELDS_Y
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,1,k)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,1,k)]->face[2][1]->facedata[cs];
#else
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j-1,k)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,j-1,k)]->face[2][1]->facedata[cs];
#endif
        }
    k = 0;
    for(i=0; i<nx; i++) for(j=0; j<ny; j++) { // -Z-wall
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j,k+1)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i,j,k+1)]->face[1][1]->facedata[cs];
        }
    k = nz-1;
    for(i=0; i<nx; i++) for(j=0; j<ny; j++) { // +Z-wall
            if(cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j,k-1)]->face[0][1]->facedata[cs];
            if(cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i,j,k-1)]->face[1][1]->facedata[cs];
        }
}

/** \brief J component in direction edir at a node
 * 
 * Calculation steps:
 *
 * 1. calculate J along a cell edge using line integral and faces that are not perpendicular to edge.
 * 2.interpolate x (y,z) values from the two neighbouring x (y,z) edges of node.
 *
 * NOTE: Not completely correct for nodes that are on the boundaries of the
 * refinements, since the different face sizes are not taken into account.
 */
real Tgrid::Tnode::jcomponent(int edir, real dx)
{
    real bf1, bf2, bf3, bf4; // B on faces for line integral.
    real invdxmu0 = 1./(dx*Params::mu_0);
    real epos=0.;
    real eneg=0.;
    real intp;
    bool posf=false;
    bool negf=false;
    switch(edir) {
    case 0 :  // x-edges
        if(cell[1][0][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][0][1]
           or cell[1][0][1]!=cell[1][1][1]) { // positive side edge exists (it might not because of refinement)?
            posf=true;
            bf1 = (cell[1][0][0]!=cell[1][0][1]) ? cell[1][0][0]->faceave(2,1,FACEDATA_B)
                  : cell[1][0][0]->celldata[CELLDATA_B][2];
            bf2 = (cell[1][0][0]!=cell[1][1][0]) ? cell[1][0][0]->faceave(1,1,FACEDATA_B)
                  : cell[1][0][0]->celldata[CELLDATA_B][1];
            bf3 = (cell[1][1][0]!=cell[1][1][1]) ? cell[1][1][0]->faceave(2,1,FACEDATA_B)
                  : cell[1][1][0]->celldata[CELLDATA_B][2];
            bf4 = (cell[1][1][1]!=cell[1][0][1]) ? cell[1][1][1]->faceave(1,0,FACEDATA_B)
                  : cell[1][1][1]->celldata[CELLDATA_B][1];
            epos = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        if(cell[0][0][0]!=cell[0][1][0] or cell[0][1][0]!=cell[0][0][1]
           or cell[0][0][1]!=cell[0][1][1]) { // negative side edge exists?
            negf=true;
            bf1 = (cell[0][0][0]!=cell[0][0][1]) ? cell[0][0][0]->faceave(2,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][2];
            bf2 = (cell[0][0][0]!=cell[0][1][0]) ? cell[0][0][0]->faceave(1,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][1];
            bf3 = (cell[0][1][0]!=cell[0][1][1]) ? cell[0][1][0]->faceave(2,1,FACEDATA_B)
                  : cell[0][1][0]->celldata[CELLDATA_B][2];
            bf4 = (cell[0][1][1]!=cell[0][0][1]) ? cell[0][1][1]->faceave(1,0,FACEDATA_B)
                  : cell[0][1][1]->celldata[CELLDATA_B][1];
            eneg = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        break;
    case 1 :  // y-edges
        if(cell[0][1][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][1][1]
           or cell[1][1][1]!=cell[0][1][1]) { // positive side edge exists?
            posf=true;
            bf1 = (cell[0][1][0]!=cell[1][1][0]) ? cell[0][1][0]->faceave(0,1,FACEDATA_B)
                  : cell[0][1][0]->celldata[CELLDATA_B][0];
            bf2 = (cell[0][1][0]!=cell[0][1][1]) ? cell[0][1][0]->faceave(2,1,FACEDATA_B)
                  : cell[0][1][0]->celldata[CELLDATA_B][2];
            bf3 = (cell[0][1][1]!=cell[1][1][1]) ? cell[0][1][1]->faceave(0,1,FACEDATA_B)
                  : cell[0][1][1]->celldata[CELLDATA_B][0];
            bf4 = (cell[1][1][1]!=cell[1][1][0]) ? cell[1][1][1]->faceave(2,0,FACEDATA_B)
                  : cell[1][1][1]->celldata[CELLDATA_B][2];
            epos = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        if(cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][0][1]
           or cell[1][0][1]!=cell[0][0][1]) { // negative side edge exists?
            negf=true;
            bf1 = (cell[0][0][0]!=cell[1][0][0]) ? cell[0][0][0]->faceave(0,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][0];
            bf2 = (cell[0][0][0]!=cell[0][0][1]) ? cell[0][0][0]->faceave(2,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][2];
            bf3 = (cell[0][0][1]!=cell[1][0][1]) ? cell[0][0][1]->faceave(0,1,FACEDATA_B)
                  : cell[0][0][1]->celldata[CELLDATA_B][0];
            bf4 = (cell[1][0][1]!=cell[1][0][0]) ? cell[1][0][1]->faceave(2,0,FACEDATA_B)
                  : cell[1][0][1]->celldata[CELLDATA_B][2];
            eneg = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        break;
    case 2 :  // z-edges
        if(cell[0][0][1]!=cell[1][0][1] or cell[1][0][1]!=cell[1][1][1]
           or cell[1][1][1]!=cell[0][1][1]) { // positive side edge exists?
            posf=true;
            bf1 = (cell[0][0][1]!=cell[0][1][1]) ? cell[0][0][1]->faceave(1,1,FACEDATA_B)
                  : cell[0][0][1]->celldata[CELLDATA_B][1];
            bf2 = (cell[0][0][1]!=cell[1][0][1]) ? cell[0][0][1]->faceave(0,1,FACEDATA_B)
                  : cell[0][0][1]->celldata[CELLDATA_B][0];
            bf3 = (cell[1][1][1]!=cell[1][0][1]) ? cell[1][0][1]->faceave(1,1,FACEDATA_B)
                  : cell[1][0][1]->celldata[CELLDATA_B][1];
            bf4 = (cell[1][1][1]!=cell[0][1][1]) ? cell[1][1][1]->faceave(0,0,FACEDATA_B)
                  : cell[1][1][1]->celldata[CELLDATA_B][0];
            epos = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        if(cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][1][0]
           or cell[1][1][0]!=cell[0][1][0]) { // negative side edge exists?
            negf=true;
            bf1 = (cell[0][0][0]!=cell[0][1][0]) ? cell[0][0][0]->faceave(1,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][1];
            bf2 = (cell[0][0][0]!=cell[1][0][0]) ? cell[0][0][0]->faceave(0,1,FACEDATA_B)
                  : cell[0][0][0]->celldata[CELLDATA_B][0];
            bf3 = (cell[1][1][0]!=cell[1][0][0]) ? cell[1][0][0]->faceave(1,1,FACEDATA_B)
                  : cell[1][0][0]->celldata[CELLDATA_B][1];
            bf4 = (cell[1][1][0]!=cell[0][1][0]) ? cell[1][1][0]->faceave(0,0,FACEDATA_B)
                  : cell[1][1][0]->celldata[CELLDATA_B][0];
            eneg = invdxmu0*(-bf1+bf2+bf3-bf4);
        }
        break;
    default :
        cerr << "Error (FATAL): jcomponent called with invalid direction." <<endl;
        exit(-1);
    }
    intp = (posf and negf) ? 0.5 : 1.0;
    return intp*(epos+eneg);
}

//! Calculate J
void Tgrid::Tnode::calc_j(real dx)
{
    nodedata[NODEDATA_J][0] = jcomponent(0,dx);
    nodedata[NODEDATA_J][1] = jcomponent(1,dx);
    nodedata[NODEDATA_J][2] = jcomponent(2,dx);
    real vw = 0.0;
    const real rho_q = nodedata[NODEDATA_ne][0];
    const real Btot = sqrt( sqr(nodedata[NODEDATA_B][0]) + sqr(nodedata[NODEDATA_B][1]) + sqr(nodedata[NODEDATA_B][2]) );
    if(rho_q > 0.0 && dx > 0.0) {
       vw = 2.0*Btot*M_PI/( Params::mu_0*rho_q*dx );
    }
    if (vw > Params::maxVw && Params::maxVw > 0.0) {
       const real d = Params::maxVw/vw;
       nodedata[NODEDATA_J][0] *= d;
       nodedata[NODEDATA_J][1] *= d;
       nodedata[NODEDATA_J][2] *= d;
#ifndef NO_DIAGNOSTICS
       Tgrid::fieldCounter.cutRateMaxVw += 1.0;
#endif
    }
}

//! Calculate J at a node (recursive)
void Tgrid::Tcell::calc_node_j_recursive()
{
    int f,f2,ch,d=0;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->calc_node_j_recursive();
    else {
        // NODE LOOP
        if (isrefined_face(d,1))
            for(f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_j(size);
        else
            face[d][1]->node[2]->calc_j(size);
        for (d=1; d<3; d++)
            if (isrefined_face(d,1))
                for (f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_j(size);
    }
}

//! Calculate J at a node
void Tgrid::calc_node_j()
{
    int i,j,k;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++)
                cells[flatindex(i,j,k)]->calc_node_j_recursive();
}

//! cell2node interpolation of ne
void Tgrid::Tnode::CN1_ne()
{
    real ne=0.;
    register Tcell *c;
    gridreal weightsum = 0.0;
    gridreal invweightsum;
    for (int a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        register const gridreal weight = c->invsize;
        ne+=weight*c->rho_q;
        weightsum+= weight;
    }
    invweightsum = (weightsum==0.) ? 0. : 1./weightsum;
    nodedata[NODEDATA_ne][0]=nodedata[NODEDATA_ne][1]=nodedata[NODEDATA_ne][2]=ne*invweightsum;
}

//! cell2node interpolation of ne (recursive)
void Tgrid::Tcell::CN_ne_recursive()
{
    int ch,d,f,f2;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->CN_ne_recursive();
    else {
        // NODE LOOP
        d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1))
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_ne();
        else
            face[d][1]->node[2]->CN1_ne(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        for (d=1; d<3; d++) if (isrefined_face(d,1))
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->CN1_ne();
    }
}

//! cell2node interpolation of ne
void Tgrid::CN_ne()
{
    int i,j,k;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++)
                cells[flatindex(i,j,k)]->CN_ne_recursive();
}

//! Calculate Ue at a node
void Tgrid::Tnode::calc_ue1()
{
    real rhoq =  nodedata[NODEDATA_ne][0];
    real invrhoq = (rhoq > 0) ? 1./rhoq : 0.;
    // Check zero field radius
    if (r2 < Params::R_zeroFields2) {
        for(int d=0; d<3; d++)
            nodedata[NODEDATA_UE][d] =0.;
        return;
    }
    real ue2=0.;
    for(int d=0; d<3; d++) {
#ifndef IGNORE_ELECTRIC_FIELD_HALL_TERM
        nodedata[NODEDATA_UE][d] = (nodedata[NODEDATA_Ji][d] - nodedata[NODEDATA_J][d])*invrhoq;
#else
        nodedata[NODEDATA_UE][d] = (nodedata[NODEDATA_Ji][d])*invrhoq;
#endif
        ue2+= sqr(nodedata[NODEDATA_UE][d]);
    }
    // Check maximum electron fluid velocity (CONSTRAINT)
    if (Params::Ue_max > 0 && ue2 > Params::Ue_max2) {
        const real norm = Params::Ue_max/sqrt(ue2);
        for(int d=0; d<3; d++)
            nodedata[NODEDATA_UE][d] *= norm;
#ifndef NO_DIAGNOSTICS
        Tgrid::fieldCounter.cutRateUe += 1.0;
#endif
    }
}

//! Calculate Ue at a node (recursive)
void Tgrid::Tcell::calc_node_ue_recursive()
{
    int ch,d,f,f2;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->calc_node_ue_recursive();
    else {
        // NODE LOOP
        d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1))
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_ue1();
        else
            face[d][1]->node[2]->calc_ue1(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        for (d=1; d<3; d++) if (isrefined_face(d,1))
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->calc_ue1();
    }
}

//! Calculate Ue at a node
void Tgrid::calc_node_ue()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->calc_node_ue_recursive();
            }
}

// =================================================================================
// ========================== SPHERICAL COORDINATE SYSTEM ==========================
// =================================================================================

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of Tgrid::init
void Tgrid::sph_init(int nx1, int ny1, int nz1, gridreal x1, gridreal y1, gridreal z1, gridreal theta1, gridreal phi1, gridreal bgdx1, gridreal bgdy1, gridreal bgdz1)
{
    MSGFUNCTIONCALL("Tgrid::sph_init");
    static bool onlyOneTgridObject = false;
    if (onlyOneTgridObject == true) {
        ERRORMSG("only one Tgrid object can be defined");
        doabort();
    } else {
        onlyOneTgridObject = true;
    }
    saved_cellptr = 0;
    n_particles = 0;
    n_pdftables = 0;
    ave_ntimes = 0;
    previous_found_cell = 0;
    nx = nx1 + 2;
    ny = ny1 + 2;
    nz = nz1 + 2;
    bgdx     = bgdx1;
    sph_bgdy = bgdy1;
    sph_bgdz = bgdz1;
    invbgdx = 1.0/bgdx;
    sph_invbgdy = 1.0/sph_bgdy;
    sph_invbgdz = 1.0/sph_bgdz;
    x_1     = x1 - bgdx;
    y_1     = y1 - sph_bgdy;
    z_1     = z1 - sph_bgdz;
    if (Params::sph_BC_type == 0) {
        sph_bgdtheta = Params::sph_box_theta/ny;
        sph_bgdphi   = Params::sph_box_phi/nz;
        sph_theta_1 = theta1;
        sph_phi_1   = phi1;
    }
    if (Params::sph_BC_type == 1) {
        sph_bgdtheta = Params::sph_dtheta;
        sph_bgdphi   = Params::sph_dphi;
        sph_theta_1 = theta1 - sph_bgdtheta;
        sph_phi_1   = phi1 - sph_bgdphi;
    }
    // "unit" is the smallest safely representable distance using gridreal
    inv_unit = invbgdx*65536.0;
    const int N = nx*ny*nz;
    mainlog << "|------------------------- GRID DETAILS -------------------------|\n";
    mainlog << "| Grid [nx,ny,nz]  = [" << nx << "," << ny << "," << nz << "] (ghosts included)\n"
            << "| Base cells       = " << N << " (with ghosts)\n"
            << "| Base cells       = " << nx1*ny1*nz1 << " (no ghosts)\n"
            << "| Base cell size   = " << bgdx/1e3 << " km = " << bgdx/Params::R_P << " R_P\n"
            << "| Box size [x,y,z] = [" << Params::box_X/1e3 << ", " << Params::box_Y/1e3 << ", " << Params::box_Z/1e3 << "] km  = [" <<
            Params::box_X/Params::R_P << "," << Params::box_Y/Params::R_P << "," << Params::box_Z/Params::R_P << "] R_P\n"
            << "| Box volume       = " << Params::box_V << " m^3 = " << Params::box_V/cube(Params::R_P) << " R_P^3\n"
            << "| Gridreal unit    = " << 1.0/inv_unit << " m\n"
            << "| Bytes/cell       ~ " << approx_bytes_per_cell() << "\n";
    mainlog << "|----------------------------------------------------------------|\n\n";
    mainlog << "|------------------------ GRID DETAILS II -----------------------|\n";
    mainlog << "| Base x cell size = " << bgdx/Params::R_P << " R_P\n"
            << "| x_1              = " << x_1/Params::R_P << " R_P\n"
            << "| x(nx)            = " << (x_1 + (nx)*bgdx)/Params::R_P << " R_P\n"
            << "| x(1)             = " << (x_1 + 1*bgdx)/Params::R_P << " R_P\n"
            << "| x(nx-1)          = " << (x_1 + (nx-1)*bgdx)/Params::R_P << " R_P\n";
    mainlog << "| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|\n";
    mainlog << "| Base y cell size = " << sph_bgdy/Params::R_P << " R_P\n"
            << "| y_1              = " << y_1/Params::R_P << " R_P\n"
            << "| y(ny)            = " << (y_1 + (ny)*sph_bgdy)/Params::R_P << " R_P\n"
            << "| y(1)             = " << (y_1 + 1*sph_bgdy)/Params::R_P << " R_P\n"
            << "| y(ny-1)          = " << (y_1 + (ny-1)*sph_bgdy)/Params::R_P << " R_P\n";
    mainlog << "| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|\n";
    mainlog << "| Base y cell size = " << sph_bgdz/Params::R_P << " R_P\n"
            << "| z_1              = " << z_1/Params::R_P << " R_P\n"
            << "| z(nz)            = " << (z_1 + (nz)*sph_bgdz)/Params::R_P << " R_P\n"
            << "| z(1)             = " << (z_1 + 1*sph_bgdz)/Params::R_P << " R_P\n"
            << "| z(nz-1)          = " << (z_1 + (nz-1)*sph_bgdz)/Params::R_P << " R_P\n";
    mainlog << "|----------------------------------------------------------------|\n\n";
    mainlog << "|-------------------- SPHERICAL GRID DETAILS --------------------|\n";
    mainlog << "| Theta cell size  = " << sph_bgdtheta << "\n"
            << "| Theta box        = " << Params::sph_box_theta << "\n"
            << "| sph_theta_1          = " << sph_theta_1 << "\n"
            << "| theta(1)         = " << sph_theta_1 + 1*sph_bgdtheta << "\n"
            << "| theta(ny-1)      = " << sph_theta_1 + (ny-1)*sph_bgdtheta << "\n";
    mainlog << "| - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -|\n";
    mainlog << "| Phi cell size    = " << sph_bgdphi << "\n"
            << "| Phi box          = " << Params::sph_box_phi << "\n"
            << "| sph_phi_1            = " << sph_phi_1 << "\n"
            << "| phi(1)           = " << sph_phi_1 + 1*sph_bgdphi << "\n"
            << "| phi(nz-1)        = " << sph_phi_1 + (nz-1)*sph_bgdphi << "\n";
    mainlog << "|----------------------------------------------------------------|\n\n";
    mainlog << "|--------------- COORDINATE LIMITS ---------------|\n";
    mainlog << "| " << Params::box_xmin/1e3 << " km (" << Params::box_xmin/Params::R_P << " R_P) < x < " << Params::box_xmax/1e3 << " km (" <<
            Params::box_xmax/Params::R_P << " R_P)\n"
            << "| " << Params::box_ymin/1e3 << " km (" << Params::box_ymin/Params::R_P << " R_P) < y < " << Params::box_ymax/1e3 << " km (" << Params::box_ymax/Params::R_P << " R_P)\n"
            << "| " << Params::box_zmin/1e3 << " km (" << Params::box_zmin/Params::R_P << " R_P) < z < " << Params::box_zmax/1e3 << " km (" << Params::box_zmax/Params::R_P << " R_P)\n"
            << "| " << Params::sph_theta_min/pi << " pi < theta < " << Params::sph_theta_max/pi << " pi\n"
            << "| " << Params::sph_phi_min/pi << " pi < phi < " << Params::sph_phi_max/pi << " pi\n";
    mainlog << "|-------------------------------------------------|\n\n";
    cells = new TCellPtr [N];
    int i,j,k,c;
    // Create cells and faces, set up cell fields but not yet face fields
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c] = new Tcell;
        cells[c]->plist.init();
        cells[c]->haschildren = false;
        cells[c]->refine_it = cells[c]->recoarsen_it = cells[c]->forbid_psplit = false;
        cells[c]->refstatus = 0;
        cells[c]->level = 0;
        cells[c]->parent = 0;
        cells[c]->running_index = -123456;
        cells[c]->face[0][0] = (i > 0) ? cells[flatindex(i-1,j,k)]->face[0][1] : TFacePtr(0);
        cells[c]->face[0][1] = (i < nx-1) ? new Tface : 0;
        cells[c]->face[1][0] = (j > 0) ? cells[flatindex(i,j-1,k)]->face[1][1] : TFacePtr(0);
        cells[c]->face[1][1] = (j < ny-1) ? new Tface : 0;
        cells[c]->face[2][0] = (k > 0) ? cells[flatindex(i,j,k-1)]->face[2][1] : TFacePtr(0);
        cells[c]->face[2][1] = (k < nz-1) ? new Tface : 0;
        cells[c]->nc = 0.0;
        cells[c]->rho_q = 0.0;
        cells[c]->rho_q_bg = 0.0;
        cells[c]->ave_nc = 0.0;
#ifdef SAVE_POPULATION_AVERAGES
        for (int l = 0; l < Params::POPULATIONS; ++l) {
            cells[c]->pop_ave_n.push_back(0.0);
            cells[c]->pop_ave_vx.push_back(0.0);
            cells[c]->pop_ave_vy.push_back(0.0);
            cells[c]->pop_ave_vz.push_back(0.0);
        }
#endif
        memset(&cells[c]->celldata[0][0],0,sizeof(datareal)*NCELLDATA*3);
        cells[c]->centroid[0] = x_1 + (i+0.5)*bgdx;
        cells[c]->centroid[1] = y_1 + (j+0.5)*sph_bgdy;
        cells[c]->centroid[2] = z_1 + (k+0.5)*sph_bgdz;
        cells[c]->size  = bgdx;
        cells[c]->sph_sizey = sph_bgdy;
        cells[c]->sph_sizez = sph_bgdz;
        cells[c]->r2 = sqr(cells[c]->centroid[0]);
        cells[c]->invsize = 1.0/bgdx;
        cells[c]->sph_invsizey = 1.0/sph_bgdy;
        cells[c]->sph_invsizez = 1.0/sph_bgdz;
        cells[c]->flatind = c;
        cells[c]->sph_CN_boundary_flag = 0;
        // SC coordinate elements!!
        cells[c]->sph_centroid[0] = x_1     + (i+0.5)*bgdx;
        cells[c]->sph_centroid[1] = sph_theta_1 + (j+0.5)*sph_bgdtheta;
        cells[c]->sph_centroid[2] = sph_phi_1   + (k+0.5)*sph_bgdphi;
        cells[c]->sph_coor[0] = x_1     + i*bgdx;
        cells[c]->sph_coor[1] = sph_theta_1 + j*sph_bgdtheta;
        cells[c]->sph_coor[2] = sph_phi_1   + k*sph_bgdphi;
        ///////////////////////////////////////////////////////////////////////////////////////////
        // Diagnostic tool. To check the limits for shpherical coordinates theta and phi //
        ///////////////////////////////////////////////////////////////////////////////////////////
        /*if (cells[c]->sph_coor[1] < 0.0 || cells[c]->sph_coor[1] > pi)
           {
        errorlog << "j = " << j << " theta(j) = " << sph_theta_1 + j*bgdtheta << " sph_theta_1 = " << sph_theta_1 << " bgdtheta = " << bgdtheta << " ny = " << ny << " theta_max = " << sph_theta_1 + ny*bgdtheta << "\n";
            doabort();
        }
            if (cells[c]->sph_coor[2] < -pi - 1e-06 || cells[c]->sph_coor[2] > pi + 1e-06)
               {
            	  errorlog << "k = " << k << " phi(k) = " << sph_phi_1 + k*bgdphi << " sph_phi_1 = " << sph_phi_1 << " bgdphi = " << bgdphi << " nz = " << nz << " phi_max = " << sph_phi_1 + ny*bgdphi << "\n";
                doabort();
         }*/
        ///////////////////////////////////////////////////////////////////////////////////////////
        /*if (cells[c]->sph_coor[1] == 0.1 || cells[c]->sph_coor[1] == pi)
           {
        errorlog << "j = " << j << " theta(j) = " << cells[c]->sph_coor[1] << "\n";
            doabort();
        }
            //if (cells[c]->sph_coor[2] == -pi || cells[c]->sph_coor[2] >= pi)
            if (k == 19)
               {
            	  errorlog << "k = " << k << " phi(k) = " << cells[c]->sph_coor[2] << " bgdphi = " << bgdphi << " nz1 = " << nz1 << " nz = " << nz << "\n";
                doabort();
         }*/
        ///////////////////////////////////////////////////////////////////////////////////////////
        cells[c]->sph_coor_next[0] = x_1     + (i+1)*bgdx;     // r
        cells[c]->sph_coor_next[1] = sph_theta_1 + (j+1)*sph_bgdtheta; // theta
        cells[c]->sph_coor_next[2] = sph_phi_1   + (k+1)*sph_bgdphi;   // phi
        cells[c]->sph_size[0] = bgdx;     // dr
        cells[c]->sph_size[1] = sph_bgdtheta; // dtheta
        cells[c]->sph_size[2] = sph_bgdphi;   // dphi
        // SC line grid elements!!
        cells[c]->sph_dl[0] = cells[c]->sph_coor_next[0] - cells[c]->sph_coor[0];                                                    // dl_r
        cells[c]->sph_dl[1] = cells[c]->sph_coor[0]*(cells[c]->sph_coor_next[1] - cells[c]->sph_coor[1]);                            // dl_theta
        cells[c]->sph_dl[2] = cells[c]->sph_coor[0]*sin(cells[c]->sph_coor[1])*(cells[c]->sph_coor_next[2] - cells[c]->sph_coor[2]); // dl_phi
        cells[c]->sph_diag  = sqrt(sqr(cells[c]->sph_dl[0]) + sqr(cells[c]->sph_dl[1]) + sqr(cells[c]->sph_dl[2]));                  // diagonal
        // SC areas!!
        // Area of cell dS _|_ dr dS(0)=r^2*sin(theta)*dtheta*dphi
        cells[c]->sph_dS[0]=sqr(cells[c]->sph_coor[0])*sin(cells[c]->sph_coor[1])*cells[c]->sph_size[1]*cells[c]->sph_size[2];
        // Area of cell dS _|_ dtheta dS(1)=rc*sin(theta)*dr*dphi
        cells[c]->sph_dS[1]=cells[c]->sph_centroid[0]*sin(cells[c]->sph_coor[1])*cells[c]->sph_size[0]*cells[c]->sph_size[2];
        // Area of cell dS _|_ dphi dS(2)=rc*dr*dtheta
        cells[c]->sph_dS[2]=cells[c]->sph_centroid[0]*cells[c]->sph_size[0]*cells[c]->sph_size[1];
        // Area of cell dS _|_ dr dS(0)=r_c^2*sin(theta)*dtheta*dphi
        cells[c]->sph_dS_centroid[0]=sqr(cells[c]->sph_centroid[0])*sin(cells[c]->sph_coor[1])*cells[c]->sph_size[1]*cells[c]->sph_size[2];
        // Area of cell dS _|_ dtheta dS(1)=rc*sin(theta_c)*dr*dphi
        cells[c]->sph_dS_centroid[1]=cells[c]->sph_centroid[0]*sin(cells[c]->sph_centroid[1])*cells[c]->sph_size[0]*cells[c]->sph_size[2];
        // Area of cell dS _|_ dphi dS(2)=rc*dr*dtheta
        cells[c]->sph_dS_centroid[2]=cells[c]->sph_centroid[0]*cells[c]->sph_size[0]*cells[c]->sph_size[1];
        // Area of next cell  dS _|_ dr dS(0)=r_next^2*sin(theta)*dtheta*dphi
        cells[c]->sph_dS_next[0]=sqr(cells[c]->sph_coor_next[0])*sin(cells[c]->sph_coor[1])*cells[c]->sph_size[1]*cells[c]->sph_size[2];
        // Area of next cell dS _|_ dtheta dS(1)=rc*sin(theta_next)*dr*dphi
        cells[c]->sph_dS_next[1]=cells[c]->sph_centroid[0]*sin(cells[c]->sph_coor_next[1])*cells[c]->sph_size[0]*cells[c]->sph_size[2];
        // Area of next cell dS _|_ dphi dS(2)=rc*dr*dtheta
        cells[c]->sph_dS_next[2]=cells[c]->sph_centroid[0]*cells[c]->sph_size[0]*cells[c]->sph_size[1];
        // SC volumes!!
        // Volume of cell dV=r^2*sin(theta)*dr*dtheta*dphi
        cells[c]->sph_dV = sqr(cells[c]->sph_centroid[0])*sin(cells[c]->sph_centroid[1])*cells[c]->sph_size[0]*cells[c]->sph_size[1]*cells[c]->sph_size[2];
        // Volume of cell dV=dphi*(r2^3-r1^3)*(cos(theta1)-cos(theta2))/3
        // cells[c]->sph_dV = cells[c]->sph_size[2]*(pow(cells[c]->sph_coor_next[0], 3) - pow(cells[c]->sph_coor[0], 3))
        //                    *(cos(cells[c]->sph_coor[1]) - cos(cells[c]->sph_coor_next[1]))/3.0;
    }

    // Set up neighbour fields
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->neighbour[0][0] = (i > 0)    ? cells[flatindex(i-1,j,k)] : TCellPtr(0);
        cells[c]->neighbour[0][1] = (i < nx-1) ? cells[flatindex(i+1,j,k)] : TCellPtr(0);
        cells[c]->neighbour[1][0] = (j > 0)    ? cells[flatindex(i,j-1,k)] : TCellPtr(0);
        cells[c]->neighbour[1][1] = (j < ny-1) ? cells[flatindex(i,j+1,k)] : TCellPtr(0);
        cells[c]->neighbour[2][0] = (k > 0)    ? cells[flatindex(i,j,k-1)] : TCellPtr(0);
        cells[c]->neighbour[2][1] = (k < nz-1) ? cells[flatindex(i,j,k+1)] : TCellPtr(0);
    }
    // Create nodes, set up node fields.
    // Must store the nodes somewhere, use temporary vector 'nodes' for it.
    TNodePtr *nodes = new TNodePtr [N];
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                const TNodePtr n = new Tnode;
                n->cell[0][0][0] = cells[c];
                n->cell[1][0][0] = cells[flatindex(i+1,j,  k  )];
                n->cell[0][1][0] = cells[flatindex(i,  j+1,k  )];
                n->cell[1][1][0] = cells[flatindex(i+1,j+1,k  )];
                n->cell[0][0][1] = cells[flatindex(i,  j,  k+1)];
                n->cell[1][0][1] = cells[flatindex(i+1,j,  k+1)];
                n->cell[0][1][1] = cells[flatindex(i,  j+1,k+1)];
                n->cell[1][1][1] = cells[flatindex(i+1,j+1,k+1)];
                n->centroid[0] = x_1 + (i+1)*bgdx;
                n->centroid[1] = y_1 + (j+1)*sph_bgdy;
                n->centroid[2] = z_1 + (k+1)*sph_bgdz;
                // SC coordinate elements!!
                n->sph_centroid[0] = x_1     + (i+1)*bgdx;
                n->sph_centroid[1] = sph_theta_1 + (j+1)*sph_bgdtheta;
                n->sph_centroid[2] = sph_phi_1   + (k+1)*sph_bgdphi;
                n->r0 = n->sph_centroid[0];
                nodes[c] = n;
            }
    // Set up pointers from faces to nodes. Use right-pointing faces (dir=1)
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                // X-directed faces:
                if (j > 0 && k > 0) cells[c]->face[0][1]->node[0] = nodes[flatindex(i,j-1,k-1)];
                if (k > 0         ) cells[c]->face[0][1]->node[1] = nodes[flatindex(i,j,  k-1)];
                if (true          ) cells[c]->face[0][1]->node[2] = nodes[flatindex(i,j,  k)];
                if (j > 0         ) cells[c]->face[0][1]->node[3] = nodes[flatindex(i,j-1,k)];
                // Y-directed faces:
                if (i > 0 && k > 0) cells[c]->face[1][1]->node[0] = nodes[flatindex(i-1,j,k-1)];
                if (i > 0         ) cells[c]->face[1][1]->node[1] = nodes[flatindex(i-1,j,k  )];
                if (true          ) cells[c]->face[1][1]->node[2] = nodes[flatindex(i,  j,k)];
                if (k > 0         ) cells[c]->face[1][1]->node[3] = nodes[flatindex(i  ,j,k-1)];
                // Z-directed faces:
                if (i > 0 && j > 0) cells[c]->face[2][1]->node[0] = nodes[flatindex(i-1,j-1,k)];
                if (j > 0         ) cells[c]->face[2][1]->node[1] = nodes[flatindex(i,  j-1,k)];
                if (true          ) cells[c]->face[2][1]->node[2] = nodes[flatindex(i,  j,  k)];
                if (i > 0         ) cells[c]->face[2][1]->node[3] = nodes[flatindex(i-1,j,  k)];
            }
    delete [] nodes;
    MSGFUNCTIONEND("Tgrid::sph_init");
}

//! (SPHERICAL) Spherical version of "FC_recursive"
void Tgrid::Tcell::sph_FC_recursive(TFaceDataSelect fs, TCellDataSelect cs, int FC_boundary_flag)
{
    int d;
    gridreal AF0[3] = {faceave(0,0,fs), faceave(1,0,fs), faceave(2,0,fs)};
    gridreal AF1[3] = {faceave(0,1,fs), faceave(1,1,fs), faceave(2,1,fs)};
    gridreal Ac[3]  = {0, 0, 0};
    gridreal dS0[3] = {sph_dS[0], sph_dS[1], sph_dS[2]};
    gridreal dS1[3] = {sph_dS_next[0], sph_dS_next[1], sph_dS_next[2]};
    gridreal dSc[3] = {sph_dS_centroid[0], sph_dS_centroid[1], sph_dS_centroid[2]};
    if (FC_boundary_flag == 0) {
        //for (d=0; d<3; d++) Ac[d] = 0.5*(AF0[d] + AF1[d]);
        for (d=0; d<3; d++) Ac[d] = 0.5*(AF0[d]*dS0[d] + AF1[d]*dS1[d])/dSc[d];
    }
    if (FC_boundary_flag == 20) {
        //Ac[0] = 0.5*(AF0[0] + AF1[0]);
        //Ac[1] =      AF1[1];
        //Ac[2] = 0.5*(AF0[2] + AF1[2]);
        Ac[0] = 0.5*(AF0[0]*dS0[0] + AF1[0]*dS1[0])/dSc[0];
        Ac[1] =      AF1[1];
        Ac[2] = 0.5*(AF0[2]*dS0[2] + AF1[2]*dS1[2])/dSc[2];
    }
    if (FC_boundary_flag == 21) {
        //Ac[0] = 0.5*(AF0[0] + AF1[0]);
        //Ac[1] =      AF0[1];
        //Ac[2] = 0.5*(AF0[2] + AF1[2]);
        Ac[0] = 0.5*(AF0[0]*dS0[0] + AF1[0]*dS1[0])/dSc[0];
        Ac[1] =      AF0[1];
        Ac[2] = 0.5*(AF0[2]*dS0[2] + AF1[2]*dS1[2])/dSc[2];
    }
    for (d=0; d<3; d++) celldata[cs][d] = Ac[d];
}

//! (SPHERICAL) Spherical version of "NC_recursive"
void Tgrid::Tcell::sph_NC_recursive(TNodeDataSelect ns,TCellDataSelect cs)
{
    int dir,d,f,f2;
    gridreal rn[3];
    gridreal A[3];
    gridreal temp_car[3] = {0.0, 0.0, 0.0};
    gridreal temp_sph[3] = {0.0, 0.0, 0.0};
    celldata[cs][0]=0.;
    celldata[cs][1]=0.;
    celldata[cs][2]=0.;
    for (dir=0; dir<3; dir++) for(d=0; d<2; d++) {
            temp_car[0]=temp_car[1]=temp_car[2]=0.;
            // without refinement we are interested only this block
            for (f=0; f<4; f++) {
                rn[0] = face[dir][d]->node[f]->sph_centroid[0];
                rn[1] = face[dir][d]->node[f]->sph_centroid[1];
                rn[2] = face[dir][d]->node[f]->sph_centroid[2];
                A[0] = face[dir][d]->node[f]->nodedata[ns][0];
                A[1] = face[dir][d]->node[f]->nodedata[ns][1];
                A[2] = face[dir][d]->node[f]->nodedata[ns][2];
                sph_transf_S2C_A1(rn, A);
                temp_car[0] += A[0];
                temp_car[1] += A[1];
                temp_car[2] += A[2];
            }
            // 6 faces * 4 nodes from each face = 24
            temp_sph[0] += temp_car[0]/24;
            temp_sph[1] += temp_car[1]/24;
            temp_sph[2] += temp_car[2]/24;
        }
    gridreal rc[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    sph_transf_C2S_A1(rc, temp_sph);
    celldata[cs][0] = temp_sph[0];
    celldata[cs][1] = temp_sph[1];
    celldata[cs][2] = temp_sph[2];
}

/*//! (SPHERICAL) Another Spherical version of "NC_recursive"
void Tgrid::Tcell::sph_NC_recursive(TNodeDataSelect ns,TCellDataSelect cs)
{
  if (haschildren) {
    int ch;
    for (ch=0; ch<8; ch++) child[0][0][ch]->sph_NC_recursive(ns,cs);
   }
else
   {
    int dir,d,f,f2;
    datareal tempx,tempy,tempz;
    gridreal weight = 0.0;
    gridreal weightsum = 0.0;
    celldata[cs][0]=0.;
    celldata[cs][1]=0.;
    celldata[cs][2]=0.;
    gridreal dS[3][2];
             dS[0][0] = sph_dS[0];
             dS[0][1] = sph_dS_next[0];
             dS[1][0] = sph_dS[1];
             dS[1][1] = sph_dS_next[1];
             dS[2][0] = sph_dS[2];
             dS[2][1] = sph_dS_next[2];

    for (dir=0; dir<3; dir++) for(d=0; d<2; d++)
      {
       tempx=tempy=tempz=0.;
       if (isrefined_face(dir,d))
        {
	 for (f=0; f<4; f++) for (f2=0; f2<4; f2++)
           {
            tempx+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][0];
            tempy+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][1];
            tempz+=refintf[dir][d]->face[f]->node[f2]->nodedata[ns][2];
	   }
	 tempx*=0.25;
	 tempy*=0.25;
	 tempz*=0.25;
	}
     else
        {
         // without refinement we are interested only this block
	 for (f=0;f<4;f++)
           {
            weight = dS[dir][d];

            tempx +=face[dir][d]->node[f]->nodedata[ns][0]*weight;
	    tempy +=face[dir][d]->node[f]->nodedata[ns][1]*weight;
	    tempz +=face[dir][d]->node[f]->nodedata[ns][2]*weight;

            weightsum+= weight;
	   }
	}
        // 6 faces * 4 nodes from each face = 24
        register const gridreal invweightsum = 1.0/weightsum;

	celldata[cs][0]+=tempx*invweightsum;
	celldata[cs][1]+=tempy*invweightsum;
	celldata[cs][2]+=tempz*invweightsum;
       }
  }
}*/

//! (SPHERICAL) Spherical version of "NC_smoothing_recursive"
void Tgrid::Tcell::sph_NC_smoothing_recursive()
{
    int dir,d,f,f2;
    datareal tempnc,temprhoq,tempvx,tempvy,tempvz;
    nc = 0.0;
    rho_q = 0.0;
    celldata[CELLDATA_Ji][0] = 0.0;
    celldata[CELLDATA_Ji][1] = 0.0;
    celldata[CELLDATA_Ji][2] = 0.0;
    for (dir=0; dir<3; dir++) for(d=0; d<2; d++) {
            tempnc = temprhoq = tempvx = tempvy = tempvz = 0.;
            for (f=0; f<4; f++) {
                tempnc   +=face[dir][d]->node[f]->nn;
                temprhoq +=face[dir][d]->node[f]->nodedata[NODEDATA_UE][0];
                // For vector values we need to use Cartesian coordinates for averaging
                // Ji is calculated in Cartesian coordinates. We don't need to make any transformation
                gridreal r[3] = {face[dir][d]->node[f]->sph_centroid[0], face[dir][d]->node[f]->sph_centroid[1],face[dir][d]->node[f]->sph_centroid[2]};
                gridreal A[3] = {celldata[CELLDATA_J][0], celldata[CELLDATA_J][1], celldata[CELLDATA_J][2]};
                tempvx   +=face[dir][d]->node[f]->nodedata[NODEDATA_J][0];
                tempvy   +=face[dir][d]->node[f]->nodedata[NODEDATA_J][1];
                tempvz   +=face[dir][d]->node[f]->nodedata[NODEDATA_J][2];
            }
            // 6 faces * 4 nodes from each face = 24
            nc+=tempnc/24.;
            rho_q+=temprhoq/24.;
            celldata[CELLDATA_Ji][0]+=tempvx/24.;
            celldata[CELLDATA_Ji][1]+=tempvy/24.;
            celldata[CELLDATA_Ji][2]+=tempvz/24.;
        }
}

//! (SPHERICAL) Spherical version of "CN1". (Very simple version, Based on Cartesian transformation.)
void Tgrid::Tnode::sph_CN1(TCellDataSelect cs, TNodeDataSelect ns)
{
    int a;
    gridreal result[3] = {0,0,0};
    gridreal temp[3];
    gridreal rn[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
        gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        sph_transf_S2C_A1(r, A);
        for (d=0; d<3; d++) result[d]+= A[d];
    }
    for (d=0; d<3; d++) temp[d] = result[d]/8;
    sph_transf_C2S_A1(rn, temp);
    for (d=0; d<3; d++) nodedata[ns][d] = temp[d];
}

/*
//! (SPHERICAL) 2nd Spherical version of "CN1"
void Tgrid::Tnode::sph_CN1(TCellDataSelect cs, TNodeDataSelect ns)
{
  int a;
  real flux[3] = {0,0,0};
  gridreal temp[3];

  gridreal areax = 0.0;
  gridreal areay = 0.0;
  gridreal areaz = 0.0;

  gridreal area_nodex = 0.0;
  gridreal area_nodey = 0.0;
  gridreal area_nodez = 0.0;

  gridreal total_areax = 0.0;
  gridreal total_areay = 0.0;
  gridreal total_areaz = 0.0;

  register Tcell *c;

  // VERY IMPORTANT
  // The scheme of cells neighboured by node

  // cell[0][0][0] = cell[0][0][0]
  // cell[0][0][1] = cell[0][0][1]
  // cell[0][1][0] = cell[0][0][2]
  // cell[0][1][1] = cell[0][0][3]
  // cell[1][0][0] = cell[0][0][4]
  // cell[1][0][1] = cell[0][0][5]
  // cell[1][1][0] = cell[0][0][6]
  // cell[1][1][1] = cell[0][0][7]

  // Calculation a magnetic flux thruogh the certain cell center in different directions (x,y,z) and accumulation of the total fluxes in result[d]
  for (a=0; a<8; a++)
      {
       c = cell[0][0][a];
       if (c == 0) continue;
       total_areax = abs(c->sph_dS_centroid[0]);
       total_areay = abs(c->sph_dS_centroid[1]);
       total_areaz = abs(c->sph_dS_centroid[2]);

       gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
       gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};

       //sph_transf_H2S_V(A);
       //sph_transf_S2C_V(r, A);

       flux[0]+= total_areax*A[0];
       flux[1]+= total_areay*A[1];
       flux[2]+= total_areaz*A[2];
      }

  // Calculation of squares sum in different direction (x,y,z). See the scheme of cells.
  for (a=0; a<8; a++)
      {
       c = cell[0][0][a];
       if (c == 0) continue;
       areax = abs(c->sph_dS_next[0]);
       areay = abs(c->sph_dS_next[1]);
       areaz = abs(c->sph_dS_next[2]);

       if (a <= 3)
          {
           // Calculation of squares sum in x direction
           area_nodex+= areax;
          }
       if (a == 0 || a == 1 || a == 4 || a == 5)
          {
           // Calculation of squares sum in y direction
           area_nodey+= areay;
          }
       if (a == 0 || a == 2 || a == 4 || a == 6)
          {
           // Calculation of squares sum in z direction
           area_nodez+= areaz;
          }
      }

  // Calculation of magnetic field from total flux
  temp[0] = flux[0]/area_nodex/2;
  temp[1] = flux[1]/area_nodey/2;
  temp[2] = flux[2]/area_nodez/2;

  //gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};

  //sph_transf_C2S_V(r, temp);
  //sph_transf_S2H_V(temp);

  nodedata[ns][0] = temp[0];
  nodedata[ns][1] = temp[1];
  nodedata[ns][2] = temp[2];
}
*/

/*
//! (SPHERICAL) 3rd Spherical version of "CN1". Based on cell volume
void Tgrid::Tnode::sph_CN1_vol(TCellDataSelect cs, TNodeDataSelect ns)
{
  int a;
  real flux[3] = {0,0,0};

  gridreal volume = 0.0;
  gridreal total_volume = 0.0;

  register Tcell *c;

  // Calculation a magnetic flux thruogh the certain cell center in different directions (x,y,z) and accumulation of the total fluxes in result[d]
  for (a=0; a<8; a++)
      {
       c = cell[0][0][a];
       if (c == 0) continue;
       volume = abs(c->sph_dV);
       flux[0]+= volume*c->celldata[cs][0];
       flux[1]+= volume*c->celldata[cs][1];
       flux[2]+= volume*c->celldata[cs][2];

       total_volume += volume;
      }

  // Calculation of magnetic field from total flux
  nodedata[ns][0] = flux[0]/total_volume;
  nodedata[ns][1] = flux[1]/total_volume;
  nodedata[ns][2] = flux[2]/total_volume;
}
*/

//! (SPHERICAL) 4th Spherical version of "CN1". Use the cartesian vectors. Based on cell volume.
void Tgrid::Tnode::sph_CN1_vol(TCellDataSelect cs, TNodeDataSelect ns)
{
    int a;
    real flux[3] = {0,0,0};
    gridreal volume = 0.0;
    gridreal total_volume = 0.0;
    register Tcell *c;
    // Calculation a magnetic flux thruogh the certain cell center in different directions (x,y,z) and accumulation of the total fluxes in result[d]
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
        gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        sph_transf_H2S_V(A);
        sph_transf_S2C_V(r, A);
        volume = abs(c->sph_dV);
        flux[0]+= volume*A[0];
        flux[1]+= volume*A[1];
        flux[2]+= volume*A[2];
        total_volume += volume;
    }
    // Calculation of magnetic field from total flux
    nodedata[ns][0] = flux[0]/total_volume;
    nodedata[ns][1] = flux[1]/total_volume;
    nodedata[ns][2] = flux[2]/total_volume;
}

//! (SPHERICAL) Spherical version of "CN1" + CN boundary calculations.
void Tgrid::Tnode::sph_CNb1(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag)
{
    int a;
    gridreal result[3] = {0, 0, 0};
    gridreal temp[3];
    gridreal rn[3]     = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
        gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        sph_transf_S2C_A1(r, A);
        // 1 internal region
        if (sph_CN_boundary_flag == 1)  {
            for (d=0; d<3; d++) result[d]+= A[d];
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6)  {
            for (d=0; d<3; d++) result[d]+= 2*A[d];
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            for (d=0; d<3; d++) result[d]+= 4*A[d];
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8)  {
            for (d=0; d<3; d++) result[d]+= 8*A[d];
        }
    }
    for (d=0; d<3; d++) temp[d] = result[d]/8;
    sph_transf_C2S_A1(rn, temp);
    for (d=0; d<3; d++) nodedata[ns][d] = temp[d];
}

/*
//! (SPHERICAL) 2nd Spherical version of "CN1" + CN boundary calculations. Based on diagonal as a weight element
void Tgrid::Tnode::sph_CNb1(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag)
{
 int a;
 gridreal result[3] = {0, 0, 0};
 gridreal temp[3];
 gridreal rn[3]     = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
 gridreal weightsum = 0.0;

 int d;
 register Tcell *c;

 for (a=0; a<8; a++)
     {
      c = cell[0][0][a];
      if (c == 0) continue;

      gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
      gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};

      sph_transf_S2C_A1(r, A);

      // 1 internal region
      if (sph_CN_boundary_flag == 1)
         {
          gridreal weight = 1.0*c->sph_size[0];
          //gridreal weight = 1.0;

          result[0] += weight*A[0];
          result[1] += weight*A[1];
          result[2] += weight*A[2];
          weightsum += weight;
         }

      // 6 faces
      if (sph_CN_boundary_flag == 6)
         {
          gridreal weight = 2.0*c->sph_size[0];
          //gridreal weight = 2.0;

          result[0] += weight*A[0];
          result[1] += weight*A[1];
          result[2] += weight*A[2];
          weightsum += weight;
         }

      // 12 edges
      if (sph_CN_boundary_flag == 12)
         {
          gridreal weight = 4.0*c->sph_size[0];
          //gridreal weight = 4.0;

          result[0] += weight*A[0];
          result[1] += weight*A[1];
          result[2] += weight*A[2];
          weightsum += weight;
         }

       // 8 nodes
      if (sph_CN_boundary_flag == 8)
         {
          gridreal weight = 8.0*c->sph_size[0];
          //gridreal weight = 8.0;

          result[0] += weight*A[0];
          result[1] += weight*A[1];
          result[2] += weight*A[2];
          weightsum += weight;
         }
     }

 gridreal invweightsum = 1.0/weightsum;
 for (d=0; d<3; d++) temp[d] = result[d]*invweightsum;
 sph_transf_C2S_A1(rn, temp);
 for (d=0; d<3; d++) nodedata[ns][d] = temp[d];
}
*/

//! (SPHERICAL) Spherical version of "CN1". For Ji only
void Tgrid::Tnode::sph_CN1_Ji(TCellDataSelect cs, TNodeDataSelect ns)
{
    int a;
    gridreal result[3] = {0,0,0};
    gridreal rn[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal r[3]  = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
        gridreal Ji[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        for (d=0; d<3; d++) result[d]+= Ji[d];
    }
    for (d=0; d<3; d++) nodedata[ns][d] = result[d]/8;
}

//! (SPHERICAL) Spherical version of "CN1_Ji" + CN boundary calculations.
void Tgrid::Tnode::sph_CNb1_Ji(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag)
{
    int a;
    gridreal result[3] = {0, 0, 0};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        // 1 internal region
        if (sph_CN_boundary_flag == 1)  {
            for (d=0; d<3; d++) result[d]+= A[d];
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6)  {
            for (d=0; d<3; d++) result[d]+= 2*A[d];
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            for (d=0; d<3; d++) result[d]+= 4*A[d];
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8)  {
            for (d=0; d<3; d++) result[d]+= 8*A[d];
        }
    }
    for (d=0; d<3; d++) nodedata[ns][d] = result[d]/8;
}


//! (SPHERICAL) Spherical version of "CN1_smoothing"
void Tgrid::Tnode::sph_CN1_smoothing()
{
    int a;
    gridreal result[3] = {0, 0, 0};
    gridreal numberdensity = 0.0, chargedensity = 0.0;
    gridreal rn[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal A[3] = {c->celldata[CELLDATA_Ji][0], c->celldata[CELLDATA_Ji][1], c->celldata[CELLDATA_Ji][2]};
        // If Ji is calculated in Cartesian coordinates we don't neet to make any transformations
        numberdensity += c->nc;
        chargedensity += c->rho_q;
        for (d=0; d<3; d++) result[d]+= A[d];
    }
    nn = numberdensity/8;
    nodedata[NODEDATA_UE][0] = chargedensity/8;
    for (d=0; d<3; d++) nodedata[NODEDATA_J][d] = result[d]/8; //NODEDATA_J is used as a temporary storage place for ion current density on node
}

//! (SPHERICAL) Spherical version of "CN1_smoothing" + CN boundary calculations.
void Tgrid::Tnode::sph_CNb1_smoothing(int sph_CN_boundary_flag)
{
// VERY IMPORTANT
// The scheme of cells neighboured by node
// cell[0][0][0] = cell[0][0][0] (i  ,j  ,k  )
// cell[0][0][1] = cell[0][0][1] (i,  j,  k+1)
// cell[0][1][0] = cell[0][0][2] (i,  j+1,k  )
// cell[0][1][1] = cell[0][0][3] (i,  j+1,k+1)
// cell[1][0][0] = cell[0][0][4] (i+1,j,  k  )
// cell[1][0][1] = cell[0][0][5] (i+1,j,  k+1)
// cell[1][1][0] = cell[0][0][6] (i+1,j+1,k  )
// cell[1][1][1] = cell[0][0][7] (i+1,j+1,k+1)
    int a;
    gridreal result[3] = {0, 0, 0};
    gridreal numberdensity = 0.0, chargedensity = 0.0;
    gridreal rn[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        gridreal A[3] = {c->celldata[CELLDATA_Ji][0], c->celldata[CELLDATA_Ji][1], c->celldata[CELLDATA_Ji][2]};
        // If Ji is calculated in Cartesian coordinates we don't neet to make any transformations
        // 1 internal region
        if (sph_CN_boundary_flag == 1)  {
            numberdensity += 1*c->nc;
            chargedensity += 1*c->rho_q;
            for (d=0; d<3; d++) result[d]+= 1*A[d];
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6)  {
            numberdensity += 2*c->nc;
            chargedensity += 2*c->rho_q;
            for (d=0; d<3; d++) result[d]+= 2*A[d];
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            numberdensity += 4*c->nc;
            chargedensity += 4*c->rho_q;
            for (d=0; d<3; d++) result[d]+= 4*A[d];
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8)  {
            numberdensity += 8*c->nc;
            chargedensity += 8*c->rho_q;
            for (d=0; d<3; d++) result[d]+= 8*A[d];
        }
    }
    nn = numberdensity/8;
    nodedata[NODEDATA_UE][0] = chargedensity/8;
    for (d=0; d<3; d++) nodedata[NODEDATA_J][d] = result[d]/8; //NODEDATA_J is used as a temporary storage place for ion current density on node
}

//! (SPHERICAL) Spherical version of "CN1_rhoq" (volume method)
void Tgrid::Tnode::sph_CN1_rhoq()
{
    int a;
    gridreal weightsum = 0.0;
    gridreal invweightsum;
    register Tcell *c;
    real chargedensity=0.0;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        register const gridreal weight = 1.0/c->sph_dV;
        chargedensity+=weight*c->rho_q; //charge density
        weightsum+= weight;
    }
    if(weightsum==0.)invweightsum=0;
    else invweightsum=1./weightsum;
    nodedata[NODEDATA_UE][0]=chargedensity*invweightsum; //NODEDATA_UE is used as a temporary storage place for charge density on node
}

//! (SPHERICAL) Spherical version of "CN1_rhoq" (simple method)
void Tgrid::Tnode::sph_CNb1_rhoq(int sph_CN_boundary_flag)
{
    int a;
    register Tcell *c;
    real chargedensity=0.0;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        // 1 internal region
        if (sph_CN_boundary_flag == 1)  {
            chargedensity += 1*c->rho_q;
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6)  {
            chargedensity += 2*c->rho_q;
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            chargedensity += 4*c->rho_q;
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8)  {
            chargedensity += 8*c->rho_q;
        }
    }
    nodedata[NODEDATA_UE][0]=chargedensity/8; //NODEDATA_UE is used as a temporary storage place for charge density on node
}

//! (SPHERICAL) Spherical version of "CN_donor1". Uses uns to approximate ns (using CN-interpolation) in upstream (dx/2) of the uns field
void Tgrid::Tnode::sph_CN_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    int d,ch;
    real result[3] = {0,0,0};
    register real weightsum = 0.0;
    gridreal displacement[3] = {nodedata[uns][0], nodedata[uns][1], nodedata[uns][2]};
    Tcell *c;
    // Find the smallest cell size among the neighbours
    gridreal basic_dx = -1;
    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!ch) continue;
        if (c->size < basic_dx || basic_dx < 0) basic_dx = c->size;
    }
    basic_dx*= 0.5;
    // let basic_dx be one half of the smallest touching cell size
    real norm = sqrt(sqr(displacement[0]) + sqr(displacement[1]) + sqr(displacement[2]));
    if (norm == 0) {                          // nodedata velocity can be exactly zero (e.g., in the initial state), thus must protect ourselves
        displacement[0] = basic_dx;             // choose displacement pointing in x-direction (arbitrarily)
        displacement[1] = displacement[2] = 0;
    } else {                                  // the usual case: displacement is directed along velocity, but length is basic_dx
        norm = basic_dx/norm;
        for (d=0; d<3; d++) displacement[d]*= norm;
    }
    const gridreal upwind_x = centroid[0] - displacement[0];
    const gridreal upwind_y = centroid[1] - displacement[1];
    const gridreal upwind_z = centroid[2] - displacement[2];
    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!c) {
            errorlog << "Tnode::CN_donor1: nonexistent cell at " << Tr3v(centroid).toString() << ", upwind_x=" << upwind_x/3e6 << "\n";
            continue;
        }
        // notice: this is not strictly correct in a T-junction node
        const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
        //      const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
        const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
        result[0]+= weight*c->celldata[cs][0];
        result[1]+= weight*c->celldata[cs][1];
        result[2]+= weight*c->celldata[cs][2];
        weightsum+= weight;
    }
    const real invweightsum = 1.0/weightsum;
    for (d=0; d<3; d++) nodedata[ns][d] = invweightsum*result[d];
}

//! (SPHERICAL) Spherical version of "CN_donor1". Uses uns to approximate ns (using CN-interpolation) in upstream (dx/2) of the uns field.
void Tgrid::Tnode::sph_CNb_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt, int sph_CN_boundary_flag)
{
    int d,ch;
    real result[3] = {0,0,0};
    register real weightsum = 0.0;
    gridreal displacement[3] = {nodedata[uns][0], nodedata[uns][1], nodedata[uns][2]};
    gridreal temp[3];
    gridreal rn[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    Tcell *c;

// Find the smallest cell size among the neighbours
    gridreal basic_dx = -1;

    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!ch) continue;
        if (c->size < basic_dx || basic_dx < 0) basic_dx = c->size;
    }

    basic_dx*= 0.5;
// let basic_dx be one half of the smallest touching cell size
    real norm = sqrt(sqr(displacement[0]) + sqr(displacement[1]) + sqr(displacement[2]));
    if (norm == 0) {
        // nodedata velocity can be exactly zero (e.g., in the initial state), thus must protect ourselves
        displacement[0] = basic_dx;             // choose displacement pointing in x-direction (arbitrarily)
        displacement[1] = displacement[2] = 0;
    } else {
        // the usual case: displacement is directed along velocity, but length is basic_dx
        norm = basic_dx/norm;
        for (d=0; d<3; d++) displacement[d]*= norm;
    }
    const gridreal upwind_x = centroid[0] - displacement[0];
    const gridreal upwind_y = centroid[1] - displacement[1];
    const gridreal upwind_z = centroid[2] - displacement[2];
    //const gridreal upwind_x = 0.0;
    //const gridreal upwind_y = 0.0;
    //const gridreal upwind_z = 0.0;
    for (ch=0; ch<8; ch++) {
        c = cell[0][0][ch];
        if (!c) {
            errorlog << "Tnode::sph_CNb_donor1: nonexistent cell at " << Tr3v(centroid).toString() << ", upwind_x=" << upwind_x/3e6 << "\n";
            continue;
        }
        gridreal r[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
        gridreal A[3] = {c->celldata[cs][0], c->celldata[cs][1], c->celldata[cs][2]};
        sph_transf_S2C_A1(r, A);
        // 1 internal region
        if (sph_CN_boundary_flag == 1) {
            // notice: this is not strictly correct in a T-junction node
            const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
            //const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
            const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
            result[0]+= weight*A[0];
            result[1]+= weight*A[1];
            result[2]+= weight*A[2];
            weightsum+= weight;
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6) {
            // notice: this is not strictly correct in a T-junction node
            const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
            //const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
            const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
            result[0]+= 2*weight*A[0];
            result[1]+= 2*weight*A[1];
            result[2]+= 2*weight*A[2];
            weightsum+= weight;
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            // notice: this is not strictly correct in a T-junction node
            const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
            //const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
            const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
            result[0]+= 4*weight*A[0];
            result[1]+= 4*weight*A[1];
            result[2]+= 4*weight*A[2];
            weightsum+= weight;
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8) {
            // notice: this is not strictly correct in a T-junction node
            const gridreal dist2 = sqr(c->centroid[0]-upwind_x) + sqr(c->centroid[1]-upwind_y) + sqr(c->centroid[2]-upwind_z);
            //const gridreal weight = 1.0/(1e-30+dist2);          // faster: weight proportinal to inverse SQUARED DISTANCE
            const gridreal weight = 1.0/(1e-30+sqrt(dist2));    // slower but maybe better: weight proportional to inverse DISTANCE
            result[0]+= 8*weight*A[0];
            result[1]+= 8*weight*A[1];
            result[2]+= 8*weight*A[2];
            weightsum+= weight;
        }
    }
    const real invweightsum = 1.0/weightsum;
    for (d=0; d<3; d++) temp[d] = result[d];
    sph_transf_C2S_A1(rn, temp);
    for (d=0; d<3; d++) nodedata[ns][d] = invweightsum*temp[d];
}

/** \brief Spherical version of "calc_E1": Calculate E = -U_e x (B + B_0) + eta*J
 *
 * Here is used cartesian method of calculation E = - Ue x B.
 * 1. Vectors Ue and B transformed to catresian coordinates.
 * 2. Calculation of E = - Ue x B.
 * 3. E transformation to spherical coordinates.
 */
void Tgrid::Tnode::sph_calc_E1(void)
{
    real B0[3] = {0.0, 0.0, 0.0};
    real B[3]  = {0.0, 0.0, 0.0};
    real Ue[3] = {0.0, 0.0, 0.0};
    real E[3]  = {0.0, 0.0, 0.0};
    // Add constant magnetic field B0 (centroid = exact node coordinates)
    addConstantMagneticField(centroid, B0);
    gridreal R     = sph_centroid[0];
    gridreal theta = sph_centroid[1];
    gridreal phi   = sph_centroid[2];
    gridreal r[3]  = {R, theta, phi};
    B[0] = nodedata[NODEDATA_B][0] + B0[0];
    B[1] = nodedata[NODEDATA_B][1] + B0[1];
    B[2] = nodedata[NODEDATA_B][2] + B0[2];
    Ue[0] = nodedata[NODEDATA_UE][0];
    Ue[1] = nodedata[NODEDATA_UE][1];
    Ue[2] = nodedata[NODEDATA_UE][2];
    // 1. Transformation of Ue and B from spherical to cartesian coordinates
    sph_transf_S2C_A(r, B);
    sph_transf_S2C_A(r, Ue);
    // 2. Calculation of E = - Ue x B.
    E[0] = B[1]*Ue[2] - B[2]*Ue[1];
    E[1] = B[2]*Ue[0] - B[0]*Ue[2];
    E[2] = B[0]*Ue[1] - B[1]*Ue[0];
    // 4. Add resistivity.
    real j[3] = {nodedata[NODEDATA_J][0], nodedata[NODEDATA_J][1], nodedata[NODEDATA_J][2]};
    sph_transf_S2C_A(r, j);
    E[0] += eta*j[0];
    E[1] += eta*j[1];
    E[2] += eta*j[2];
    // 3. Transformation of E to spherical coordinates.
    sph_transf_C2S_A(r, E);
    nodedata[NODEDATA_E][0] = E[0];
    nodedata[NODEDATA_E][1] = E[1];
    nodedata[NODEDATA_E][2] = E[2];
    // Check electric field cut value (CONSTRAINT)
    if (Params::Ecut > 0) {
        sph_transf_S2C_A(r, E);
        real absE = sqrt(sqr(E[0]) + sqr(E[1]) + sqr(E[2]));
        if (absE > Params::Ecut) {
            real scaling = Params::Ecut/absE;
            E[0] *= scaling;
            E[1] *= scaling;
            E[2] *= scaling;
            sph_transf_C2S_A(r, E);
            nodedata[NODEDATA_E][0] = E[0];
            nodedata[NODEDATA_E][1] = E[1];
            nodedata[NODEDATA_E][2] = E[2];
#ifndef NO_DIAGNOSTICS
            // Increase counter
            Tgrid::fieldCounter.cutRateE += 1.0;
#endif
        }
    }
}

//! (SPHERICAL) Spherical node boundary conditions
void Tgrid::Tnode::sph_Node_BC1(TNodeDataSelect ns)
{
    real t = Params::t;
    real V = abs(Params::V);
    real x_min = Params::box_xmin;
    real x_max = Params::box_xmax;
    real S = V*t;
    real A[3] = {0.0, 0.0, 0.0};
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    gridreal z    = r[0]*cos(r[1]);
    gridreal x    = r[0]*sin(r[1])*cos(r[2]);
// Radial propagation
    if (Params::sph_propagation_type == 0) {
        A[0] = Params::boundary_Bx;
        A[1] = Params::boundary_By;
        A[2] = Params::boundary_Bz*sin(r[1])*sin(r[1]);
        nodedata[ns][0] = A[0];
        nodedata[ns][1] = A[1];
        nodedata[ns][2] = A[2];
    }
// Flat front propagation
    if (Params::sph_propagation_type == 1) {
        //Propagation along x-axis
        if (Params::sph_propagation_dir == 0) {
            if (x >= 0.0 && x >= x_max - S) {
                A[0] = Params::boundary_Bx;
                A[1] = Params::boundary_By;
                A[2] = Params::boundary_Bz;
                sph_transf_C2S_A(r, A);
                nodedata[ns][0] = A[0];
                nodedata[ns][1] = A[1];
                nodedata[ns][2] = A[2];
            }
        }
        //Propagation along z-axis
        if (Params::sph_propagation_dir == 1) {
            if (z >= 0.0 && z >= x_max - S) {
                A[0] = Params::boundary_Bx;
                A[1] = Params::boundary_By;
                A[2] = Params::boundary_Bz;
                sph_transf_C2S_A(r, A);
                nodedata[ns][0] = A[0];
                nodedata[ns][1] = A[1];
                nodedata[ns][2] = A[2];
            }
        } else {
            return;
        }
    }
}

//! (SPHERICAL) Spherical version of "CN_recursive"
void Tgrid::Tcell::sph_CN_recursive(TCellDataSelect cs, TNodeDataSelect ns)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_CN_recursive(cs,ns);
    } else {
        // NODE LOOP
        // Idea: To pass through all nodes, pass through all right-pointing faces,
        // and all upper-right corner points thereof
        // IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1(cs,ns);
        } else {
            face[d][1]->node[2]->sph_CN1(cs,ns); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1(cs,ns);
            }
    }
}

//! (SPHERICAL) Spherical version of "CN_recursive"
void Tgrid::Tcell::sph_CNb_recursive(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag)
{
// NODE LOOP
// Idea: To pass through all nodes, pass through all right-pointing faces,
// and all upper-right corner points thereof
// IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
    int d = 0;          // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb1(cs,ns,sph_CN_boundary_flag);
}

//! (SPHERICAL) Spherical version of "CN_recursive" based on cell volume
void Tgrid::Tcell::sph_CN_vol_recursive(TCellDataSelect cs, TNodeDataSelect ns)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_CN_vol_recursive(cs,ns);
    } else {
        // NODE LOOP
        // Idea: To pass through all nodes, pass through all right-pointing faces,
        // and all upper-right corner points thereof
        // IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_vol(cs,ns);
        } else {
            face[d][1]->node[2]->sph_CN1_vol(cs,ns); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_vol(cs,ns);
            }
    }
}

//! (SPHERICAL) Spherical version of "CN_recursive". For Ji only
void Tgrid::Tcell::sph_CN_Ji_recursive(TCellDataSelect cs, TNodeDataSelect ns)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_CN_Ji_recursive(cs,ns);
    } else {
        // NODE LOOP
        // Idea: To pass through all nodes, pass through all right-pointing faces,
        // and all upper-right corner points thereof
        // IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_Ji(cs,ns);
        } else {
            face[d][1]->node[2]->sph_CN1_Ji(cs,ns); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_Ji(cs,ns);
            }
    }
}

//! (SPHERICAL) Spherical version of "CN_recursive" + CN boundary condition. For Ji only
void Tgrid::Tcell::sph_CNb_Ji_recursive(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag)
{
// NODE LOOP
// Idea: To pass through all nodes, pass through all right-pointing faces,
// and all upper-right corner points thereof
// IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
    int d = 0;          // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb1_Ji(cs,ns,sph_CN_boundary_flag);
}

//! (SPHERICAL) Spherical version of "CN_smoothing_recursive"
void Tgrid::Tcell::sph_CN_smoothing_recursive()
{
// NODE LOOP
    int d = 0;                                // choose to pass in X-direction
    face[d][1]->node[2]->sph_CN1_smoothing(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
}

//! (SPHERICAL) Spherical version of "CN_smoothing_recursive"
void Tgrid::Tcell::sph_CNb_smoothing_recursive(int sph_CN_boundary_flag)
{
// NODE LOOP
    int d = 0;                                                 // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb1_smoothing(sph_CN_boundary_flag); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
}

//! (SPHERICAL) Spherical version of "CN_rhoq_recursive"
void Tgrid::Tcell::sph_CN_rhoq_recursive()
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_CN_rhoq_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_rhoq();
        } else {
            face[d][1]->node[2]->sph_CN1_rhoq(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_rhoq();
            }
    }
}

//! (SPHERICAL) Spherical version of "CN_rhoq_recursive"
void Tgrid::Tcell::sph_CNb_rhoq_recursive(int sph_CN_boundary_flag)
{
// NODE LOOP
    int d = 0;                                                 // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb1_rhoq(sph_CN_boundary_flag); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
}

//! (SPHERICAL) Spherical version of "CN_rhoq_recursive"
void Tgrid::Tcell::sph_CN_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
// NODE LOOP
    int d = 0;          // choose to pass in X-direction
    face[d][1]->node[2]->sph_CN_donor1(cs,ns,uns,dt); //assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
}

//! (SPHERICAL) Spherical version of "CN_rhoq_recursive" + Boundary calculations without ghost cells
void Tgrid::Tcell::sph_CNb_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt, int sph_CN_boundary_flag)
{
// NODE LOOP
    int d = 0;          // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb_donor1(cs,ns,uns,dt,sph_CN_boundary_flag); //assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
}

//! (SPHERICAL) Spherical version of "calc_node_E_recursive". See CN_recursive above. This is almost the same, just the call calc_E1() is different.
void Tgrid::Tcell::sph_calc_node_E_recursive(void)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_node_E_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_E1();
        } else {
            // without refinement we are interesting only this line
            face[d][1]->node[2]->sph_calc_E1();
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_E1();
            }
    }
}

//! (SPHERICAL) Spherical node boundary conditions
void Tgrid::Tcell::sph_Node_BC_recursive(TNodeDataSelect ns)
{
    int d = 0; // choose to pass in X-direction

    face[d][1]->node[2]->sph_Node_BC1(ns);
}

//! (SPHERICAL) Spherical version of "calc_cell_E_recursive": calculate E field inside the cell and save the result temporarily in CELLDATA_TEMP2.
void Tgrid::Tcell::sph_calc_cell_E_recursive(void)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_cell_E_recursive();
    } else {
        celldata[CELLDATA_TEMP2][0] = 0.;
        celldata[CELLDATA_TEMP2][1] = 0.;
        celldata[CELLDATA_TEMP2][2] = 0.;
        //electron pressure contribution
        if (sph_centroid[0] > Params::R_zeroPolarizationField) {
            //if the cell is outside of the background ionosphere density peak, include the polarization electric field.
            // E = - kT/e*nabla(rho_q)/rho_q
            const real pressure_coef=Params::k_B*Params::Te/Params::e; //KTe/e for isothermal plasma.
            celldata[CELLDATA_TEMP2][0] -= pressure_coef*celldata[CELLDATA_TEMP1][0]/rho_q;
            celldata[CELLDATA_TEMP2][1] -= pressure_coef*celldata[CELLDATA_TEMP1][1]/rho_q;
            celldata[CELLDATA_TEMP2][2] -= pressure_coef*celldata[CELLDATA_TEMP1][2]/rho_q;
        }
    }
}

//! (SPHERICAL) Spherical version of "calc_ue_recursive"
void Tgrid::Tcell::sph_calc_ue_recursive(void)
{
    if (haschildren) {
        for (int ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_ue_recursive();
    } else {
        // We add one r layer to be sure that at least one sph_centroin inside Ue = 0 region
        // Check zero field radius. It must be used only for flat front propagation
        if (Params::sph_propagation_type == 1) {
            if (sph_centroid[0] < Params::R_zeroFields) {
                for (int d=0; d<3; d++) {
                    celldata[CELLDATA_UE][d] = 0;
                }
                return;
            }
        }
        // Charge density in thel cell
        real chargedensity = rho_q;
        const real invrho_q = 1.0/chargedensity;
        gridreal ue[3];
        real ue2 = 0.0;
        // Ji comes in Cartesian coordinates and J in spherical, so we need to transform J to Cartesian
        gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
        gridreal J[3] = {celldata[CELLDATA_J][0], celldata[CELLDATA_J][1], celldata[CELLDATA_J][2]};
        sph_transf_S2C_A1(r, J);
        for (int d=0; d<3; d++) {
#ifndef IGNORE_ELECTRIC_FIELD_HALL_TERM
            ue[d] = (celldata[CELLDATA_Ji][d] - J[d])*invrho_q;
#else
            ue[d] = celldata[CELLDATA_Ji][d]*invrho_q;
#endif
            ue2+= sqr(ue[d]);
        }
        // Check maximum electron fluid velocity (CONSTRAINT)
        if (Params::Ue_max > 0 && ue2 > Params::Ue_max2) {
            const real norm = Params::Ue_max/sqrt(ue2);
            for (int d=0; d<3; d++) ue[d]*= norm;
            // Increase counter
#ifndef NO_DIAGNOSTICS
            Tgrid::fieldCounter.cutRateUe += 1.0;
#endif
        }
        //gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
        //sph_transf_H2S_V(ue);
        //sph_transf_S2C_V(r, ue);
        sph_transf_C2S_A1(r, ue);
        for (int d=0; d<3; d++) celldata[CELLDATA_UE][d] = ue[d];
    }
}

//! (SPHERICAL) Spherical version of "FaceCurl_recursive"
void Tgrid::Tcell::sph_FaceCurl_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_FaceCurl_recursive(ns,fs,d,factor);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->sph_Curl1(ns,fs,d,factor);
        } else {
            face[d][dir]->sph_Curl1(ns,fs,d,factor); // without refinement we need only this line
        }
    }
}

//! (SPHERICAL) New Spherical version of "FaceCurl_recursive"
void Tgrid::Tcell::sph_FaceCurl_2_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_FaceCurl_2_recursive(ns,fs,d,factor);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->sph_Curl1_2(ns,fs,d,factor);
        } else {
            face[d][dir]->sph_Curl1_2(ns,fs,d,factor);
        }
    }
}

//! (SPHERICAL) Spherical version of "NF_recursive"
void Tgrid::Tcell::sph_NF_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_NF_recursive(ns,fs,d);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->sph_NF1(ns,fs,d);
        } else {
            face[d][dir]->sph_NF1(ns,fs,d);
        }
    }
}

//! (SPHERICAL) Spherical version of "NF_rhoq_recursive"
void Tgrid::Tcell::sph_NF_rhoq_recursive(int d)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_NF_rhoq_recursive(d);
    } else {
        const int dir = 1;
        if (isrefined_face(d,dir)) {
            int f;
            for (f=0; f<4; f++) refintf[d][dir]->face[f]->sph_NF_rhoq1(d);
        } else {
            face[d][dir]->sph_NF_rhoq1(d);
        }
    }
}

//! (SPHERICAL) Spherical version of "calc_facediv": Compute magnetic log values
void Tgrid::sph_calc_facediv(TFaceDataSelect fs, TCellDataSelect cs, MagneticLog& result) const
{
    MagneticLog BLog;
    // Loop thru the grid cells
    int i,j,k;
    ForInterior(i,j,k) {
        int c = flatindex(i,j,k);
        cells[c]->sph_calc_facediv_recursive(fs,cs,BLog);
        // Comparate current max values to new ones
        if(BLog.maxDivB > result.maxDivB) {
            result.maxDivB = BLog.maxDivB;
        }
        if(BLog.maxB > result.maxB) {
            result.maxB = BLog.maxB;
            result.posMaxB[0] = BLog.posMaxB[0];
            result.posMaxB[1] = BLog.posMaxB[1];
            result.posMaxB[2] = BLog.posMaxB[2];
        }
        if(BLog.maxDxDivBperB > result.maxDxDivBperB) {
            result.maxDxDivBperB = BLog.maxDxDivBperB;
        }
        // Sum average values and energy
        result.avgDivB += BLog.avgDivB;
        result.avgBx += BLog.avgBx;
        result.avgBy += BLog.avgBy;
        result.avgBz += BLog.avgBz;
        result.avgB += BLog.avgB;
        result.energyB += BLog.energyB;
    }
    // Normalize average values
    real normCoeff = Params::nx * Params::ny * Params::nz;
    result.avgBx /= normCoeff;
    result.avgBy /= normCoeff;
    result.avgBz /= normCoeff;
    result.avgB /= normCoeff;
    result.avgDivB /= normCoeff;
    // Energy dimension
    result.energyB *= cube(bgdx)/(2 * Params::mu_0);
}

//! (SPHERICAL) Spherical version of "calc_facediv_recursive": Magnetic log calculation
void Tgrid::Tcell::sph_calc_facediv_recursive(TFaceDataSelect fs, TCellDataSelect cs, MagneticLog& result) const
{
    if (haschildren) {
        MagneticLog BLog;
        // Loop thru the child cells
        for (int ch=0; ch<8; ch++) {
            child[0][0][ch]->sph_calc_facediv_recursive(fs,cs,BLog);
            // Comparate current max values to new ones
            if(BLog.maxDivB > result.maxDivB) {
                result.maxDivB = BLog.maxDivB;
            }
            if(BLog.maxB > result.maxB) {
                result.maxB = BLog.maxB;
                result.posMaxB[0] = BLog.posMaxB[0];
                result.posMaxB[1] = BLog.posMaxB[1];
                result.posMaxB[2] = BLog.posMaxB[2];
            }
            if(BLog.maxDxDivBperB > result.maxDxDivBperB) {
                result.maxDxDivBperB = BLog.maxDxDivBperB;
            }
            // Sum average values and energy
            result.avgBx += BLog.avgBx;
            result.avgBy += BLog.avgBy;
            result.avgBz += BLog.avgBz;
            result.avgB += BLog.avgB;
            result.avgDivB += BLog.avgDivB;
            result.energyB += BLog.energyB;
        }
        // Child cell normalization
        result.avgDivB *= 0.125;
        result.avgBx *= 0.125;
        result.avgBy *= 0.125;
        result.avgBz *= 0.125;
        result.avgB *= 0.125;
        result.energyB *= 0.125;
    } else {
        // Flux out of cell, aveBxyz and posMaxB
        real flux = 0.0;
        real tempAvgB[3];
        real A[3];
        gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
        for (int d=0; d<3; d++) {
            flux += faceave(d,1,fs) - faceave(d,0,fs);
            A[d] = celldata[cs][d];
            //tempAvgB[d] = 0.5 *(faceave(d,1,fs) + faceave(d,0,fs));
            result.posMaxB[d] = centroid[d];
        }
        sph_transf_S2C_A(r, A);
        tempAvgB[0] = A[0];
        tempAvgB[1] = A[1];
        tempAvgB[2] = A[2];
        result.avgBx = tempAvgB[0];
        result.avgBy = tempAvgB[1];
        result.avgBz = tempAvgB[2];
        // B^2 in this cell
        result.energyB = sqr(result.avgBx) + sqr(result.avgBy) + sqr(result.avgBz);
        // B in this cell
        result.maxB = result.avgB = sqrt(result.energyB);
        // divB in this cell
        result.avgDivB = result.maxDivB = fabs(flux)*invsize;
        // DxDivBperB in this cell
        if (result.avgB > 1e-20) {
            result.maxDxDivBperB = size*result.avgDivB/result.avgB;
        } else {
            result.maxDxDivBperB = 0.0;
        }
    }
}

//! (SPHERICAL) Spherical version of "CalcGradient_rhoq_recursive"
void Tgrid::Tcell::sph_CalcGradient_rhoq_recursive(void)
{
    if (haschildren) {
        for (int ch=0; ch<8; ch++) child[0][0][ch]->sph_CalcGradient_rhoq_recursive();
    } else {
        // dLr     = dr
        // dLtheta = r_c*dtheta
        // dLphi   =   r_c*sin(theta_c)*dphi
        gridreal dL[3];
        dL[0] = sph_size[0];
        dL[1] = sph_centroid[0]*sph_size[1];
        dL[2] = sph_centroid[0]*sin(sph_centroid[1])*sph_size[2];
        for (int d=0; d<3; d++) {
            celldata[CELLDATA_TEMP1][d]= (faceave(d,1,Tgrid::FACEDATA_MINUSDB) - faceave(d,0,Tgrid::FACEDATA_MINUSDB))/dL[d];
        }
    }
}

//! (SPHERICAL) Spherical version of "Curl1": This is modified version which deal with the counterclockwise vector circulation
void Tgrid::Tface::sph_Curl1(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor)
{
    // interpolation from nodes to sidenodes
    /* int dd = 3;
    gridreal A0[3] = {0, 0, 0};
    gridreal A1[3] = {0, 0, 0};
    gridreal A2[3] = {0, 0, 0};
    gridreal A3[3] = {0, 0, 0};*/
    gridreal node0_F[3] = {node[0]->nodedata[ns][0], node[0]->nodedata[ns][1], node[0]->nodedata[ns][2]};
    gridreal node1_F[3] = {node[1]->nodedata[ns][0], node[1]->nodedata[ns][1], node[1]->nodedata[ns][2]};
    gridreal node2_F[3] = {node[2]->nodedata[ns][0], node[2]->nodedata[ns][1], node[2]->nodedata[ns][2]};
    gridreal node3_F[3] = {node[3]->nodedata[ns][0], node[3]->nodedata[ns][1], node[3]->nodedata[ns][2]};
    gridreal r_node0[3] = {node[0]->sph_centroid[0], node[0]->sph_centroid[1], node[0]->sph_centroid[2]};
    gridreal r_node1[3] = {node[1]->sph_centroid[0], node[1]->sph_centroid[1], node[1]->sph_centroid[2]};
    gridreal r_node2[3] = {node[2]->sph_centroid[0], node[2]->sph_centroid[1], node[2]->sph_centroid[2]};
    gridreal r_node3[3] = {node[3]->sph_centroid[0], node[3]->sph_centroid[1], node[3]->sph_centroid[2]};
    //  Note: the number 0-3 inside [] is the index used in hybrid model
    //  [3] --- line23 ----  [2] ------ phi2      length of line01 =  r*dtheta
    //   |                    |                   length of line23 = -r*dtheta (dtheta < 0)
    //   |                    |                   length of line12 =  r*sin(theta2)*dphi
    // line30    Br(o)      line12       \/       length of line30 = -r*sin(theta1)*dphi (dphi < 0)
    //   |                    |
    //   |                    |
    //  [0] --- line01 ----  [1] ------ phi1
    //   |                    |
    //   |                    |
    // theta1      <        theta2
    real dS = 0.0;
    if (d==0) {
        const gridreal r       =  r_node0[0];
        const gridreal theta1  =  r_node0[1];
        const gridreal theta2  =  r_node1[1];
        const gridreal phi1    =  r_node0[2];
        const gridreal phi2    =  r_node2[2];
        const gridreal dtheta  =  theta2 - theta1;
        const gridreal dphi    =  phi2 - phi1;
        // interpolation from nodes to sidenodes
        /*gridreal r01c[3] = {r, theta1 + dtheta/2, phi1};
        gridreal r12c[3] = {r, theta2,            phi1 + dphi/2};
        gridreal r23c[3] = {r, theta1 + dtheta/2, phi2};
        gridreal r30c[3] = {r, theta1,            phi1 + dphi/2};
        //------------------- line 01 -------------------------
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        sph_transf_S2C_A1(r01c, A0);
        sph_transf_S2C_A1(r01c, A1);
        sph_transf_C2S_A1(r01c, A0);
        sph_transf_C2S_A1(r01c, A1);
        const real L01 =  0.5*(A0[1] + A1[1])*r*dtheta;
        //------------------- line 12 -------------------------
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        sph_transf_S2C_A1(r12c, A1);
        sph_transf_S2C_A1(r12c, A2);
        sph_transf_C2S_A1(r12c, A1);
        sph_transf_C2S_A1(r12c, A2);
        const real L12 =  0.5*(A1[2] + A2[2])*r*sin(theta2)*dphi;
        //------------------- line 23 -------------------------
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        sph_transf_S2C_A1(r23c, A2);
        sph_transf_S2C_A1(r23c, A3);
        sph_transf_C2S_A1(r23c, A2);
        sph_transf_C2S_A1(r23c, A3);
        const real L23 = -0.5*(A2[1] + A3[1])*r*dtheta;
        //------------------- line 30 -------------------------
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        sph_transf_S2C_A1(r30c, A3);
        sph_transf_S2C_A1(r30c, A0);
        sph_transf_C2S_A1(r30c, A3);
        sph_transf_C2S_A1(r30c, A0);
        const real L30 = -0.5*(A3[2] + A0[2])*r*sin(theta1)*dphi;*/
        const real L01 =  0.5*(node0_F[1] + node1_F[1])*r*dtheta;
        const real L12 =  0.5*(node1_F[2] + node2_F[2])*r*sin(theta2)*dphi;
        const real L23 = -0.5*(node2_F[1] + node3_F[1])*r*dtheta;
        const real L30 = -0.5*(node3_F[2] + node0_F[2])*r*sin(theta1)*dphi;
        dS = r*r*dphi*(cos(theta1)-cos(theta2));
        if (dS<0.0) dS = -dS;
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
    // Note: the number 0-3 inside [] is the index used in hybrid model
    //    [2]________
    //     | line12  \________[1] ----- phi2     length of line01 =  r1*sin(theta)*dphi
    //     |                   |                 length of line23 = -r2*sin(theta)*dphi (dphi < 0)
    //     |                   |                 length of line12 =  dr
    //   line23   Btheta(o)  line01              length of line30 = -dr (dr < 0)
    //     |                   |
    //     |         _________ |
    //     | _______/ line30  [0] ----- phi1
    //    [3]
    //     |                   |
    //    r2                   r1
    if (d==1) {
        const gridreal r1     =  r_node0[0];
        const gridreal r2     =  r_node2[0];
        const gridreal theta  =  r_node0[1];
        const gridreal phi1   =  r_node0[2];
        const gridreal phi2   =  r_node1[2];
        const gridreal dphi   =  phi2 - phi1;
        const gridreal dr     =  r2 - r1;
        /*// interpolation from nodes to sidenodes
        gridreal r01c[3] = {r1,        theta, phi1 + dphi/2};
        gridreal r12c[3] = {r1 + dr/2, theta, phi2};
        gridreal r23c[3] = {r2,        theta, phi1 + dphi/2};
        gridreal r30c[3] = {r1 + dr/2, theta, phi1};
        //------------------- line 01 -------------------------
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        sph_transf_S2C_A1(r01c, A0);
        sph_transf_S2C_A1(r01c, A1);
        sph_transf_C2S_A1(r01c, A0);
        sph_transf_C2S_A1(r01c, A1);
        const real L01 =  0.5*(A0[2] + A1[2])*r1*sin(theta)*dphi;
        //------------------- line 12 -------------------------
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        sph_transf_S2C_A1(r12c, A1);
        sph_transf_S2C_A1(r12c, A2);
        sph_transf_C2S_A1(r12c, A1);
        sph_transf_C2S_A1(r12c, A2);
        const real L12 =  0.5*(A1[0] + A2[0])*dr;
        //------------------- line 23 -------------------------
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        sph_transf_S2C_A1(r23c, A2);
        sph_transf_S2C_A1(r23c, A3);
        sph_transf_C2S_A1(r23c, A2);
        sph_transf_C2S_A1(r23c, A3);
        const real L23 = -0.5*(A2[2] + A3[2])*r2*sin(theta)*dphi;
        //------------------- line 30 -------------------------
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        sph_transf_S2C_A1(r30c, A3);
        sph_transf_S2C_A1(r30c, A0);
        sph_transf_C2S_A1(r30c, A3);
        sph_transf_C2S_A1(r30c, A0);
        const real L30 = -0.5*(A3[0] + A0[0])*dr;*/
        const real L01 =  0.5*(node0_F[2] + node1_F[2])*r1*sin(theta)*dphi;
        const real L12 =  0.5*(node1_F[0] + node2_F[0])*dr;
        const real L23 = -0.5*(node2_F[2] + node3_F[2])*r2*sin(theta)*dphi;
        const real L30 = -0.5*(node3_F[0] + node0_F[0])*dr;
        dS = 0.5*(r2*r2 - r1*r1)*sin(theta)*dphi;
        if (dS<0.0) dS = -dS;
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
    // Note: the number 0-3 inside [] is the index used in hybrid model
    /*          [0]- line30 -[3]       ----- r1    length of line01 =  dr
                /             \                    length of line12 =  r2*dtheta
               /               \                   length of line23 = -dr (dr < 0)
              /                 \                  length of line30 = -r1*dtheta (dtheta < 0)
           line01    Bphi(o)   line23
            /                     \
           /                       \
          /                         \
        [1]-------- line12 --------[2] ----- r2
         |                          |
       theta1                     theta2
    */
    if (d==2) {
        const gridreal r1      =  r_node0[0];
        const gridreal r2      =  r_node1[0];
        const gridreal theta1  =  r_node0[1];
        const gridreal theta2  =  r_node2[1];
        const gridreal phi     =  r_node0[2];
        gridreal dr      =  r2 - r1;         // Note: d_r does NOT have to be a constant !
        gridreal dtheta  =  theta2 - theta1; // Note: d_theta does NOT have to be a constant !
        /*// interpolation from nodes to sidenodes
        gridreal r01c[3] = {r1 + dr/2, theta1,            phi};
        gridreal r12c[3] = {r2,        theta1 + dtheta/2, phi};
        gridreal r23c[3] = {r1 + dr/2, theta2,            phi};
        gridreal r30c[3] = {r1,        theta1 + dtheta/2, phi};
        //------------------- line 01 -------------------------
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        sph_transf_S2C_A1(r01c, A0);
        sph_transf_S2C_A1(r01c, A1);
        sph_transf_C2S_A1(r01c, A0);
        sph_transf_C2S_A1(r01c, A1);
        const real L01 =  0.5*(A0[0] + A1[0])*dr;
        //------------------- line 12 -------------------------
        for (d=0; d<3; d++) A1[d] =  node1_F[d];
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        sph_transf_S2C_A1(r12c, A1);
        sph_transf_S2C_A1(r12c, A2);
        sph_transf_C2S_A1(r12c, A1);
        sph_transf_C2S_A1(r12c, A2);
        const real L12 =  0.5*(A1[1] + A2[1])*r2*dtheta;
        //------------------- line 23 -------------------------
        for (d=0; d<3; d++) A2[d] =  node2_F[d];
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        sph_transf_S2C_A1(r23c, A2);
        sph_transf_S2C_A1(r23c, A3);
        sph_transf_C2S_A1(r23c, A2);
        sph_transf_C2S_A1(r23c, A3);
        const real L23 = -0.5*(A2[0] + A3[0])*dr;
        //------------------- line 30 -------------------------
        for (d=0; d<3; d++) A3[d] =  node3_F[d];
        for (d=0; d<3; d++) A0[d] =  node0_F[d];
        sph_transf_S2C_A1(r30c, A3);
        sph_transf_S2C_A1(r30c, A0);
        sph_transf_C2S_A1(r30c, A3);
        sph_transf_C2S_A1(r30c, A0);
        const real L30 = -0.5*(A3[1] + A0[1])*r1*dtheta;*/
        const real L01 =  0.5*(node0_F[0]+node1_F[0])*dr;
        const real L12 =  0.5*(node1_F[1]+node2_F[1])*r2*dtheta;
        const real L23 = -0.5*(node2_F[0]+node3_F[0])*dr;
        const real L30 = -0.5*(node3_F[1]+node0_F[1])*r1*dtheta;
        dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
        if (dS<0.0) dS = -dS; // Just in case ...
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
}

//! (SPHERICAL) New Spherical version of "Curl1". Cartesian approach.
void Tgrid::Tface::sph_Curl1_2(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor)
{
    gridreal r_node0[3] = {node[0]->sph_centroid[0], node[0]->sph_centroid[1], node[0]->sph_centroid[2]};
    gridreal r_node1[3] = {node[1]->sph_centroid[0], node[1]->sph_centroid[1], node[1]->sph_centroid[2]};
    gridreal r_node2[3] = {node[2]->sph_centroid[0], node[2]->sph_centroid[1], node[2]->sph_centroid[2]};
    gridreal r_node3[3] = {node[3]->sph_centroid[0], node[3]->sph_centroid[1], node[3]->sph_centroid[2]};
    gridreal A_node0[3] = {node[0]->nodedata[ns][0], node[0]->nodedata[ns][1], node[0]->nodedata[ns][2]};
    gridreal A_node1[3] = {node[1]->nodedata[ns][0], node[1]->nodedata[ns][1], node[1]->nodedata[ns][2]};
    gridreal A_node2[3] = {node[2]->nodedata[ns][0], node[2]->nodedata[ns][1], node[2]->nodedata[ns][2]};
    gridreal A_node3[3] = {node[3]->nodedata[ns][0], node[3]->nodedata[ns][1], node[3]->nodedata[ns][2]};
    sph_transf_S2C_A1(r_node0, A_node0);
    sph_transf_S2C_A1(r_node1, A_node1);
    sph_transf_S2C_A1(r_node2, A_node2);
    sph_transf_S2C_A1(r_node3, A_node3);
    real Ax0 = A_node0[0];
    real Ay0 = A_node0[1];
    real Az0 = A_node0[2];
    real Ax1 = A_node1[0];
    real Ay1 = A_node1[1];
    real Az1 = A_node1[2];
    real Ax2 = A_node2[0];
    real Ay2 = A_node2[1];
    real Az2 = A_node2[2];
    real Ax3 = A_node3[0];
    real Ay3 = A_node3[1];
    real Az3 = A_node3[2];
    real dS = 0.0;
    //  Note: the number 0-3 inside [] is the index used in hybrid model
    //  [3] --- line23 ----  [2] ------ phi2      length of line01 =  r*dtheta
    //   |                    |                   length of line23 = -r*dtheta (dtheta < 0)
    //   |                    |                   length of line12 =  r*sin(theta2)*dphi
    // line30    Br(o)      line12       \/       length of line30 = -r*sin(theta1)*dphi (dphi < 0)
    //   |                    |
    //   |                    |
    //  [0] --- line01 ----  [1] ------ phi1
    //   |                    |
    //   |                    |
    // theta1      <        theta2
    if (d==0) {
        const gridreal r       =  r_node0[0];
        const gridreal theta1  =  r_node0[1];
        const gridreal theta2  =  r_node1[1];
        const gridreal phi1    =  r_node0[2];
        const gridreal phi2    =  r_node2[2];
        gridreal dtheta  =  theta2 - theta1;
        gridreal dphi    =  phi2 - phi1;
        real L01 =  0.5*(Ax0 + Ax1)*r*cos(phi1)*(sin(theta2)-sin(theta1)) + 0.5*(Ay0 + Ay1)*r*sin(phi1)*(sin(theta2)-sin(theta1))
                    +0.5*(Az0 + Az1)*r*(cos(theta2)-cos(theta1));
        real L23 =  0.5*(Ax2 + Ax3)*r*cos(phi2)*(sin(theta1)-sin(theta2)) + 0.5*(Ay2 + Ay3)*r*sin(phi2)*(sin(theta1)-sin(theta2))
                    +0.5*(Az2 + Az3)*r*(cos(theta1)-cos(theta2));
        real L12 =  0.5*(Ax1 + Ax2)*r*sin(theta2)*(cos(phi2)-cos(phi1)) + 0.5*(Ay1 + Ay2)*r*sin(theta2)*(sin(phi2)-sin(phi1));
        real L30 =  0.5*(Ax3 + Ax0)*r*sin(theta1)*(cos(phi1)-cos(phi2)) + 0.5*(Ay3 + Ay0)*r*sin(theta1)*(sin(phi1)-sin(phi2));
        dS = r*r*dphi*(cos(theta1)-cos(theta2));
        if (dS<0.0) dS = -dS;
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
    // Note: the number 0-3 inside [] is the index used in hybrid model
    //    [2]________
    //     | line12  \________[1] ----- phi2     length of line01 =  r1*sin(theta)*dphi
    //     |                   |                 length of line23 = -r2*sin(theta)*dphi (dphi < 0)
    //     |                   |                 length of line12 =  dr
    //   line23   Btheta(o)  line01              length of line30 = -dr (dr < 0)
    //     |                   |
    //     |         _________ |
    //     | _______/ line30  [0] ----- phi1
    //    [3]
    //     |                   |
    //    r2                   r1
    if (d==1) {
        const gridreal r1     =  r_node0[0];
        const gridreal r2     =  r_node2[0];
        const gridreal theta  =  r_node0[1];
        const gridreal phi1   =  r_node0[2];
        const gridreal phi2   =  r_node1[2];
        gridreal dphi   =  phi2 - phi1;
        gridreal dr     =  r2 - r1;
        real Er0 = Ax0*sin(theta)*cos(phi1) + Ay0*sin(theta)*sin(phi1) + Az0*cos(theta);
        real Er1 = Ax1*sin(theta)*cos(phi2) + Ay1*sin(theta)*sin(phi2) + Az1*cos(theta);
        real Er2 = Ax2*sin(theta)*cos(phi2) + Ay2*sin(theta)*sin(phi2) + Az2*cos(theta);
        real Er3 = Ax3*sin(theta)*cos(phi1) + Ay3*sin(theta)*sin(phi1) + Az3*cos(theta);
        real L01 =  0.5*(Ax0 + Ax1)*r1*sin(theta)*(cos(phi2)-cos(phi1)) + 0.5*(Ay0 + Ay1)*r1*sin(theta)*(sin(phi2)-sin(phi1));
        real L23 =  0.5*(Ax2 + Ax3)*r2*sin(theta)*(cos(phi1)-cos(phi2)) + 0.5*(Ay2 + Ay3)*r2*sin(theta)*(sin(phi1)-sin(phi2));
        real L12 =  0.5*(Er1 + Er2)*dr;
        real L30 = -0.5*(Er3 + Er0)*dr;
        dS = 0.5*dr*(r2 + r1)*sin(theta)*dphi;
        if (dS<0.0) dS = -dS;
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
    // Note: the number 0-3 inside [] is the index used in hybrid model
    /*          [0]- line30 -[3]       ----- r1    length of line01 =  dr
                /             \                    length of line12 =  r2*dtheta
               /               \                   length of line23 = -dr (dr < 0)
              /                 \                  length of line30 = -r1*dtheta (dtheta < 0)
           line01    Bphi(o)   line23
            /                     \
           /                       \
          /                         \
        [1]-------- line12 --------[2] ----- r2
         |                          |
       theta1                     theta2
    */
    if (d==2) {
        const gridreal r1      =  r_node0[0];
        const gridreal r2      =  r_node1[0];
        const gridreal theta1  =  r_node0[1];
        const gridreal theta2  =  r_node2[1];
        const gridreal phi     =  r_node0[2];
        gridreal dr      =  r2 - r1;         // Note: d_r does NOT have to be a constant !
        gridreal dtheta  =  theta2 - theta1; // Note: d_theta does NOT have to be a constant !
        real Er0 = Ax0*sin(theta1)*cos(phi) + Ay0*sin(theta1)*sin(phi) + Az0*cos(theta1);
        real Er1 = Ax1*sin(theta1)*cos(phi) + Ay1*sin(theta1)*sin(phi) + Az1*cos(theta1);
        real Er2 = Ax2*sin(theta2)*cos(phi) + Ay2*sin(theta2)*sin(phi) + Az2*cos(theta2);
        real Er3 = Ax3*sin(theta2)*cos(phi) + Ay3*sin(theta2)*sin(phi) + Az3*cos(theta2);
        real L01 =  0.5*(Er0 + Er1)*dr;
        real L23 = -0.5*(Er2 + Er3)*dr;
        real L12 =  0.5*(Ax1 + Ax2)*r2*cos(phi)*(sin(theta2)-sin(theta1)) + 0.5*(Ay1 + Ay2)*r2*sin(phi)*(sin(theta2)-sin(theta1))
                    +0.5*(Az1 + Az2)*r2*(cos(theta2)-cos(theta1));
        real L30 =  0.5*(Ax3 + Ax0)*r1*cos(phi)*(sin(theta1)-sin(theta2)) + 0.5*(Ay3 + Ay0)*r1*sin(phi)*(sin(theta1)-sin(theta2))
                    +0.5*(Az3 + Az0)*r1*(cos(theta1)-cos(theta2));
        dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
        if (dS<0.0) dS = -dS; // Just in case ...
        const real twice_integral = (L01 + L12 + L23 + L30)/dS;
        facedata[fs] = factor*twice_integral;
    }
}

//! (SPHERICAL) Spherical version of "NF1".
void Tgrid::Tface::sph_NF1(TNodeDataSelect ns, TFaceDataSelect fs, int d)
{
    // Three cases how to interpolate quantities from nodes to faces.
    // I.  No interpolation. We just use analytical formula to calculate values on faces. For tesing reason.
    // II. Cartesian way of interpolation.
    // II. Spherical way of interpolation
    // Analytical version. We are setting magnetic field on spherical faces.
    /*fastreal r_node0[3] = {node[0]->sph_centroid[0], node[0]->sph_centroid[1], node[0]->sph_centroid[2]};
    fastreal r_node1[3] = {node[1]->sph_centroid[0], node[1]->sph_centroid[1], node[1]->sph_centroid[2]};
    fastreal r_node2[3] = {node[2]->sph_centroid[0], node[2]->sph_centroid[1], node[2]->sph_centroid[2]};
    fastreal r_node3[3] = {node[3]->sph_centroid[0], node[3]->sph_centroid[1], node[3]->sph_centroid[2]};
    real A[3] = {0.0, 0.0, 1e-10};
    if (d==0)
       {
        const real r       =  r_node0[0];
        const real theta1  =  r_node0[1];
        const real theta2  =  r_node1[1];
        const real phi1    =  r_node0[2];
        const real phi2    =  r_node2[2];
              real theta   = (theta1 + theta2)/2;
              real phi     = (phi1 + phi2)/2;
        gridreal R[3] = {r, theta2, phi2};
        sph_transf_C2S_A(R, A);
        facedata[fs] = A[0];
        //facedata[fs] = A[2]*cos(theta2);
       }
    if (d==1)
       {
        const real r2     =  r_node0[0];
        const real r1     =  r_node3[0];
        const real theta  =  r_node0[1];
        const real phi1   =  r_node0[2];
        const real phi2   =  r_node2[2];
              real r      = (r1 + r2)/2;
              real phi    = (phi1 + phi2)/2;
        gridreal R[3] = {r2, theta, phi2};
        sph_transf_C2S_A(R, A);
        facedata[fs] = A[1];
        //facedata[fs] = -A[2]*sin(theta);
       }
    if (d==2)
       {
        const real r2      =  r_node0[0];
        const real r1      =  r_node1[0];  // Note: r1 < r2
        const real theta1  =  r_node1[1];
        const real theta2  =  r_node2[1];           // Note: the1 < the2
        const real phi     =  r_node0[2];           // Note: the1 < the2
              real r       = (r1 + r2)/2;
              real theta   = (theta1 + theta2)/2;
        gridreal R[3] = {r, theta2, phi};
        sph_transf_C2S_A(R, A);
        facedata[fs] = A[2];
        //facedata[fs] = 0.0;
       }*/
    // II. Cartesian version of NF
    /*// Face number (0,1,2)
    //         (0)  F01_dL01  (1)
    //          -------x-------
    //          |             |
    //          |             |
    // F30_dL30 x             x F12_dL12
    //          |             |
    //          |             |
    //          -------x-------
    //         (3)  F23_dL23  (2)
    real node_0_nodedata = node[0]->nodedata[ns][d];
    real node_1_nodedata = node[1]->nodedata[ns][d];
    real node_2_nodedata = node[2]->nodedata[ns][d];
    real node_3_nodedata = node[3]->nodedata[ns][d];
    const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata);
    const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata);
    const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata);
    const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata);
    const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/4;
    facedata[fs] = integral;*/
    // III. Spherical version of NF. That is just average of vector over face loop in the direction [d]
    real node_0_nodedata = node[0]->nodedata[ns][d];
    real node_1_nodedata = node[1]->nodedata[ns][d];
    real node_2_nodedata = node[2]->nodedata[ns][d];
    real node_3_nodedata = node[3]->nodedata[ns][d];
    fastreal r_node0[3] = {node[0]->sph_centroid[0], node[0]->sph_centroid[1], node[0]->sph_centroid[2]};
    fastreal r_node1[3] = {node[1]->sph_centroid[0], node[1]->sph_centroid[1], node[1]->sph_centroid[2]};
    fastreal r_node2[3] = {node[2]->sph_centroid[0], node[2]->sph_centroid[1], node[2]->sph_centroid[2]};
    fastreal r_node3[3] = {node[3]->sph_centroid[0], node[3]->sph_centroid[1], node[3]->sph_centroid[2]};
    // The face 1 goes through the nodes (in this order): 1-2-3-4. Direction of face 1 is -e_r (so to the origo)=> dBr is parall e_r
    // Note: the number 0-3 inside [] is the index used in hybrid model
    //  2[3] --- line30  ---- 1[0]     ------ theta1 (<theta2) length = r*sin(theta1)*dphi
    //  |                     |
    //  |                     |
    // line23      Br       line01
    //  |                     |
    //  |                     |
    //  3[2] --- line12  ---- 4[1]      ----- theta2           length = r*sin(theta2)*dphi
    //  |                     |
    //  |                     |
    // phi2                  phi1 (<phi2)
    if (d==0) {
        const real r       =  r_node0[0];
        const real theta1  =  r_node0[1];
        const real theta2  =  r_node1[1];
        const real phi1    =  r_node0[2];
        const real phi2    =  r_node2[2];
        const real dtheta  =  abs(theta2 - theta1);
        const real dphi    =  abs(phi2 - phi1);
        const real dL01    = r*dtheta;
        const real dL12    = r*sin(theta2)*dphi;
        const real dL23    = r*dtheta;
        const real dL30    = r*sin(theta1)*dphi;
        const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
        const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
        const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
        const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
        real L = dL01 + dL12 + dL23 + dL30;
        const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
        facedata[fs] = integral;
    }
    // The face 2 goes through the nodes (in this order): 1-5-6-2. Direction of face 1 is -e_theta => dB_theta parall e_theta
    // Note: the number 0-3 inside [] is the index used in hybrid model
    //  [0] ----- line01 ------  [1]    - r2 (> r1) length = r2*sin(theta)*dphi
    //    \                       /
    //     \                     /
    //   line30     Btheta   line12
    //       \                 /
    //       [3] --line23-- [2]         -r1  length = r1*sin(theta)*dphi
    //        |               |
    //       phi2            phi1
    if (d==1) {
        const real r2     =  r_node0[0];
        const real r1     =  r_node3[0];
        const real theta  =  r_node0[1];
        const real phi1   =  r_node0[2];
        const real phi2   =  r_node2[2];
        const real dphi   =  abs(phi2 - phi1);
        const real dr     =  abs(r2 - r1);
        const real dL01     = r2*sin(theta)*dphi;
        const real dL12     = dr;
        const real dL23     = r1*sin(theta)*dphi;
        const real dL30     = dr;
        const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
        const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
        const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
        const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
        real L = dL01 + dL12 + dL23 + dL30;
        const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
        facedata[fs] = integral;
    }
    // Note: the number 0-3 inside [] is the index used in hybrid model
    // The face 3 goes through the nodes (in this order): 1-4-8-5. Direction of face 1 is -e_phi
    //                _________ 5[0]    - theta1    length dr
    //       ________/ line51   |
    //      1[1]                |
    //      |                   |
    //     line14             line85
    //      |                   |
    //      4[2]______          |
    //                \________ 8[3]    - theta2    length dr
    //                  line 48
    //      |                  |
    //      r1                 r2
    if (d==2) {
        const real r2      =  r_node0[0];
        const real r1      =  r_node1[0];  // Note: r1 < r2
        const real theta1  =  r_node1[1];
        const real theta2  =  r_node2[1];           // Note: the1 < the2
        const real dr      =  abs(r2 - r1);         // Note: d_r does NOT have to be a constant !
        real dtheta  =  abs(theta2 - theta1); // Note: d_theta does NOT have to be a constant !
        const real dL01     = dr;
        const real dL12     = r1*dtheta;
        const real dL23     = dr;
        const real dL30     = r2*dtheta;
        const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
        const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
        const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
        const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
        real L = dL01 + dL12 + dL23 + dL30;
        const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
        facedata[fs] = integral;
    }
}

//! (SPHERICAL) Spherical version of "NF_rhoq1".
void Tgrid::Tface::sph_NF_rhoq1(int d)
{
    Tcell *c;
    real node_0_nodedata = node[0]->nodedata[NODEDATA_UE][0];
    real node_1_nodedata = node[1]->nodedata[NODEDATA_UE][0];
    real node_2_nodedata = node[2]->nodedata[NODEDATA_UE][0];
    real node_3_nodedata = node[3]->nodedata[NODEDATA_UE][0];
    // I. Cartesian version of NF
    // Face number (0,1,2)
    //         (0)  F01_dL01  (1)
    //          -------x-------
    //          |             |
    //          |             |
    // F30_dL30 x             x F12_dL12
    //          |             |
    //          |             |
    //          -------x-------
    //         (3)  F23_dL23  (2)
    const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata);
    const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata);
    const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata);
    const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata);
    const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/4;
    facedata[FACEDATA_MINUSDB] = integral;
    // II. Spherical version of NF
    /*fastreal r_node0[3] = {node[0]->sph_centroid[0], node[0]->sph_centroid[1], node[0]->sph_centroid[2]};
    fastreal r_node1[3] = {node[1]->sph_centroid[0], node[1]->sph_centroid[1], node[1]->sph_centroid[2]};
    fastreal r_node2[3] = {node[2]->sph_centroid[0], node[2]->sph_centroid[1], node[2]->sph_centroid[2]};
    fastreal r_node3[3] = {node[3]->sph_centroid[0], node[3]->sph_centroid[1], node[3]->sph_centroid[2]};
     // The face 1 goes through the nodes (in this order): 1-2-3-4. Direction of face 1 is -e_r (so to the origo)=> dBr is parall e_r
          // See Vihko p. 132:
          // Note: the number 0-3 inside [] is the index used in hybrid model
          //  2[3] --- line30  ---- 1[0]     ------ theta1 (<theta2) length = r*sin(theta1)*dphi
          //  |                     |
          //  |                     |
          // line23      Br       line01
          //  |                     |
          //  |                     |
          //  3[2] --- line12  ---- 4[1]      ----- theta2           length = r*sin(theta2)*dphi
          //  |                     |
          //  |                     |
          // phi2                  phi1 (<phi2)
    if (d==0)
     {
      const real r       =  r_node0[0];
      const real theta1  =  r_node0[1];
      const real theta2  =  r_node1[1];
      const real phi1    =  r_node0[2];
      const real phi2    =  r_node2[2];
      const real dtheta  =  theta2 - theta1;
      const real dphi    =  phi2 - phi1;
      const real dL01    = r*dtheta;
      const real dL12    = r*sin(theta2)*dphi;
      const real dL23    = r*dtheta;
      const real dL30    = r*sin(theta1)*dphi;
      const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
      const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
      const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
      const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
                   real L = dL01 + dL12 + dL23 + dL30;
      const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
      facedata[FACEDATA_MINUSDB] = integral;
     }
         // The face 2 goes through the nodes (in this order): 1-5-6-2. Direction of face 1 is -e_theta => dB_theta parall e_theta
         // Note: the number 0-3 inside [] is the index used in hybrid model
         //   [0] ----- line01 ------  [1]    - r2 (> r1) length = r2*sin(theta)*dphi
         //    \                       /
         //     \                     /
         //    line30     Btheta   line12
         //       \                 /
         //       [3] --line23-- [2]         -r1  length = r1*sin(theta)*dphi
         //        |              |
         //       phi2           phi1
    if (d==1)
     {
      const real r2     =  r_node0[0];
      const real r1     =  r_node3[0];
      const real theta1 =  r_node0[1];
      const real phi1   =  r_node0[2];
      const real phi2   =  r_node2[2];
      const real dphi   =  phi2 - phi1;
      const real dr     =  r2 - r1;
      const real dL01     = r2*sin(theta1)*dphi;
      const real dL12     = dr;
      const real dL23     = r1*sin(theta1)*dphi;
      const real dL30     = dr;
      const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
      const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
      const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
      const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
                   real L = dL01 + dL12 + dL23 + dL30;
      const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
      facedata[FACEDATA_MINUSDB] = integral;
     }
         // Note: the number 0-3 inside [] is the index used in hybrid model
         // The face 3 goes through the nodes (in this order): 1-4-8-5. Direction of face 1 is -e_phi
         //            _________ 5[0]     - theta1    length dr
         //   ________/ line51   |
         //  1[1]                |
         //  |                   |
         // line14             line85
         //  |                   |
         //  4[2]______          |
         //            \________ 8[3]    - theta2    length dr
         //              line 48
         //  |                  |
         //  r1                 r2
    if (d==2)
     {
      const real r2      =  r_node0[0];
      const real r1      =  r_node1[0];  // Note: r1 < r2
      const real theta1  =  r_node1[1];
      const real theta2  =  r_node2[1];  // Note: the1 < the2
      const real dr      =  r2 - r1;     // Note: d_r does NOT have to be a constant !
            real dtheta  =  theta2 - theta1; // Note: d_theta does NOT have to be a constant !
      const real dL01     = dr;
      const real dL12     = r1*dtheta;
      const real dL23     = dr;
      const real dL30     = r2*dtheta;
      const real F01_dL01 = 0.5*(node_0_nodedata + node_1_nodedata)*dL01;
      const real F12_dL12 = 0.5*(node_1_nodedata + node_2_nodedata)*dL12;
      const real F23_dL23 = 0.5*(node_2_nodedata + node_3_nodedata)*dL23;
      const real F30_dL30 = 0.5*(node_3_nodedata + node_0_nodedata)*dL30;
                   real L = dL01 + dL12 + dL23 + dL30;
      const real integral = (F01_dL01 + F12_dL12 + F23_dL23 + F30_dL30)/L;
      facedata[FACEDATA_MINUSDB] = integral;
     }*/
}

//! (SPHERICAL) Spherical version of "faceintpol". Here r[3] must be in spherical coordinates
void Tgrid::sph_faceintpol(shortreal r[3], TFaceDataSelect s, real result[3])
{
    int d;
    gridreal t,lowercorner[3];
    Tcell *const c = findcell(r,lowercorner); // findcell in Hybrid coordinates
    if (!c) {
        errorlog << "ERROR [Tgrid::sph_faceintpol]: findcell returned null for r=" << Tr3v(r).toString() << "\n";
        doabort();
    }
    if (c->anyrefined_face()) {
        for (d=0; d<3; d++) {
            // test: t = (r[d] - lowercorner[d])/c->sph_size[d];
            t = (r[d] - lowercorner[d])*c->invsize;
            result[d] = (1-t)*c->faceave(d,0,s) + t*c->faceave(d,1,s);
        }
    } else {
        // faster branch for fully regular case
        // This is face interpolation for spherical coordinates. We use linear interpolation F(t) = (1-t)*F1 + t*F2 but for normalized values with restpect to square F(t) = ((1-t)*F1*S1 + t*F2*S2)/S(t)
        const gridreal r_sph  = c->sph_coor[0];
        const gridreal theta1 = c->sph_coor[1];
        const gridreal theta2 = c->sph_coor_next[1];
        const gridreal dr     = c->sph_size[0];
        const gridreal dtheta = c->sph_size[1];
        const gridreal dphi   = c->sph_size[2];
        const gridreal rc     = c->sph_centroid[0];
        const gridreal dS1[3] = {c->sph_dS[0], c->sph_dS[1], c->sph_dS[2]};
        const gridreal dS2[3] = {c->sph_dS_next[0], c->sph_dS_next[1], c->sph_dS_next[2]};
        gridreal dSt[3];
        t = (r[0] - lowercorner[0])*c->invsize;
        //dSt[0] = sqr(r_sph + t*dr)*sin(theta)*dtheta*dphi;
        dSt[0] = sqr(r_sph + t*dr)*(cos(theta2)-cos(theta1))*dphi;
        result[0] = ((1-t)*c->face[0][0]->facedata[s]*dS1[0] + t*c->face[0][1]->facedata[s]*dS2[0])/dSt[0];
        t = (r[1] - lowercorner[1])*c->sph_invsizey;
        dSt[1] = rc*sin(theta1 + t*dtheta)*dr*dphi;
        result[1] = ((1-t)*c->face[1][0]->facedata[s]*dS1[1] + t*c->face[1][1]->facedata[s]*dS2[1])/dSt[1];
        t = (r[2] - lowercorner[2])*c->sph_invsizez;
        dSt[2] = rc*dr*dtheta;
        result[2] = ((1-t)*c->face[2][0]->facedata[s]*dS1[2] + t*c->face[2][1]->facedata[s]*dS2[2])/dSt[2];
        //t = (r[0] - lowercorner[0])*c->invsize;
        //dSt[0] = sqr(r_sph + t*dr)*sin(theta)*dtheta*dphi;
        //result[0] = ((1-t)*c->face[0][0]->facedata[s]/dS1[0] + t*c->face[0][1]->facedata[s]/dS2[0])*dSt[0];
        //t = (r[1] - lowercorner[1])*c->sph_invsizey;
        //dSt[1] = rc*sin(theta + t*dtheta)*dr*dphi;
        //result[1] = ((1-t)*c->face[1][0]->facedata[s]/dS1[1] + t*c->face[1][1]->facedata[s]/dS2[1])*dSt[1];
        //t = (r[2] - lowercorner[2])*c->sph_invsizez;
        //dSt[2] = rc*dr*dtheta;
        //result[2] = ((1-t)*c->face[2][0]->facedata[s]/dS1[2] + t*c->face[2][1]->facedata[s]/dS2[2])*dSt[2];
    }
    saved_cellptr = c;
}


//! (SPHERICAL) Spherical version of "FC":  Interpolate cell quantity from face quantity
void Tgrid::sph_FC(TFaceDataSelect fs, TCellDataSelect cs)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        if (j > 1 && j < ny-2) {
            c = flatindex(i,j,k);
            cells[c]->sph_FC_recursive(fs,cs,0);
        }
        if (j == 1)            {
            c = flatindex(i,j,k);
            cells[c]->sph_FC_recursive(fs,cs,20);
        }
        if (j == ny-2)         {
            c = flatindex(i,j,k);
            cells[c]->sph_FC_recursive(fs,cs,21);
        }
    }
}

//! (SPHERICAL) Another Spherical version of "FC": Interpolate cell quantity from face quantity
/*void Tgrid::sph_FC(TFaceDataSelect fs, TCellDataSelect cs)
{
	int i,j,k,c;
	ForInterior(i,j,k) {
		c = flatindex(i,j,k);
		cells[c]->sph_FC_recursive(fs,cs);
	}
}*/

//! (SPHERICAL) Spherical version of "NC_smoothing"
void Tgrid::sph_NC_smoothing()
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->sph_NC_smoothing_recursive();
    }
}

//! (SPHERICAL) Spherical version of "NC"
void Tgrid::sph_NC(TNodeDataSelect ns,TCellDataSelect cs)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->sph_NC_recursive(ns,cs);
    }
}

//! (SPHERICAL) Spherical version of "CN"
void Tgrid::sph_CN(TCellDataSelect cs, TNodeDataSelect ns)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_recursive(cs,ns);
            }
}

//! (SPHERICAL) Spherical version of "CN" + CN boundary calcultion
void Tgrid::sph_CNb(TCellDataSelect cs, TNodeDataSelect ns) // copy(to, from)
{
    int i,j,k,c;
// interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_recursive(cs,ns,1);
            }
// interpolation for the internal side nodes which touch 4 internal cells
// 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_recursive(cs,ns,6);
        }
// interpolation for the internal side nodes which touch 2 internal cells
// 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_recursive(cs,ns,12);
    }
// interpolation for the internal side nodes which touch 1 internal cell
// 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_recursive(cs,ns,8);
}

//! (SPHERICAL) Spherical version of "CN" based on cell volume
void Tgrid::sph_CN_vol(TCellDataSelect cs, TNodeDataSelect ns)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_vol_recursive(cs,ns);
            }
}

//! (SPHERICAL) Spherical version of "CN". For Ji only
void Tgrid::sph_CN_Ji(TCellDataSelect cs, TNodeDataSelect ns)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_Ji_recursive(cs,ns);
            }
}

//! (SPHERICAL) Spherical version of "CN" + CN boundary calcultion. For Ji only
void Tgrid::sph_CNb_Ji(TCellDataSelect cs, TNodeDataSelect ns) // copy(to, from)
{
    int i,j,k,c;
// interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_Ji_recursive(cs,ns,1);
            }
// interpolation for the internal side nodes which touch 4 internal cells
// 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_Ji_recursive(cs,ns,6);
        }
// interpolation for the internal side nodes which touch 2 internal cells
// 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_Ji_recursive(cs,ns,12);
    }
// interpolation for the internal side nodes which touch 1 internal cell
// 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_Ji_recursive(cs,ns,8);
}

//! (SPHERICAL) Cells averaging procedure, only for theta axis
void Tgrid::sph_YCave(TCellDataSelect cs)
{
    int i,j,k,c,d;
    int dir=0;
    for (i=0; i<nx; i++) {
        real A_from[3] = {0,0,0};
        real A_to[3]   = {0,0,0};
        real A_temp[3] = {0,0,0};
        //////////////// North pole (theta = 0) //////////////////
        // Accumulating the cell values from preboundary cell layer (for certain i and over all k)
        j = 0;
        //j = 1;
        for (k=0; k<nz; k++) {
            c = flatindex(i,j,k);
            gridreal r_from[3] = {cells[c]->sph_centroid[0], cells[c]->sph_centroid[1], cells[c]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = cells[c]->celldata[cs][d];
            sph_transf_S2C_A(r_from, A_temp);
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average cell value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz); // Use for CN type of interpolation
        //for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-2); // Use for CNb type of interpolation
        // Recording an average cell value to the boundary cell layer (for certain i and over all k)
        j = 0;
        //j = 1;
        for (k=0; k<nz; k++) {
            c = flatindex(i,j,k);
            gridreal r_to[3] = {cells[c]->sph_centroid[0], cells[c]->sph_centroid[1], cells[c]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            sph_transf_C2S_A(r_to, A_temp);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = A_temp[d];
        }
        //////////////// South pole (theta = pi) //////////////////
        for (d=0; d<3; d++) A_from[d] = 0.0;
        for (d=0; d<3; d++) A_to[d]   = 0.0;
        for (d=0; d<3; d++) A_temp[d] = 0.0;
        // Accumulating the cell values from preboundary cell layer (for certain i and over all k)
        j = ny-1;
        //j = ny-2;
        for (k=0; k<nz; k++) {
            c = flatindex(i,j,k);
            gridreal r_from[3] = {cells[c]->sph_centroid[0], cells[c]->sph_centroid[1], cells[c]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = cells[c]->celldata[cs][d];
            sph_transf_S2C_A(r_from, A_temp);
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average cell value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz); // Use for CN type of interpolation
        //for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-2); // Use for CNb type of interpolation
        // Recording an average cell value to the boundary cell layer (for certain i and over all k)
        j = ny-1;
        //j = ny-2;
        for (k=0; k<nz; k++) {
            c = flatindex(i,j,k);
            gridreal r_to[3] = {cells[c]->sph_centroid[0], cells[c]->sph_centroid[1], cells[c]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            sph_transf_C2S_A(r_to, A_temp);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = A_temp[d];
        }
    }// End of the loop "for (i=0; i<nx-1; i++)"
}

//! (SPHERICAL) Node averaging procedure, only for theta axis
void Tgrid::sph_YNave(TNodeDataSelect ns)
{
    int i,j,k,c,d;
    int dir=0;
    for (i=0; i<nx-1; i++) {
        real A_from[3] = {0,0,0};
        real A_to[3]   = {0,0,0};
        real A_temp[3] = {0,0,0};
        //////////////// North pole (theta = 0) //////////////////
        // Accumulating the nodes values from preboundary node layer (for certain i and over all k)
        //j = 1;
        j = 0;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            gridreal r_from[3] = {cells[c]->face[0][1]->node[2]->sph_centroid[0], cells[c]->face[0][1]->node[2]->sph_centroid[1], cells[c]->face[0][1]->node[2]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = cells[c]->face[dir][1]->node[2]->nodedata[ns][d];
            sph_transf_S2C_A(r_from, A_temp);
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average node value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-1);
        // Recording an average node value to the boundary node layer (for certain i and over all k)
        j = 0;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            gridreal r_to[3] = {cells[c]->face[0][1]->node[2]->sph_centroid[0], cells[c]->face[0][1]->node[2]->sph_centroid[1], cells[c]->face[0][1]->node[2]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            sph_transf_C2S_A(r_to, A_temp);
            for (d=0; d<3; d++) cells[c]->face[0][1]->node[2]->nodedata[ns][d] = A_temp[d];
        }
        //////////////// South pole (theta = pi) //////////////////
        for (d=0; d<3; d++) A_from[d] = 0.0;
        for (d=0; d<3; d++) A_to[d]   = 0.0;
        for (d=0; d<3; d++) A_temp[d] = 0.0;
        // Accumulating the nodes values from preboundary node layer (for certain i and over all k)
        j = ny-2;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            gridreal r_from[3] = {cells[c]->face[0][1]->node[2]->sph_centroid[0], cells[c]->face[0][1]->node[2]->sph_centroid[1], cells[c]->face[0][1]->node[2]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = cells[c]->face[dir][1]->node[2]->nodedata[ns][d];
            sph_transf_S2C_A(r_from, A_temp);
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average node value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-1);
        // Recording an average node value to the boundary node layer (for certain i and over all k)
        j = ny-2;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            gridreal r_to[3] = {cells[c]->face[0][1]->node[2]->sph_centroid[0], cells[c]->face[0][1]->node[2]->sph_centroid[1], cells[c]->face[0][1]->node[2]->sph_centroid[2]};
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            sph_transf_C2S_A(r_to, A_temp);
            for (d=0; d<3; d++) cells[c]->face[0][1]->node[2]->nodedata[ns][d] = A_temp[d];
        }
    }// End of the loop "for (i=0; i<nx-1; i++)"
}

//! (SPHERICAL) Node averaging procedure, only for theta axis. For Ji only.
void Tgrid::sph_YNave_Ji(TNodeDataSelect ns)
{
    int i,j,k,c,d;
    int dir=0;
    for (i=0; i<nx-1; i++) {
        real A_from[3] = {0,0,0};
        real A_to[3]   = {0,0,0};
        real A_temp[3] = {0,0,0};
        //////////////// North pole (theta = 0) //////////////////
        // Accumulating the nodes values from preboundary node layer (for certain i and over all k)
        //j = 1;
        j = 0;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) A_temp[d] = cells[c]->face[dir][1]->node[2]->nodedata[ns][d];
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average node value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-1);
        // Recording an average node value to the boundary node layer (for certain i and over all k)
        j = 0;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            for (d=0; d<3; d++) cells[c]->face[0][1]->node[2]->nodedata[ns][d] = A_temp[d];
        }
        //////////////// South pole (theta = pi) //////////////////
        for (d=0; d<3; d++) A_from[d] = 0.0;
        for (d=0; d<3; d++) A_to[d]   = 0.0;
        for (d=0; d<3; d++) A_temp[d] = 0.0;
        // Accumulating the nodes values from preboundary node layer (for certain i and over all k)
        //j = ny-3;
        j = ny-2;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) A_temp[d] = cells[c]->face[dir][1]->node[2]->nodedata[ns][d];
            for (d=0; d<3; d++) A_from[d] += A_temp[d];
        }
        // Average node value calculation
        for (d=0; d<3; d++) A_to[d] = A_from[d]/(nz-1);
        // Recording an average node value to the boundary node layer (for certain i and over all k)
        j = ny-2;
        for (k=0; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) A_temp[d] = A_to[d];
            for (d=0; d<3; d++) cells[c]->face[0][1]->node[2]->nodedata[ns][d] = A_temp[d];
        }
    }// End of the loop "for (i=0; i<nx-1; i++)"
}

//! (SPHERICAL) Interpolate node particle related variables from cell variables
void Tgrid::sph_CN_smoothing()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_smoothing_recursive();
            }
}

//! (SPHERICAL) Spherical version of "CN_cmoothing" + CN boundary calcultion
void Tgrid::sph_CNb_smoothing()
{
    int i,j,k,c;
    // interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_smoothing_recursive(1);
            }
    // interpolation for the internal side nodes which touch 4 internal cells
    // 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_smoothing_recursive(6);
        }
    // interpolation for the internal side nodes which touch 2 internal cells
    // 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_smoothing_recursive(12);
    }
    // interpolation for the internal side nodes which touch 1 internal cell
    // 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_smoothing_recursive(8);
}

//! (SPHERICAL) Spherical version of "CN_rhoq"
void Tgrid::sph_CN_rhoq()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_rhoq_recursive();
            }
}

//! (SPHERICAL) Spherical version of "CN_rhoq"
void Tgrid::sph_CNb_rhoq()
{
    int i,j,k,c;
    // interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_rhoq_recursive(1);
            }
    // interpolation for the internal side nodes which touch 4 internal cells
    // 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_rhoq_recursive(6);
        }
    // interpolation for the internal side nodes which touch 2 internal cells
    // 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_rhoq_recursive(12);
    }
    // interpolation for the internal side nodes which touch 1 internal cell
    // 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_rhoq_recursive(8);
}

//! (SPHERICAL) Cells, Nodes and Faces maps for certain r and t: used to debug boundary and pole things
void Tgrid::sph_CNF_map(TCellDataSelect cs, TNodeDataSelect ns, TFaceDataSelect fs)
{
    int i,j,k,c,d;
    unsigned long int n;
    fastreal t_cur = Params::t;
    fastreal dt    = Params::dt;
    fastreal t_map = 100.0;
    fastreal R_map = 2.0*Params::R_P;
    real Nd[3];
    real Cd[3];
    real absNd;
    real absCd;
    ofstream cnflog;
    cnflog.open("CNF_map.dat",ios_base::out);
    ///////////////////////////////////////////////////////////////////////
    if (t_map < 0.0 || t_map > Params::t_max) {
        errorlog << "ERROR [Tgrid::sph_CNF_map]: t_map = " << t_map << " that is out of limit (" << 0.0 << "," << Params::t_max << ")" << "\n";
        doabort();
    }
    if (R_map < Params::box_xmin || R_map > Params::box_xmax) {
        errorlog << "ERROR [Tgrid::sph_CNF_map]: R_map = " << R_map << " that is out of limit (" << Params::box_xmin << "," << Params::box_xmax << ")" << "\n";
        doabort();
    }
    ///////////////////////////////////////////////////////////////////////
    cnflog << scientific << showpos;
    //cnflog.precision(2);
    for (n=0; n<1; n++) {
        if (t_cur == t_map + n*dt) {
            cnflog << "\n";
            cnflog << "Cell values - B" << "\n";
            cnflog << "Node values - E" << "\n\n";
            //cnflog << "Face values - B" << "\n\n";
            cnflog << "R_map = " << R_map << "\n";
            cnflog.precision(3);
            cnflog << "t_map = " << t_cur << "\n\n";
            cnflog.precision(2);
            ///////////////////////////////////////////////////////////////////////
            i = int((R_map - x_1)/bgdx); //Cell number in r direction
            if (R_map != x_1 + i*bgdx) {
                errorlog << "ERROR [Tgrid::sph_CNF_map]: R_map(i) = " << x_1 + i*bgdx << " != R_map" << "\n";    //Checking i calculation
                doabort();
            }
            ///////////////////////////////////////////////////////////////////////
            //Marking nodes by label "thata"
            for (j=0; j<ny; j++) {
                if (j==0)             cnflog << "                 theta ";
                if (j!=0 && j < ny-1) cnflog << "          theta ";
                if (j==ny-2)          cnflog << "\n";
            }
            //Marking nodes by theta value
            for (j=0; j<ny-1; j++) {
                k=0; // there is no difference which k is used here
                c = flatindex(i,j,k);
                fastreal theta = cells[c]->face[0][1]->node[2]->sph_centroid[1];
                if (j==0)           cnflog << "               " << theta;
                if (j!=0 && j<ny-2) cnflog << "       " << theta;
                if (j==ny-2)        cnflog << "       " << theta << "\n\n";
            }
            for (k=0; k<nz; k++) {
                // Two lines without datas NOTE: n < 100
                for (j=0; j<ny; j++) {
                    if (j==0) {
                        if (n < 10) {
                            cnflog << "map0" << n << ":             |";
                        } else {
                            cnflog << "map" << n << ":             |";
                        }
                    }
                    if (j!=0 && j<ny-1) cnflog << "               |";
                    if (j == ny-1)   cnflog << "\n";
                }
                for (j=0; j<ny; j++) {
                    if (j==0) {
                        if (n < 10) {
                            cnflog << "map0" << n << ":             |";
                        } else {
                            cnflog << "map" << n << ":             |";
                        }
                    }
                    if (j!=0 && j<ny-1) cnflog << "               |";
                    if (j == ny-1)   cnflog << "\n";
                }
                // Cell data line
                for (j=0; j<ny; j++) {
                    c = flatindex(i,j,k);
                    for (d=0; d<3; d++) Cd[d] = cells[c]->celldata[cs][d];
                    fastreal phi = cells[c]->sph_centroid[2];
                    absCd = norm(Cd[0], Cd[1], Cd[2]);
                    // This block without ghost cell marker
                    //if (j == 0)              cnflog << absCd << "  |";
                    //if (j != 0 && j < ny-1)  cnflog << "   " << absCd << "   |";
                    //if (j == ny-1)           cnflog << "   " << absCd << "\n";
                    // This block with ghost cell marker
                    if (k == 0             ) {
                        if (j == 0) {
                            if (n < 10) {
                                cnflog << "map0" << n << ": G" << absCd << "  |";
                            } else {
                                cnflog << "map" << n << ": G" << absCd << "  |";
                            }
                        }
                    }
                    if (k == 0             ) {
                        if (j != 0 && j < ny-1)  cnflog << "  G" << absCd << "   |";
                    }
                    if (k == 0             ) {
                        if (j == ny-1)           cnflog << "  G" << absCd << " --- phi " << phi << "\n";
                    }
                    if (k != 0 && k != nz-1) {
                        if (j == 0) {
                            if (n < 10) {
                                cnflog << "map0" << n << ": G" << absCd << "  |";
                            } else {
                                cnflog << "map" << n << ": G" << absCd << "  |";
                            }
                        }
                    }
                    if (k != 0 && k != nz-1) {
                        if (j != 0 && j < ny-1)  cnflog << "   " << absCd << "   |";
                    }
                    if (k != 0 && k != nz-1) {
                        if (j == ny-1)           cnflog << "  G" << absCd << " --- phi " << phi << "\n";
                    }
                    if (k == nz-1          ) {
                        if (j == 0) {
                            if (n < 10) {
                                cnflog << "map0" << n << ": G" << absCd << "  |";
                            } else {
                                cnflog << "map" << n << ": G" << absCd << "  |";
                            }
                        }
                    }
                    if (k == nz-1          ) {
                        if (j != 0 && j < ny-1)  cnflog << "  G" << absCd << "   |";
                    }
                    if (k == nz-1          ) {
                        if (j == ny-1)           cnflog << "  G" << absCd << " --- phi " << phi << "\n";
                    }
                }
                // Two lines without datas NOTE: n < 100
                for (j=0; j<ny; j++) {
                    if (j==0) {
                        if (n < 10) {
                            cnflog << "map0" << n << ":             |";
                        } else {
                            cnflog << "map" << n << ":             |";
                        }
                    }
                    if (j!=0 && j<ny-1) cnflog << "               |";
                    if (j == ny-1)   cnflog << "\n";
                }
                for (j=0; j<ny; j++) {
                    if (j==0) {
                        if (n < 10) {
                            cnflog << "map0" << n << ":             |";
                        } else {
                            cnflog << "map" << n << ":             |";
                        }
                    }
                    if (j!=0 && j<ny-1) cnflog << "               |";
                    if (j == ny-1)   cnflog << "\n";
                }
                // Node data line
                if (k < nz - 1) {
                    for (j=0; j<ny-1; j++) {
                        c = flatindex(i,j,k);
                        for (d=0; d<3; d++) Nd[d] = cells[c]->face[0][1]->node[2]->nodedata[ns][d];
                        fastreal phi = cells[c]->face[0][1]->node[2]->sph_centroid[2];
                        absNd = norm(Nd[0], Nd[1], Nd[2]);
                        if (j == 0) {
                            if (n < 10) {
                                cnflog << "map0" << n << ":   ----- ";
                            } else {
                                cnflog << "map" << n << ":   ----- ";
                            }
                        }
                        if (j != 0)    cnflog << absNd << " ----- ";
                        if (j == ny-2) cnflog << absNd << " -----   ---" << " phi " << phi << "\n";
                    }
                } else return;
            }//end of k loop
        }//end of "t_cur == t_map + n*dt"
    }//end of for (n=0; n<10; n++)
    cnflog.close();
}

//! (SPHERICAL) Spherical version of "CN_donor"
void Tgrid::sph_CN_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CN_donor_recursive(cs,ns,uns,dt);
            }
}

//! (SPHERICAL) Spherical version of "CN" + CN boundary calcultion
void Tgrid::sph_CNb_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt)
{
    int i,j,k,c;
    // interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,1);
            }
    // interpolation for the internal side nodes which touch 4 internal cells
    // 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,6);
        }
    // interpolation for the internal side nodes which touch 2 internal cells
    // 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,12);
    }
    // interpolation for the internal side nodes which touch 1 internal cell
    // 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_donor_recursive(cs,ns,uns,dt,8);
}

//! (SPHERICAL) Spherical version of "accumulate_PIC_recursive"
inline gridreal Tgrid::sph_accumulate_PIC_recursive(const TBoxDef& cloudbox, Tcell *c, gridreal invvol, const shortreal v[3], real w, int popid)
{
    gridreal accum = 0.0;
    TBoxDef cellbox;
    const gridreal halfdx = 0.5*size(c);
    const gridreal halfdy = 0.5*sph_sizey(c);
    const gridreal halfdz = 0.5*sph_sizez(c);
    cellbox.lowx  = c->centroid[0] - halfdx;
    cellbox.lowy  = c->centroid[1] - halfdy;
    cellbox.lowz  = c->centroid[2] - halfdz;
    cellbox.size  = size(c);
    cellbox.sph_sizey = sph_sizey(c);
    cellbox.sph_sizez = sph_sizez(c);
    //accum = intersection_volume(cloudbox,cellbox)*invvol;
    accum = sph_intersection_volume(cloudbox,cellbox)*invvol;
    // Particle number contribution from a macroparticle to this cell
    register const real accum_w = w*accum;
    // Charge contribution from a macroparticle to this cell
    datareal charge = accum_w*Params::pops[popid]->q;
    // Add particle number contribution to the cell
    c->nc += accum_w;
    // Add charge contribution to the cell
    c->rho_q += charge;
    // Add ion current contribution to the cell
    // Here is very important to know in which coordinata are calculated velocities
    for (int d=0; d<3; d++) {
        c->celldata[CELLDATA_Ji][d] += charge*v[d];
    }
#ifdef SAVE_POPULATION_AVERAGES
    if (Params::averaging == true) {
        c->pop_ave_n[popid] += accum_w;
        c->pop_ave_vx[popid] += accum_w*v[0];
        c->pop_ave_vy[popid] += accum_w*v[1];
        c->pop_ave_vz[popid] += accum_w*v[2];
    }
#endif
    return accum;
}

//! (SPHERICAL) Spherical version of "accumulate_PIC"
void Tgrid::sph_accumulate_PIC(shortreal r[3], const shortreal v[3], real w, int popid)
{
    gridreal centroid[3];
    Tcell *const c = findcell(r,centroid);
    // This block corrects the position of macroparticle in the simulation box. If
    // cloud center (not particle) leave the region
    // [x_min + 0.5dx, x_max - 0.5dx; y_min + 0.5dy, y_max - 0.5dy; z_min + 0.5sph_dz, z_max - 0.5sph_dz]
    // we move cloud (not particle) inside the region. This algorithm is only true for the clouds, which
    // have the same shape as the particles/
    if (r[0] < Params::box_xmin + Params::dx/2) {
        r[0] = Params::box_xmin + Params::dx/2;   //near x_min
    }
    if (r[0] > Params::box_xmax - Params::dx/2) {
        r[0] = Params::box_xmax - Params::dx/2;   //near x_max
    }
    if (r[1] < Params::box_ymin + Params::sph_dy/2) {
        r[1] = Params::box_ymin + Params::sph_dy/2;   //near y_min
    }
    if (r[1] > Params::box_ymax - Params::sph_dy/2) {
        r[1] = Params::box_ymax - Params::sph_dy/2;   //near y_max
    }
    if (r[2] < Params::box_zmin + Params::sph_dz/2) {
        r[2] = Params::box_zmin + Params::sph_dz/2;   //near z_min
    }
    if (r[2] > Params::box_zmax - Params::sph_dz/2) {
        r[2] = Params::box_zmax - Params::sph_dz/2;   //near z_max
    }
    if (!c) {
        // Do not abort if no cell is found. Particle removed in pass_with_relocate afterwards.
        errorlog << "ERROR [Tgrid::sph_accumulate_PIC]: findcell returned null for r="
                 << Tr3v(r).toString() << "\n"
                 << "   v=" << Tr3v(v).toString()
                 << " popIdStr=" << Params::pops[popid]->getIdStr() << endl;
        return; //doabort();
    }
    //const gridreal halfdx = 0.5*size(c);
    gridreal halfd[3] = {0.5*size(c), 0.5*sph_sizey(c), 0.5*sph_sizez(c)};
    bool movetoright[3];        // false if move to left, true if move to right, along dimension d
    for (int d=0; d<3; d++) {
        //centroid[d]+= halfdx;
        centroid[d]+= halfd[d];
        //if (r[d]==centroid[d]) errorlog << "Problem:::"<<endl;
        movetoright[d] = (r[d] > centroid[d]);
    }
    Tcell *C[8];        // eight cells that define the stencil, some of them can be refined
    C[0] = c;
    C[1] = moveto(c,0,movetoright[0]);      // move along x, right or left
    C[2] = moveto(c,1,movetoright[1]);      // move along y, right or left
    C[3] = moveto(C[1],1,movetoright[1]);
    C[4] = moveto(c,2,movetoright[2]);      // move along z, right or left
    C[5] = moveto(C[4],0,movetoright[0]);
    C[6] = moveto(C[4],1,movetoright[1]);
    C[7] = moveto(C[5],1,movetoright[1]);
    // Check if all cells are of the same level. This is a necessary requirement for "regularity" (see below)
    int a;
    //bool all_same_level = true;
    const int clevel = c->level;
    for (a=1; a<8; a++)
        if (C[a]->level != clevel) {
            //all_same_level = false;
            break;
        }
    // Replace duplicate cells by null pointers. Duplicates can occur if the neighbour is coarser.
    // Cells 0,1,2,4 cannot be duplicates anyway so exclude them from the search.
    int b;
    bool found;
    // We want 'a' to range 3,5,6,7. It is most efficient to write it all out explicitly.
    a = 3;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 5;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 6;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    a = 7;
    found = false;
    for (b=0; b<a; b++) if (C[a] == C[b]) {
            found = true;
            break;
        }
    if (found) C[a] = 0;
    //const gridreal cloudsize0 = size(c);
    // initially, the cloud is of the same size as cell c
    const gridreal cloudsizex0 = size(c);
    const gridreal cloudsizey0 = sph_sizey(c);
    const gridreal cloudsizez0 = sph_sizez(c);
    gridreal cloudsize[3];
    cloudsize[0] = cloudsizex0;
    cloudsize[1] = cloudsizey0;
    cloudsize[2] = cloudsizez0;
    const gridreal halfcloudsize[3] = {0.5*cloudsize[0], 0.5*cloudsize[1], 0.5*cloudsize[2]};
    const gridreal invvol = 1.0/(cloudsize[0]*cloudsize[1]*cloudsize[2]);
    //const gridreal invvol = 1.0/sph_dV(c);
    gridreal accum = 0.0;
    TBoxDef cloudbox;
    cloudbox.lowx  = r[0] - halfcloudsize[0];
    cloudbox.lowy  = r[1] - halfcloudsize[1];
    cloudbox.lowz  = r[2] - halfcloudsize[2];
    cloudbox.size  = cloudsize[0];
    cloudbox.sph_sizey = cloudsize[1];
    cloudbox.sph_sizez = cloudsize[2];
    for (a=0; a<8; a++) if (C[a]) {
            accum+= sph_accumulate_PIC_recursive(cloudbox,C[a],invvol,v,w,popid);
        }
}

//! (SPHERICAL) Spherical version of "finalize_accum_recursive"
void Tgrid::sph_finalize_accum_recursive(Tcell *c)
{
    if (c->haschildren) {
        for (int ch=0; ch<8; ch++) sph_finalize_accum_recursive(c->child[0][0][ch]);
    } else {
        // Now volume is calculated in spherical coordinates. It is described in Grid.h
        // Cell volume
        gridreal invvol = 1.0/sph_dV(c);
        //gridreal invvol = 1.0/size(c)/size(c)/size(c);
        // Normalize particle number in the cell by the cell volume -> particle number density
        c->nc *= invvol;
        // Normalize charge in the cell by the cell volume -> charge density
        c->rho_q *= invvol;
        //////////////////////// Ji accumulation /////////////////////////
        // Normalize ion current in the cell by the cell volume -> ion current density
        // For flat front propagation averaging of vector values should be done in Cartesian
        // coordinates, otherwise mistake ocuures.
        // V value comes from sph_PropagateV where it could be stored either in Hybrid or in
        // Cartesian form. It is accumulated in sph_accumulate_PIC, if it is stored in
        // Cartesian coordinates we need to transform vector back to Hybrid coordinates to
        // use it in another function as a Hybrid vector.
        gridreal Ji[3] = {c->celldata[CELLDATA_Ji][0], c->celldata[CELLDATA_Ji][1], c->celldata[CELLDATA_Ji][2]};
        // Ji comes in Cartesian coordinates
        for (int d=0; d<3; d++) {
            c->celldata[CELLDATA_Ji][d] = Ji[d]*invvol;
        }
        // Add background charge density
        c->rho_q += c->rho_q_bg;
        // Check charge density min value (CONSTRAINT)
        if ( c->rho_q < Params::rho_q_min ) {
            c->rho_q = Params::rho_q_min;
#ifndef NO_DIAGNOSTICS
            // Increase counter
            Tgrid::fieldCounter.cutRateRhoQ += 1.0;
#endif
        }
        // Averaging
        if (Params::averaging == true) {
            c->ave_nc+= c->nc;
            for (int d=0; d<3; d++) {
                if (c->neighbour[d][1] == 0) {
                    continue;
                }
                if (c->isrefined_face(d,1)) {
                    for (int f=0; f<4; f++) {
                        c->refintf[d][1]->face[f]->facedata[FACEDATA_AVEB] += c->refintf[d][1]->face[f]->facedata[FACEDATA_B];
                    }
                } else {
                    c->face[d][1]->facedata[FACEDATA_AVEB] += c->face[d][1]->facedata[FACEDATA_B];
                }
            }
        }
    }
}

//! (SPHERICAL) Spherical version of "finalize_accum": Divide each nc, rho_q, CELLDATA_Ji field in the cells by the volume of the cell.
void Tgrid::sph_finalize_accum()
{
    if (Params::sph_BC_use_ghost_cell == 0) sph_Neumann_rhoq_0();
    if (Params::sph_BC_use_ghost_cell == 1) sph_Neumann_rhoq();
    int i,j,k,c;
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        sph_finalize_accum_recursive(cells[c]);
    }
// Averaging
    if (Params::averaging == true) {
        ave_ntimes++;
    }
}

//! (SPHERICAL) copy celldata with SCS transformations.
void Tgrid::copy_celldata_SCS(int cTo, int cFrom, TCellDataSelect cs)
{
    int d;
    gridreal A_from[3];
    gridreal A_to[3];
    //spherical coordinates of the points to and from
    // NB! Here must be centoids, not coordinates but centroids, overwise error then H2S and S2H transformations
    gridreal r_from[3] = {cells[cFrom]->sph_centroid[0], cells[cFrom]->sph_centroid[1], cells[cFrom]->sph_centroid[2]};
    gridreal r_to[3]   = {cells[cTo]->sph_centroid[0], cells[cTo]->sph_centroid[1], cells[cTo]->sph_centroid[2]};
    for (d=0; d<3; d++) A_from[d] = cells[cFrom]->celldata[cs][d];
    //transformatio from spherical coordinates to Cartesian for A_from vectors
    sph_transf_S2C_A1(r_from, A_from);
    for (d=0; d<3; d++) A_to[d] = A_from[d];
    //transformatio from spherical coordinates to hybrid for A_to vectors
    sph_transf_C2S_A1(r_to, A_to);
    for (d=0; d<3; d++) cells[cTo]->celldata[cs][d] = A_to[d];
}

//! (SPHERICAL) Spherical version of "copy_rhoq"
void Tgrid::sph_copy_rhoq(int cTo, int cFrom)
{
    gridreal V_from = cells[cFrom]->sph_dV;
    gridreal V_to   = cells[cTo]->sph_dV;
    gridreal koeff  = V_from/V_to;

//cells[cTo]->rho_q = cells[cFrom]->rho_q*koeff;
    cells[cTo]->rho_q = cells[cFrom]->rho_q;
}

//! (SPHERICAL) Spherical version of "copy_smoothing"
void Tgrid::sph_copy_smoothing(int cTo, int cFrom)
{
    int d;
    // Ji - Cartesian vector, we don't need to transform them
    for (d=0; d<3; d++) cells[cTo]->celldata[CELLDATA_Ji][d] = cells[cFrom]->celldata[CELLDATA_Ji][d];//ion current density
    cells[cTo]->rho_q = cells[cFrom]->rho_q;//charge density
    cells[cTo]->nc = cells[cFrom]->nc;//number density
}

//! (SPHERICAL) Spherical version of "copy_smoothing". Set 0 to the ghost cells
void Tgrid::sph_copy_smoothing_0(int cTo)
{
    int d;
    for (d=0; d<3; d++) {
        cells[cTo]->celldata[CELLDATA_Ji][d] = 0.0; //ion current density
    }
    cells[cTo]->rho_q = 0.0; //charge density
    cells[cTo]->nc    = 0.0; //number density
}

//! (SPHERICAL) Spherical version of "calc_ue"
void Tgrid::sph_calc_ue(void)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->sph_calc_ue_recursive();
    }
}

//! (SPHERICAL) Temporal function. Spherical version of "calc_ue". We set Ue manualy
void Tgrid::sph_calc_ue_app(void)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        //ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->sph_calc_ue_app_recursive();
    }
}

//! (SPHERICAL) Spherical version of "calc_node_E": Calculate electric field at the nodes (E = -U_e x (B + B_0) + eta*J)
void Tgrid::sph_calc_node_E(void)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_calc_node_E_recursive();
            }
}

//! (SPHERICAL) Spherical version of "calc_cell_E": Calculate electric field at the center of the cell (E = -U_e x (B + B_0) (- Delta P) to be added
void Tgrid::sph_calc_cell_E(void)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->sph_calc_cell_E_recursive();
    }
}

//! (SPHERICAL) Spherical node boundary conditions
void Tgrid::sph_Node_BC(TNodeDataSelect ns)
{
    int i,j,k,c;

// Radial propagation
    if (Params::sph_propagation_type == 0) {
        i=2;
        for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_Node_BC_recursive(ns);
            }
    }
// Flat front propagation
    if (Params::sph_propagation_type == 1) {
        i=nx-2;
        for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_Node_BC_recursive(ns);
            }
    }
}

//! (SPHERICAL) Spherical version of "FaceCurl"
void Tgrid::sph_FaceCurl(TNodeDataSelect ns, TFaceDataSelect fs, real factor)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->sph_FaceCurl_recursive(ns,fs,0, factor);
                if (i > 0 && k > 0) cells[c]->sph_FaceCurl_recursive(ns,fs,1, factor);
                if (i > 0 && j > 0) cells[c]->sph_FaceCurl_recursive(ns,fs,2, factor);
            }
}

//! (SPHERICAL) New Spherical version of "FaceCurl". Cartesian approach.
void Tgrid::sph_FaceCurl_2(TNodeDataSelect ns, TFaceDataSelect fs, real factor)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->sph_FaceCurl_2_recursive(ns,fs,0, factor);
                if (i > 0 && k > 0) cells[c]->sph_FaceCurl_2_recursive(ns,fs,1, factor);
                if (i > 0 && j > 0) cells[c]->sph_FaceCurl_2_recursive(ns,fs,2, factor);
            }
}

//! (SPHERICAL) Spherical version of "NF"
void Tgrid::sph_NF(TNodeDataSelect ns, TFaceDataSelect fs)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->sph_NF_recursive(ns,fs,0);
                if (i > 0 && k > 0) cells[c]->sph_NF_recursive(ns,fs,1);
                if (i > 0 && j > 0) cells[c]->sph_NF_recursive(ns,fs,2);
            }
}

//! (SPHERICAL) Spherical version of "NF_rhoq"
void Tgrid::sph_NF_rhoq()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                if (j > 0 && k > 0) cells[c]->sph_NF_rhoq_recursive(0);
                if (i > 0 && k > 0) cells[c]->sph_NF_rhoq_recursive(1);
                if (i > 0 && j > 0) cells[c]->sph_NF_rhoq_recursive(2);
            }
}

//! (SPHERICAL) Spherical version of "CalcGradient_rhoq"
void Tgrid::sph_CalcGradient_rhoq(void)
{
    int i,j,k;
    ForInterior(i,j,k) {
        cells[flatindex(i,j,k)]->sph_CalcGradient_rhoq_recursive();
    }
}

//! (SPHERICAL) Periodic version of "Neumann". Periodic boundary condition along y and z
void Tgrid::sph_Neumann_periodic(TCellDataSelect cs) // copy(to, from)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i+1,j,k),cs);
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i-1,j,k),cs);
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,0,k),flatindex(i,ny-2,k),cs); // periodic BC
    // +Y boundary
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,ny-1,k),flatindex(i,1,k),cs); // periodic BC
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata(flatindex(i,j,0),flatindex(i,j,nz-2),cs); // periodic BC
    // +Z boundary
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata(flatindex(i,j,nz-1),flatindex(i,j,1),cs); // periodic BC
}

//! (SPHERICAL) Spherical version of "Neumann"
void Tgrid::sph_Neumann(TCellDataSelect cs)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i+1,j,k),cs);
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i-1,j,k),cs);
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i,j+1,k),cs);
    //copy_celldata(flatindex(i,0,k),flatindex(i,ny-2,k),cs);
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata(flatindex(i,j,k),flatindex(i,j-1,k),cs);
    //copy_celldata(flatindex(i,ny-1,k),flatindex(i,1,k),cs);
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    //k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            //copy_celldata(flatindex(i,j,0),flatindex(i,j,1),cs);
            copy_celldata(flatindex(i,j,0),flatindex(i,j,nz-2),cs); // cyclic BC
    // +Z boundary
    //k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            //copy_celldata(flatindex(i,j,nz-1),flatindex(i,j,nz-2),cs);
            copy_celldata(flatindex(i,j,nz-1),flatindex(i,j,1),cs); // cyclic BC
}

//! (SPHERICAL) Spherical version of "Neumann" with transformations SCS
void Tgrid::sph_Neumann_SCS(TCellDataSelect cs) // copy(to, from)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i+1,j,k),cs);
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i-1,j,k),cs);
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i,j+1,k),cs);
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i,j-1,k),cs);
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i,j,k+1),cs);
    //copy_celldata_SCS(flatindex(i,j,0),flatindex(i,j,nz-2),cs); // cyclic BC
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            copy_celldata_SCS(flatindex(i,j,k),flatindex(i,j,k-1),cs);
    //copy_celldata_SCS(flatindex(i,j,nz-1),flatindex(i,j,1),cs); // cyclic BC
}

//! (SPHERICAL) Set 0 to the ghost cells
void Tgrid::sph_Neumann_0(TCellDataSelect cs)
{
    int i,j,k,d,c;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->celldata[cs][d] = 0.0;
        }
}

//! (SPHERICAL) Spherical version of "Neumann_rhoq"
void Tgrid::sph_Neumann_rhoq() // copy (to, from)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i+1,j,k));
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i-1,j,k));
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i,j+1,k));
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i,j-1,k));
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i,j,nz-2)); // cyclic BC
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_rhoq(flatindex(i,j,k),flatindex(i,j,1));    // cyclic BC
}

//! (SPHERICAL) Spherical version of "Neumann_rhoq". Set 0 to the ghost cells
void Tgrid::sph_Neumann_rhoq_0()
{
    int i,j,k,d,c;
    // -X boundary
    i = 0;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
    // +X boundary
    i = nx-1;
    for (j=1; j<ny-1; j++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
    //--------------------------------------------------------
    // -Y boundary. From now on, i extends over all points.
    j = 0;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=1; k<nz-1; k++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
    //--------------------------------------------------------
    // -Z boundary. Fron now on, both i and j extend over all points.
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
            c = flatindex(i,j,k);
            for (d=0; d<3; d++) cells[c]->rho_q = 0.0;
        }
}

//! (SPHERICAL) Spherical version of "Neumann_smoothing".
void Tgrid::sph_Neumann_smoothing() // copy (to, from)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=0; j<ny; j++) for (k=0; k<nz; k++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i+1,j,k));
    // +X boundary
    i = nx-1;
    for (j=0; j<ny; j++) for (k=0; k<nz; k++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i-1,j,k));
    //--------------------------------------------------------
    // -Y boundary
    j = 0;
    for (i=0; i<nx; i++) for (k=0; k<nz; k++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i,j+1,k));
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=0; k<nz; k++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i,j-1,k));
    //--------------------------------------------------------
    // -Z boundary
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i,j,nz-2)); // cyclic BC
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_smoothing(flatindex(i,j,k),flatindex(i,j,1));    // cyclic BC
}

//! (SPHERICAL) Spherical version of "Neumann_smoothing". Set 0 to the ghost cells
void Tgrid::sph_Neumann_smoothing_0() // copy (to, from)
{
    int i,j,k;
    // -X boundary
    i = 0;
    for (j=0; j<ny; j++) for (k=0; k<nz; k++)
            sph_copy_smoothing_0(flatindex(i,j,k));
    // +X boundary
    i = nx-1;
    for (j=0; j<ny; j++) for (k=0; k<nz; k++)
            sph_copy_smoothing_0(flatindex(i,j,k));
    //--------------------------------------------------------
    // -Y boundary
    j = 0;
    for (i=0; i<nx; i++) for (k=0; k<nz; k++)
            sph_copy_smoothing_0(flatindex(i,j,k));
    // +Y boundary
    j = ny-1;
    for (i=0; i<nx; i++) for (k=0; k<nz; k++)
            sph_copy_smoothing_0(flatindex(i,j,k));
    //--------------------------------------------------------
    // -Z boundary
    k = 0;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_smoothing_0(flatindex(i,j,k));
    // +Z boundary
    k = nz-1;
    for (i=0; i<nx; i++) for (j=0; j<ny; j++)
            sph_copy_smoothing_0(flatindex(i,j,k));
}

//! (SPHERICAL) Spherical version of "smoothing".
void Tgrid::sph_smoothing()
{
    for (int n=0; n<Params::densitySmoothingNumber; n++) {
        if(Params::sph_BC_use_ghost_cell == 0) {
            // if we use sph_CNb interpolation we must not use Neumann boundary function. We olso neet to comment Neumann_rhoq line in sph_finalize_accum
            //set up Neumann boundary(can also be other boundary condition) for the particle related quantities
            sph_Neumann_smoothing_0();
            //Cell to Node interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
            sph_CNb_smoothing();
            //Node to Cell interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
            sph_NC_smoothing();
        } else if(Params::sph_BC_use_ghost_cell == 1) {
            // if we use sph_CNb interpolation we must not use Neumann boundary function. We olso neet to comment Neumann_rhoq line in sph_finalize_accum
            //set up Neumann boundary(can also be other boundary condition) for the particle related quantities
            sph_Neumann_smoothing();
            //Cell to Node interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
            sph_CN_smoothing();
            //Node to Cell interpolation. NODEDATA_UE[0] and NODEDATA_J are used as temporary storage place for node values of rho_q and VQ
            sph_NC_smoothing();
        }
    }
    if (Params::densitySmoothingNumber > 0) {
        if (Params::sph_BC_use_ghost_cell == 0) sph_Neumann_smoothing_0(); //set up Neumann boundary(can also be other boundary condition) for the particle related quantities
        else if (Params::sph_BC_use_ghost_cell == 1) sph_Neumann_smoothing(); //set up Neumann boundary(can also be other boundary condition) for the particle related quantities
    }
    return;
}

//! (SPHERICAL) Spherical version of "smoothing_E". Note!! CELLDATA_TEMP1 is used as temporary storage space.
void Tgrid::sph_smoothing_E()
{
    for (int n=0; n<Params::electricFieldSmoothingNumber; n++) {
        if (Params::sph_BC_use_ghost_cell == 0) {
            sph_NC(NODEDATA_E, CELLDATA_TEMP1);////Node to Cell interpolation. NODEDATA_E is interpolated to CELLDATA_TEMP1. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
            sph_Neumann_0(CELLDATA_TEMP1);//set up Neumann boundary(can also be other boundary condition) for CELLDATA_TEMP1 (actually saves CELLDATA_E)
            sph_CNb(CELLDATA_TEMP1, NODEDATA_E);//Cell to Node interpolation for E. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
        } else if(Params::sph_BC_use_ghost_cell == 1) {
            sph_NC(NODEDATA_E, CELLDATA_TEMP1);////Node to Cell interpolation. NODEDATA_E is interpolated to CELLDATA_TEMP1. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
            sph_Neumann_SCS(CELLDATA_TEMP1);//set up Neumann boundary(can also be other boundary condition) for CELLDATA_TEMP1 (actually saves CELLDATA_E)
            sph_CN(CELLDATA_TEMP1, NODEDATA_E);//Cell to Node interpolation for E. CELLDATA_TEMP1 is used as temporary storage place for cell values of E.
        }
    }
    return;
}

//! (SPHERICAL) Spherical version of "boundarypass": Here we add the possibility to be dependent on coordinates
void Tgrid::sph_boundarypass(int dim, bool toRight, void (*funcB)(datareal cdata[NCELLDATA][3], gridreal sph_centroid[3], int))
{
    int i,j,k;
    switch (dim) {
    case 0: // x-boundaries
        i = toRight ? nx-1 : 0;
        for (j=0; j<ny; j++) {
            for (k=0; k<nz; k++) {
                TCellPtr c = cells[flatindex(i,j,k)];
                gridreal centroid_temp[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
                (*funcB)(cells[flatindex(i,j,k)]->celldata, centroid_temp, 0);
            }
        }
        break;
    case 1: // y-boundaries
        j = toRight ? ny-1 : 0;
        for (i=0; i<nx; i++) {
            for (k=0; k<nz; k++) {
                TCellPtr c = cells[flatindex(i,j,k)];
                gridreal centroid_temp[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
                (*funcB)(cells[flatindex(i,j,k)]->celldata, centroid_temp, 1);
            }
        }
        break;
    case 2: // z-boundaries
        k = toRight ? nz-1 : 0;
        for (i=0; i<nx; i++) {
            for (j=0; j<ny; j++) {
                TCellPtr c = cells[flatindex(i,j,k)];
                gridreal centroid_temp[3] = {c->sph_centroid[0], c->sph_centroid[1], c->sph_centroid[2]};
                (*funcB)(cells[flatindex(i,j,k)]->celldata, centroid_temp, 2);
            }
        }
        break;
    }
}

//! (SPHERICAL) Spherical version of "cellintpol_fluid". Calculate fluid parameters in the cell: n, avg(v) and P=nT (use particle mass of the population popID[0])
void Tgrid::Tcell::sph_cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, vector<int> popId)
{
    n = plist.calc_mass(popId)/(Params::pops[popId[0]]->m*sph_dV);
    vx = vy = vz = 0;
    plist.calc_avev(vx, vy, vz, popId);
    const real T = 0.5*plist.calc_avemv2(vx, vy, vz, popId);
    P = n*T;
}

//! (SPHERICAL) Spherical version of "addparticle". Add particle into the grid with given parameters
void Tgrid::sph_addparticle(shortreal x, shortreal y, shortreal z, shortreal vx,shortreal vy,shortreal vz, shortreal w, int popid)
{
    const shortreal r[3] = {x,y,z};
    TCellPtr c = findcell(r);
    if (!c) {
        errorlog << "WARNING: Tgrid::addparticle" << Tr3v(r).toString()
                 << " idStr=" << Params::pops[popid]->getIdStr()
                 << " out of box (not created)\n";
        return;
    }
    //gridreal theta = c->sph_centroid[1];
    //c->plist.add(x,y,z,vx,vy,vz,w*sin(theta),popid);
    c->plist.add(x,y,z,vx,vy,vz,w,popid);
    n_particles++;
}

//! (SPHERICAL) Spherical version of "generate_random_point". Generate random point inside cell
void Tgrid::Tcell::sph_generate_random_point(gridreal r[3])
{
    int d;
    gridreal siz[3] = {size, sph_sizey, sph_sizez};
    for (d=0; d<3; d++) {
        r[d]= centroid[d] + (1-Params::box_eps)*siz[d]*(uniformrnd()-0.5);
    }
}

//! (SPHERICAL) Spherical version of "boundary_faces". Copies J to ghost faces, so J can be calculated for nodes in contact with ghosts in the same way as for inner nodes.
void Tgrid::sph_boundary_faces(TFaceDataSelect cs)
{
    int i,j,k;
// X - wall ---------------------------------------
    i = 0;
    for (j=0; j<ny; j++) for(k=0; k<nz; k++) { // -X-wall
            if (cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i+1,j,k)]->face[1][1]->facedata[cs];
            if (cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i+1,j,k)]->face[2][1]->facedata[cs];
        }
    i = nx-1;
    for (j=0; j<ny; j++) for(k=0; k<nz; k++) { // +X-wall
            if (cells[flatindex(i,j,k)]->face[1][1])
                cells[flatindex(i,j,k)]->face[1][1]->facedata[cs] = cells[flatindex(i-1,j,k)]->face[1][1]->facedata[cs];
            if (cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i-1,j,k)]->face[2][1]->facedata[cs];
        }
// Y - wall ---------------------------------------
    j = 0;
    for (i=0; i<nx; i++) for(k=0; k<nz; k++) { // -Y-wall
            if (cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j+1,k)]->face[0][1]->facedata[cs];
            if (cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,j+1,k)]->face[2][1]->facedata[cs];
        }
    j = ny-1;
    for (i=0; i<nx; i++) for(k=0; k<nz; k++) { // +Y-wall
            if (cells[flatindex(i,j,k)]->face[0][1])
                cells[flatindex(i,j,k)]->face[0][1]->facedata[cs] = cells[flatindex(i,j-1,k)]->face[0][1]->facedata[cs];
            if (cells[flatindex(i,j,k)]->face[2][1])
                cells[flatindex(i,j,k)]->face[2][1]->facedata[cs] = cells[flatindex(i,j-1,k)]->face[2][1]->facedata[cs];
        }
// Z - wall ---------------------------------------
//k = 0;
    for (i=0; i<nx; i++) for(j=0; j<ny; j++) { // -Z-wall // cyclic BC
            if (cells[flatindex(i,j,0)]->face[0][1])
                cells[flatindex(i,j,0)]->face[0][1]->facedata[cs] = cells[flatindex(i,j,nz-2)]->face[0][1]->facedata[cs];
            if (cells[flatindex(i,j,0)]->face[1][1])
                cells[flatindex(i,j,0)]->face[1][1]->facedata[cs] = cells[flatindex(i,j,nz-2)]->face[1][1]->facedata[cs];
        }
//k = nz-1;
    for (i=0; i<nx; i++) for(j=0; j<ny; j++) { // +Z-wall // cyclic BC
            if (cells[flatindex(i,j,nz-1)]->face[0][1])
                cells[flatindex(i,j,nz-1)]->face[0][1]->facedata[cs] = cells[flatindex(i,j,1)]->face[0][1]->facedata[cs];
            if (cells[flatindex(i,j,nz-1)]->face[1][1])
                cells[flatindex(i,j,nz-1)]->face[1][1]->facedata[cs] = cells[flatindex(i,j,1)]->face[1][1]->facedata[cs];
        }
}

//! (SPHERICAL) Spherical version of "jcomponent"
real Tgrid::Tnode::sph_jcomponent(int edir)
{
    // n->cell[0][0][0] = cells[c];
    // n->cell[1][0][0] = cells[flatindex(i+1,j,  k  )];
    // n->cell[0][1][0] = cells[flatindex(i,  j+1,k  )];
    // n->cell[1][1][0] = cells[flatindex(i+1,j+1,k  )];
    // n->cell[0][0][1] = cells[flatindex(i,  j,  k+1)];
    // n->cell[1][0][1] = cells[flatindex(i+1,j,  k+1)];
    // n->cell[0][1][1] = cells[flatindex(i,  j+1,k+1)];
    // n->cell[1][1][1] = cells[flatindex(i+1,j+1,k+1)];
    real bf1, bf2, bf3, bf4; // B on faces for line integral.
    real invmu0 = 1./Params::mu_0;
    real epos=0.;
    real eneg=0.;
    real intp;
    bool posf=false;
    bool negf=false;
    switch (edir) {
    case 0 :  // x-edges
        if (cell[1][0][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][0][1] or cell[1][0][1]!=cell[1][1][1]) { // positive side edge exists (it might not because of refinement)?
            // length of line 4 = -r*dtheta
            // length of line 2 =  r*dtheta (dtheta < 0)
            // length of line 3 = -r*sin(theta2)*dphi
            // length of line 1 =  r*sin(theta1)*dphi (dphi < 0)
            real r      = cell[1][0][0]->sph_centroid[0];
            real theta1 = cell[1][0][0]->sph_centroid[1];
            real theta2 = cell[1][1][0]->sph_centroid[1];
            real phi1   = cell[1][0][0]->sph_centroid[2];
            real phi2   = cell[1][0][1]->sph_centroid[2];
            real dtheta = theta2 - theta1;
            real dphi   = phi2 - phi1;
            real dL1 = -r*sin(theta1)*dphi;
            real dL2 =  r*dtheta;
            real dL3 =  r*sin(theta2)*dphi;
            real dL4 = -r*dtheta;
            //real dS00 = cell[1][0][0]->sph_dS_centroid[0];
            //real dS01 = cell[1][0][1]->sph_dS_centroid[0];
            //real dS10 = cell[1][1][0]->sph_dS_centroid[0];
            //real dS11 = cell[1][1][1]->sph_dS_centroid[0];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS = r*r*sin(theta1)*dtheta*dphi;
            if (dS<0.0) dS = -dS;
            posf=true;
            bf1 = (cell[1][0][0]!=cell[1][0][1]) ? cell[1][0][0]->faceave(2,1,FACEDATA_B)*dL1 : cell[1][0][0]->celldata[CELLDATA_B][2]*dL1;
            bf2 = (cell[1][0][0]!=cell[1][1][0]) ? cell[1][0][0]->faceave(1,1,FACEDATA_B)*dL2 : cell[1][0][0]->celldata[CELLDATA_B][1]*dL2;
            bf3 = (cell[1][1][0]!=cell[1][1][1]) ? cell[1][1][0]->faceave(2,1,FACEDATA_B)*dL3 : cell[1][1][0]->celldata[CELLDATA_B][2]*dL3;
            bf4 = (cell[1][1][1]!=cell[1][0][1]) ? cell[1][1][1]->faceave(1,0,FACEDATA_B)*dL4 : cell[1][1][1]->celldata[CELLDATA_B][1]*dL4;
            epos = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        if (cell[0][0][0]!=cell[0][1][0] or cell[0][1][0]!=cell[0][0][1] or cell[0][0][1]!=cell[0][1][1]) { // negative side edge exists?
            // length of line 1 =  r2*sin(theta)*dphi
            // length of line 3 = -r1*sin(theta)*dphi (dphi < 0)
            // length of line 2 = -dr (dr < 0)
            // length of line 4 =  dr
            real r      = cell[0][0][0]->sph_centroid[0];
            real theta1 = cell[0][0][0]->sph_centroid[1];
            real theta2 = cell[0][1][0]->sph_centroid[1];
            real phi1   = cell[0][0][0]->sph_centroid[2];
            real phi2   = cell[0][0][1]->sph_centroid[2];
            real dtheta = theta2 - theta1;
            real dphi   = phi2 - phi1;
            real dL1 = -r*sin(theta1)*dphi;
            real dL2 =  r*dtheta;
            real dL3 =  r*sin(theta2)*dphi;
            real dL4 = -r*dtheta;
            //real dS00 = cell[0][0][0]->sph_dS_centroid[0];
            //real dS01 = cell[0][0][1]->sph_dS_centroid[0];
            //real dS10 = cell[0][1][0]->sph_dS_centroid[0];
            //real dS11 = cell[0][1][1]->sph_dS_centroid[0];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS = r*r*sin(theta1)*dtheta*dphi;
            if (dS<0.0) dS = -dS;
            negf=true;
            bf1 = (cell[0][0][0]!=cell[0][0][1]) ? cell[0][0][0]->faceave(2,1,FACEDATA_B)*dL1 : cell[0][0][0]->celldata[CELLDATA_B][2]*dL1;
            bf2 = (cell[0][0][0]!=cell[0][1][0]) ? cell[0][0][0]->faceave(1,1,FACEDATA_B)*dL2 : cell[0][0][0]->celldata[CELLDATA_B][1]*dL2;
            bf3 = (cell[0][1][0]!=cell[0][1][1]) ? cell[0][1][0]->faceave(2,1,FACEDATA_B)*dL3 : cell[0][1][0]->celldata[CELLDATA_B][2]*dL3;
            bf4 = (cell[0][1][1]!=cell[0][0][1]) ? cell[0][1][1]->faceave(1,0,FACEDATA_B)*dL4 : cell[0][1][1]->celldata[CELLDATA_B][1]*dL4;
            eneg = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        break;
    case 1 :  // y-edges
        if (cell[0][1][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][1][1] or cell[1][1][1]!=cell[0][1][1]) { // positive side edge exists?
            // length of line 1 =  r2*sin(theta)*dphi
            // length of line 3 = -r1*sin(theta)*dphi (dphi < 0)
            // length of line 2 = -dr (dr < 0)
            // length of line 4 =  dr
            real r1     = cell[0][1][0]->sph_centroid[0];
            real r2     = cell[1][1][0]->sph_centroid[0];
            real theta  = cell[0][1][0]->sph_centroid[1];
            real phi1   = cell[0][1][0]->sph_centroid[2];
            real phi2   = cell[0][1][1]->sph_centroid[2];
            real dr     = r2 - r1;
            real dphi   = phi2 - phi1;
            real dL1 = -dr;
            real dL2 =  r1*sin(theta)*dphi;
            real dL3 =  dr;
            real dL4 = -r2*sin(theta)*dphi;
            //real dS00 = cell[0][1][0]->sph_dS_centroid[1];
            //real dS01 = cell[0][1][1]->sph_dS_centroid[1];
            //real dS10 = cell[1][1][0]->sph_dS_centroid[1];
            //real dS11 = cell[1][1][1]->sph_dS_centroid[1];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS = 0.5*dr*(r2 + r1)*sin(theta)*dphi;
            if (dS<0.0) dS = -dS;
            posf=true;
            bf1 = (cell[0][1][0]!=cell[1][1][0]) ? cell[0][1][0]->faceave(0,1,FACEDATA_B)*dL1 : cell[0][1][0]->celldata[CELLDATA_B][0]*dL1;
            bf2 = (cell[0][1][0]!=cell[0][1][1]) ? cell[0][1][0]->faceave(2,1,FACEDATA_B)*dL2 : cell[0][1][0]->celldata[CELLDATA_B][2]*dL2;
            bf3 = (cell[0][1][1]!=cell[1][1][1]) ? cell[0][1][1]->faceave(0,1,FACEDATA_B)*dL3 : cell[0][1][1]->celldata[CELLDATA_B][0]*dL3;
            bf4 = (cell[1][1][1]!=cell[1][1][0]) ? cell[1][1][1]->faceave(2,0,FACEDATA_B)*dL4 : cell[1][1][1]->celldata[CELLDATA_B][2]*dL4;
            epos = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        if (cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][0][1] or cell[1][0][1]!=cell[0][0][1]) { // negative side edge exists?
            real r1     = cell[0][0][0]->sph_centroid[0];
            real r2     = cell[1][0][0]->sph_centroid[0];
            real theta  = cell[0][0][0]->sph_centroid[1];
            real phi1   = cell[0][0][0]->sph_centroid[2];
            real phi2   = cell[0][0][1]->sph_centroid[2];
            real dr     = r2 - r1;
            real dphi   = phi2 - phi1;
            real dL1 = -dr;
            real dL2 =  r1*sin(theta)*dphi;
            real dL3 =  dr;
            real dL4 = -r2*sin(theta)*dphi;
            //real dS00 = cell[0][0][0]->sph_dS_centroid[1];
            //real dS01 = cell[0][0][1]->sph_dS_centroid[1];
            //real dS10 = cell[1][0][0]->sph_dS_centroid[1];
            //real dS11 = cell[1][0][1]->sph_dS_centroid[1];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS = 0.5*dr*(r2 + r1)*sin(theta)*dphi;
            if (dS<0.0) dS = -dS;
            negf=true;
            bf1 = (cell[0][0][0]!=cell[1][0][0]) ? cell[0][0][0]->faceave(0,1,FACEDATA_B)*dL1 : cell[0][0][0]->celldata[CELLDATA_B][0]*dL1;
            bf2 = (cell[0][0][0]!=cell[0][0][1]) ? cell[0][0][0]->faceave(2,1,FACEDATA_B)*dL2 : cell[0][0][0]->celldata[CELLDATA_B][2]*dL2;
            bf3 = (cell[0][0][1]!=cell[1][0][1]) ? cell[0][0][1]->faceave(0,1,FACEDATA_B)*dL3 : cell[0][0][1]->celldata[CELLDATA_B][0]*dL3;
            bf4 = (cell[1][0][1]!=cell[1][0][0]) ? cell[1][0][1]->faceave(2,0,FACEDATA_B)*dL4 : cell[1][0][1]->celldata[CELLDATA_B][2]*dL4;
            eneg = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        break;
    case 2 :  // z-edges
        if (cell[0][0][1]!=cell[1][0][1] or cell[1][0][1]!=cell[1][1][1] or cell[1][1][1]!=cell[0][1][1]) { // positive side edge exists?
            // length of line01 =  dr
            // length of line12 =  r2*dtheta
            // length of line23 = -dr (dr < 0)
            // length of line30 = -r1*dtheta (dtheta < 0)
            real r1     = cell[0][0][1]->sph_centroid[0];
            real r2     = cell[1][0][1]->sph_centroid[0];
            real theta1 = cell[0][0][1]->sph_centroid[1];
            real theta2 = cell[0][1][1]->sph_centroid[1];
            real dr     = r2 - r1;
            real dtheta = theta2 - theta1;
            real dL1 = -r1*dtheta;
            real dL2 =  dr;
            real dL3 =  r2*dtheta;
            real dL4 = -dr;
            //real dS00 = cell[0][0][1]->sph_dS_centroid[2];
            //real dS01 = cell[0][1][1]->sph_dS_centroid[2];
            //real dS10 = cell[1][0][1]->sph_dS_centroid[2];
            //real dS11 = cell[1][1][1]->sph_dS_centroid[2];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
            if (dS<0.0) dS = -dS;
            posf=true;
            bf1 = (cell[0][0][1]!=cell[0][1][1]) ? cell[0][0][1]->faceave(1,1,FACEDATA_B)*dL1 : cell[0][0][1]->celldata[CELLDATA_B][1]*dL1;
            bf2 = (cell[0][0][1]!=cell[1][0][1]) ? cell[0][0][1]->faceave(0,1,FACEDATA_B)*dL2 : cell[0][0][1]->celldata[CELLDATA_B][0]*dL2;
            bf3 = (cell[1][1][1]!=cell[1][0][1]) ? cell[1][0][1]->faceave(1,1,FACEDATA_B)*dL3 : cell[1][0][1]->celldata[CELLDATA_B][1]*dL3;
            bf4 = (cell[1][1][1]!=cell[0][1][1]) ? cell[1][1][1]->faceave(0,0,FACEDATA_B)*dL4 : cell[1][1][1]->celldata[CELLDATA_B][0]*dL4;
            epos = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        if (cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][1][0] or cell[1][1][0]!=cell[0][1][0]) { // negative side edge exists?
            // length of line01 =  dr
            // length of line12 =  r2*dtheta
            // length of line23 = -dr (dr < 0)
            // length of line30 = -r1*dtheta (dtheta < 0)
            real r1     = cell[0][0][0]->sph_centroid[0];
            real r2     = cell[1][0][0]->sph_centroid[0];
            real theta1 = cell[0][0][0]->sph_centroid[1];
            real theta2 = cell[0][1][0]->sph_centroid[1];
            real dr     = r2 - r1;
            real dtheta = theta2 - theta1;
            real dL1 = -r1*dtheta;
            real dL2 =  dr;
            real dL3 =  r2*dtheta;
            real dL4 = -dr;
            //real dS00 = cell[0][0][0]->sph_dS_centroid[2];
            //real dS01 = cell[0][1][0]->sph_dS_centroid[2];
            //real dS10 = cell[1][0][0]->sph_dS_centroid[2];
            //real dS11 = cell[1][1][0]->sph_dS_centroid[2];
            //real dS = (abs(dS00) + abs(dS01) + abs(dS10) + abs(dS11))/4;
            real dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
            if (dS<0.0) dS = -dS;
            negf=true;
            bf1 = (cell[0][0][0]!=cell[0][1][0]) ? cell[0][0][0]->faceave(1,1,FACEDATA_B)*dL1 : cell[0][0][0]->celldata[CELLDATA_B][1]*dL1;
            bf2 = (cell[0][0][0]!=cell[1][0][0]) ? cell[0][0][0]->faceave(0,1,FACEDATA_B)*dL2 : cell[0][0][0]->celldata[CELLDATA_B][0]*dL2;
            bf3 = (cell[1][1][0]!=cell[1][0][0]) ? cell[1][0][0]->faceave(1,1,FACEDATA_B)*dL3 : cell[1][0][0]->celldata[CELLDATA_B][1]*dL3;
            bf4 = (cell[1][1][0]!=cell[0][1][0]) ? cell[1][1][0]->faceave(0,0,FACEDATA_B)*dL4 : cell[1][1][0]->celldata[CELLDATA_B][0]*dL4;
            eneg = invmu0*(bf1+bf2+bf3+bf4)/dS;
        }
        break;
    default :
        cerr << "Error (FATAL): jcomponent called with invalid direction." <<endl;
        exit(-1);
    }
    intp = (posf and negf) ? 0.5 : 1.0;
    return intp*(epos+eneg);
}

//! (SPHERICAL) New Spherical version of "jcomponent". Based on Cartesian approach.
real Tgrid::Tnode::sph_jcomponent_2(int edir)
{
    // n->cell[0][0][0] = cells[c];
    // n->cell[1][0][0] = cells[flatindex(i+1,j,  k  )];
    // n->cell[0][1][0] = cells[flatindex(i,  j+1,k  )];
    // n->cell[1][1][0] = cells[flatindex(i+1,j+1,k  )];
    // n->cell[0][0][1] = cells[flatindex(i,  j,  k+1)];
    // n->cell[1][0][1] = cells[flatindex(i+1,j,  k+1)];
    // n->cell[0][1][1] = cells[flatindex(i,  j+1,k+1)];
    // n->cell[1][1][1] = cells[flatindex(i+1,j+1,k+1)];
    real bf1, bf2, bf3, bf4; // B on faces for line integral.
    real invmu0 = 1./Params::mu_0;
    real epos=0.;
    real eneg=0.;
    real intp;
    bool posf=false;
    bool negf=false;
    real dS = 0.0;
////////////////////////////////////////////////////////////////////
    gridreal rc_000[3] = {cell[0][0][0]->sph_centroid[0], cell[0][0][0]->sph_centroid[1], cell[0][0][0]->sph_centroid[2]};
    gridreal rc_100[3] = {cell[1][0][0]->sph_centroid[0], cell[1][0][0]->sph_centroid[1], cell[1][0][0]->sph_centroid[2]};
    gridreal rc_010[3] = {cell[0][1][0]->sph_centroid[0], cell[0][1][0]->sph_centroid[1], cell[0][1][0]->sph_centroid[2]};
    gridreal rc_110[3] = {cell[1][1][0]->sph_centroid[0], cell[1][1][0]->sph_centroid[1], cell[1][1][0]->sph_centroid[2]};
    gridreal rc_001[3] = {cell[0][0][1]->sph_centroid[0], cell[0][0][1]->sph_centroid[1], cell[0][0][1]->sph_centroid[2]};
    gridreal rc_101[3] = {cell[1][0][1]->sph_centroid[0], cell[1][0][1]->sph_centroid[1], cell[1][0][1]->sph_centroid[2]};
    gridreal rc_011[3] = {cell[0][1][1]->sph_centroid[0], cell[0][1][1]->sph_centroid[1], cell[0][1][1]->sph_centroid[2]};
    gridreal rc_111[3] = {cell[1][1][1]->sph_centroid[0], cell[1][1][1]->sph_centroid[1], cell[1][1][1]->sph_centroid[2]};
    gridreal Ac_000[3] = {cell[0][0][0]->celldata[CELLDATA_B][0], cell[0][0][0]->celldata[CELLDATA_B][1], cell[0][0][0]->celldata[CELLDATA_B][2]};
    gridreal Ac_100[3] = {cell[1][0][0]->celldata[CELLDATA_B][0], cell[1][0][0]->celldata[CELLDATA_B][1], cell[1][0][0]->celldata[CELLDATA_B][2]};
    gridreal Ac_010[3] = {cell[0][1][0]->celldata[CELLDATA_B][0], cell[0][1][0]->celldata[CELLDATA_B][1], cell[0][1][0]->celldata[CELLDATA_B][2]};
    gridreal Ac_110[3] = {cell[1][1][0]->celldata[CELLDATA_B][0], cell[1][1][0]->celldata[CELLDATA_B][1], cell[1][1][0]->celldata[CELLDATA_B][2]};
    gridreal Ac_001[3] = {cell[0][0][1]->celldata[CELLDATA_B][0], cell[0][0][1]->celldata[CELLDATA_B][1], cell[0][0][1]->celldata[CELLDATA_B][2]};
    gridreal Ac_101[3] = {cell[1][0][1]->celldata[CELLDATA_B][0], cell[1][0][1]->celldata[CELLDATA_B][1], cell[1][0][1]->celldata[CELLDATA_B][2]};
    gridreal Ac_011[3] = {cell[0][1][1]->celldata[CELLDATA_B][0], cell[0][1][1]->celldata[CELLDATA_B][1], cell[0][1][1]->celldata[CELLDATA_B][2]};
    gridreal Ac_111[3] = {cell[1][1][1]->celldata[CELLDATA_B][0], cell[1][1][1]->celldata[CELLDATA_B][1], cell[1][1][1]->celldata[CELLDATA_B][2]};
    sph_transf_S2C_A1(rc_000, Ac_000);
    sph_transf_S2C_A1(rc_001, Ac_001);
    sph_transf_S2C_A1(rc_010, Ac_010);
    sph_transf_S2C_A1(rc_011, Ac_011);
    sph_transf_S2C_A1(rc_100, Ac_100);
    sph_transf_S2C_A1(rc_101, Ac_101);
    sph_transf_S2C_A1(rc_110, Ac_110);
    sph_transf_S2C_A1(rc_111, Ac_111);
////////////////////////////////////////////////////////////////////
//switch (edir)
//       {
//        case 0 :  // x-edges
//        if (cell[1][0][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][0][1]
//            or cell[1][0][1]!=cell[1][1][1]) // positive side edge exists (it might not because of refinement)?
    if (edir == 0) {
        real Ax0 = Ac_100[0];
        real Ay0 = Ac_100[1];
        real Az0 = Ac_100[2];
        real Ax1 = Ac_110[0];
        real Ay1 = Ac_110[1];
        real Az1 = Ac_110[2];
        real Ax2 = Ac_111[0];
        real Ay2 = Ac_111[1];
        real Az2 = Ac_111[2];
        real Ax3 = Ac_101[0];
        real Ay3 = Ac_101[1];
        real Az3 = Ac_101[2];
        //  [3] --- line23 ----  [2] ------ phi2      length of line01 =  r*dtheta
        //   |                    |                   length of line23 = -r*dtheta (dtheta < 0)
        //   |                    |                   length of line12 =  r*sin(theta2)*dphi
        // line30    Br(o)      line12       \/       length of line30 = -r*sin(theta1)*dphi (dphi < 0)
        //   |                    |
        //   |                    |
        //  [0] --- line01 ----  [1] ------ phi1
        //   |                    |
        //   |                    |
        // theta1      <        theta2
        real r      = cell[1][0][0]->sph_centroid[0];
        real theta1 = cell[1][0][0]->sph_centroid[1];
        real theta2 = cell[1][1][0]->sph_centroid[1];
        real phi1   = cell[1][0][0]->sph_centroid[2];
        real phi2   = cell[1][0][1]->sph_centroid[2];
        real dtheta = theta2 - theta1;
        real dphi   = phi2 - phi1;
        real L01 =  0.5*(Ax0 + Ax1)*r*cos(phi1)*(sin(theta2)-sin(theta1)) + 0.5*(Ay0 + Ay1)*r*sin(phi1)*(sin(theta2)-sin(theta1))
                    +0.5*(Az0 + Az1)*r*(cos(theta2)-cos(theta1));
        real L23 =  0.5*(Ax2 + Ax3)*r*cos(phi2)*(sin(theta1)-sin(theta2)) + 0.5*(Ay2 + Ay3)*r*sin(phi2)*(sin(theta1)-sin(theta2))
                    +0.5*(Az2 + Az3)*r*(cos(theta1)-cos(theta2));
        real L12 =  0.5*(Ax1 + Ax2)*r*sin(theta2)*(cos(phi2)-cos(phi1)) + 0.5*(Ay1 + Ay2)*r*sin(theta2)*(sin(phi2)-sin(phi1));
        real L30 =  0.5*(Ax3 + Ax0)*r*sin(theta1)*(cos(phi1)-cos(phi2)) + 0.5*(Ay3 + Ay0)*r*sin(theta1)*(sin(phi1)-sin(phi2));
        real dS = r*r*sin(theta1)*dtheta*dphi;
        if (dS<0.0) dS = -dS;
        posf=true;
        epos = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //if (cell[0][0][0]!=cell[0][1][0] or cell[0][1][0]!=cell[0][0][1]
    //    or cell[0][0][1]!=cell[0][1][1]) // negative side edge exists?
    if (edir == 0) {
        real Ax0 = Ac_000[0];
        real Ay0 = Ac_000[1];
        real Az0 = Ac_000[2];
        real Ax1 = Ac_010[0];
        real Ay1 = Ac_010[1];
        real Az1 = Ac_010[2];
        real Ax2 = Ac_011[0];
        real Ay2 = Ac_011[1];
        real Az2 = Ac_011[2];
        real Ax3 = Ac_001[0];
        real Ay3 = Ac_001[1];
        real Az3 = Ac_001[2];
        //  [3] --- line23 ----  [2] ------ phi2      length of line01 =  r*dtheta
        //   |                    |                   length of line23 = -r*dtheta (dtheta < 0)
        //   |                    |                   length of line12 =  r*sin(theta2)*dphi
        // line30    Br(o)      line12       \/       length of line30 = -r*sin(theta1)*dphi (dphi < 0)
        //   |                    |
        //   |                    |
        //  [0] --- line01 ----  [1] ------ phi1
        //   |                    |
        //   |                    |
        // theta1      <        theta2
        real r      = cell[0][0][0]->sph_centroid[0];
        real theta1 = cell[0][0][0]->sph_centroid[1];
        real theta2 = cell[0][1][0]->sph_centroid[1];
        real phi1   = cell[0][0][0]->sph_centroid[2];
        real phi2   = cell[0][0][1]->sph_centroid[2];
        real dtheta = theta2 - theta1;
        real dphi   = phi2 - phi1;
        real L01 =  0.5*(Ax0 + Ax1)*r*cos(phi1)*(sin(theta2)-sin(theta1)) + 0.5*(Ay0 + Ay1)*r*sin(phi1)*(sin(theta2)-sin(theta1))
                    +0.5*(Az0 + Az1)*r*(cos(theta2)-cos(theta1));
        real L23 =  0.5*(Ax2 + Ax3)*r*cos(phi2)*(sin(theta1)-sin(theta2)) + 0.5*(Ay2 + Ay3)*r*sin(phi2)*(sin(theta1)-sin(theta2))
                    +0.5*(Az2 + Az3)*r*(cos(theta1)-cos(theta2));
        real L12 =  0.5*(Ax1 + Ax2)*r*sin(theta2)*(cos(phi2)-cos(phi1)) + 0.5*(Ay1 + Ay2)*r*sin(theta2)*(sin(phi2)-sin(phi1));
        real L30 =  0.5*(Ax3 + Ax0)*r*sin(theta1)*(cos(phi1)-cos(phi2)) + 0.5*(Ay3 + Ay0)*r*sin(theta1)*(sin(phi1)-sin(phi2));
        real dS = r*r*sin(theta1)*dtheta*dphi;
        if (dS<0.0) dS = -dS;
        negf=true;
        eneg = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //break;
    //case 1 :  // y-edges
    //if (cell[0][1][0]!=cell[1][1][0] or cell[1][1][0]!=cell[1][1][1]
    //    or cell[1][1][1]!=cell[0][1][1]) // positive side edge exists?
    if (edir == 1) {
        real Ax0 = Ac_010[0];
        real Ay0 = Ac_010[1];
        real Az0 = Ac_010[2];
        real Ax1 = Ac_011[0];
        real Ay1 = Ac_011[1];
        real Az1 = Ac_011[2];
        real Ax2 = Ac_111[0];
        real Ay2 = Ac_111[1];
        real Az2 = Ac_111[2];
        real Ax3 = Ac_110[0];
        real Ay3 = Ac_110[1];
        real Az3 = Ac_110[2];
        //    [2]________
        //     | line12  \________[1] ----- phi2     length of line01 =  r1*sin(theta)*dphi
        //     |                   |                 length of line23 = -r2*sin(theta)*dphi (dphi < 0)
        //     |                   |                 length of line12 =  dr
        //   line23   Btheta(o)  line01              length of line30 = -dr (dr < 0)
        //     |                   |
        //     |         _________ |
        //     | _______/ line30  [0] ----- phi1
        //    [3]
        //     |                   |
        //    r2                   r1
        real r1     = cell[0][1][0]->sph_centroid[0];
        real r2     = cell[1][1][0]->sph_centroid[0];
        real theta  = cell[0][1][0]->sph_centroid[1];
        real phi1   = cell[0][1][0]->sph_centroid[2];
        real phi2   = cell[0][1][1]->sph_centroid[2];
        real dr     = r2 - r1;
        real dphi   = phi2 - phi1;
        real Er0 = Ax0*sin(theta)*cos(phi1) + Ay0*sin(theta)*sin(phi1) + Az0*cos(theta);
        real Er1 = Ax1*sin(theta)*cos(phi2) + Ay1*sin(theta)*sin(phi2) + Az1*cos(theta);
        real Er2 = Ax2*sin(theta)*cos(phi2) + Ay2*sin(theta)*sin(phi2) + Az2*cos(theta);
        real Er3 = Ax3*sin(theta)*cos(phi1) + Ay3*sin(theta)*sin(phi1) + Az3*cos(theta);
        real L01 =  0.5*(Ax0 + Ax1)*r1*sin(theta)*(cos(phi2)-cos(phi1)) + 0.5*(Ay0 + Ay1)*r1*sin(theta)*(sin(phi2)-sin(phi1));
        real L23 =  0.5*(Ax2 + Ax3)*r2*sin(theta)*(cos(phi1)-cos(phi2)) + 0.5*(Ay2 + Ay3)*r2*sin(theta)*(sin(phi1)-sin(phi2));
        real L12 =  0.5*(Er1 + Er2)*dr;
        real L30 = -0.5*(Er3 + Er0)*dr;
        real dS = 0.5*dr*(r2 + r1)*sin(theta)*dphi;
        if (dS<0.0) dS = -dS;
        posf=true;
        epos = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //if (cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][0][1]
    //    or cell[1][0][1]!=cell[0][0][1]) // negative side edge exists?
    if (edir == 1) {
        real Ax0 = Ac_000[0];
        real Ay0 = Ac_000[1];
        real Az0 = Ac_000[2];
        real Ax1 = Ac_001[0];
        real Ay1 = Ac_001[1];
        real Az1 = Ac_001[2];
        real Ax2 = Ac_101[0];
        real Ay2 = Ac_101[1];
        real Az2 = Ac_101[2];
        real Ax3 = Ac_100[0];
        real Ay3 = Ac_100[1];
        real Az3 = Ac_100[2];
        //    [2]________
        //     | line12  \________[1] ----- phi2     length of line01 =  r1*sin(theta)*dphi
        //     |                   |                 length of line23 = -r2*sin(theta)*dphi (dphi < 0)
        //     |                   |                 length of line12 =  dr
        //   line23   Btheta(o)  line01              length of line30 = -dr (dr < 0)
        //     |                   |
        //     |         _________ |
        //     | _______/ line30  [0] ----- phi1
        //    [3]
        //     |                   |
        //    r2                   r1
        real r1     = cell[0][0][0]->sph_centroid[0];
        real r2     = cell[1][0][0]->sph_centroid[0];
        real theta  = cell[0][0][0]->sph_centroid[1];
        real phi1   = cell[0][0][0]->sph_centroid[2];
        real phi2   = cell[0][0][1]->sph_centroid[2];
        real dr     = r2 - r1;
        real dphi   = phi2 - phi1;
        real Er0 = Ax0*sin(theta)*cos(phi1) + Ay0*sin(theta)*sin(phi1) + Az0*cos(theta);
        real Er1 = Ax1*sin(theta)*cos(phi2) + Ay1*sin(theta)*sin(phi2) + Az1*cos(theta);
        real Er2 = Ax2*sin(theta)*cos(phi2) + Ay2*sin(theta)*sin(phi2) + Az2*cos(theta);
        real Er3 = Ax3*sin(theta)*cos(phi1) + Ay3*sin(theta)*sin(phi1) + Az3*cos(theta);
        real L01 =  0.5*(Ax0 + Ax1)*r1*sin(theta)*(cos(phi2)-cos(phi1)) + 0.5*(Ay0 + Ay1)*r1*sin(theta)*(sin(phi2)-sin(phi1));
        real L23 =  0.5*(Ax2 + Ax3)*r2*sin(theta)*(cos(phi1)-cos(phi2)) + 0.5*(Ay2 + Ay3)*r2*sin(theta)*(sin(phi1)-sin(phi2));
        real L12 =  0.5*(Er1 + Er2)*dr;
        real L30 = -0.5*(Er3 + Er0)*dr;
        real dS = 0.5*dr*(r2 + r1)*sin(theta)*dphi;
        if (dS<0.0) dS = -dS;
        negf=true;
        eneg = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //break;
    //case 2 :  // z-edges
    //if (cell[0][0][1]!=cell[1][0][1] or cell[1][0][1]!=cell[1][1][1]
    //    or cell[1][1][1]!=cell[0][1][1]) // positive side edge exists?
    if (edir == 2) {
        real Ax0 = Ac_001[0];
        real Ay0 = Ac_001[1];
        real Az0 = Ac_001[2];
        real Ax1 = Ac_101[0];
        real Ay1 = Ac_101[1];
        real Az1 = Ac_101[2];
        real Ax2 = Ac_111[0];
        real Ay2 = Ac_111[1];
        real Az2 = Ac_111[2];
        real Ax3 = Ac_011[0];
        real Ay3 = Ac_011[1];
        real Az3 = Ac_011[2];
        /*          [0]- line30 -[3]       ----- r1    length of line01 =  dr
                    /             \                    length of line12 =  r2*dtheta
                   /               \                   length of line23 = -dr (dr < 0)
                  /                 \                  length of line30 = -r1*dtheta (dtheta < 0)
               line01    Bphi(o)   line23
                /                     \
               /                       \
              /                         \
            [1]-------- line12 --------[2] ----- r2
             |                          |
           theta1                     theta2
        */
        real r1     = cell[0][0][1]->sph_centroid[0];
        real r2     = cell[1][0][1]->sph_centroid[0];
        real theta1 = cell[0][0][1]->sph_centroid[1];
        real theta2 = cell[0][1][1]->sph_centroid[1];
        real phi    = cell[0][0][1]->sph_centroid[2];
        real dr     = r2 - r1;
        real dtheta = theta2 - theta1;
        real Er0 = Ax0*sin(theta1)*cos(phi) + Ay0*sin(theta1)*sin(phi) + Az0*cos(theta1);
        real Er1 = Ax1*sin(theta1)*cos(phi) + Ay1*sin(theta1)*sin(phi) + Az1*cos(theta1);
        real Er2 = Ax2*sin(theta2)*cos(phi) + Ay2*sin(theta2)*sin(phi) + Az2*cos(theta2);
        real Er3 = Ax3*sin(theta2)*cos(phi) + Ay3*sin(theta2)*sin(phi) + Az3*cos(theta2);
        real L01 =  0.5*(Er0 + Er1)*dr;
        real L23 = -0.5*(Er2 + Er3)*dr;
        real L12 =  0.5*(Ax1 + Ax2)*r2*cos(phi)*(sin(theta2)-sin(theta1)) + 0.5*(Ay1 + Ay2)*r2*sin(phi)*(sin(theta2)-sin(theta1))
                    +0.5*(Az1 + Az2)*r2*(cos(theta2)-cos(theta1));
        real L30 =  0.5*(Ax3 + Ax0)*r1*cos(phi)*(sin(theta1)-sin(theta2)) + 0.5*(Ay3 + Ay0)*r1*sin(phi)*(sin(theta1)-sin(theta2))
                    +0.5*(Az3 + Az0)*r1*(cos(theta1)-cos(theta2));
        real dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
        if (dS<0.0) dS = -dS; // Just in case ...
        posf=true;
        epos = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //if (cell[0][0][0]!=cell[1][0][0] or cell[1][0][0]!=cell[1][1][0]
    //    or cell[1][1][0]!=cell[0][1][0]) // negative side edge exists?
    if (edir == 2) {
        real Ax0 = Ac_000[0];
        real Ay0 = Ac_000[1];
        real Az0 = Ac_000[2];
        real Ax1 = Ac_100[0];
        real Ay1 = Ac_100[1];
        real Az1 = Ac_100[2];
        real Ax2 = Ac_110[0];
        real Ay2 = Ac_110[1];
        real Az2 = Ac_110[2];
        real Ax3 = Ac_010[0];
        real Ay3 = Ac_010[1];
        real Az3 = Ac_010[2];
        /*          [0]- line30 -[3]       ----- r1    length of line01 =  dr
                    /             \                    length of line12 =  r2*dtheta
                   /               \                   length of line23 = -dr (dr < 0)
                  /                 \                  length of line30 = -r1*dtheta (dtheta < 0)
               line01    Bphi(o)   line23
                /                     \
               /                       \
              /                         \
            [1]-------- line12 --------[2] ----- r2
             |                          |
           theta1                     theta2
        */
        real r1     = cell[0][0][0]->sph_centroid[0];
        real r2     = cell[1][0][0]->sph_centroid[0];
        real theta1 = cell[0][0][0]->sph_centroid[1];
        real theta2 = cell[0][1][0]->sph_centroid[1];
        real phi    = cell[0][0][0]->sph_centroid[2];
        real dr     = r2 - r1;
        real dtheta = theta2 - theta1;
        real Er0 = Ax0*sin(theta1)*cos(phi) + Ay0*sin(theta1)*sin(phi) + Az0*cos(theta1);
        real Er1 = Ax1*sin(theta1)*cos(phi) + Ay1*sin(theta1)*sin(phi) + Az1*cos(theta1);
        real Er2 = Ax2*sin(theta2)*cos(phi) + Ay2*sin(theta2)*sin(phi) + Az2*cos(theta2);
        real Er3 = Ax3*sin(theta2)*cos(phi) + Ay3*sin(theta2)*sin(phi) + Az3*cos(theta2);
        real L01 =  0.5*(Er0 + Er1)*dr;
        real L23 = -0.5*(Er2 + Er3)*dr;
        real L12 =  0.5*(Ax1 + Ax2)*r2*cos(phi)*(sin(theta2)-sin(theta1)) + 0.5*(Ay1 + Ay2)*r2*sin(phi)*(sin(theta2)-sin(theta1))
                    +0.5*(Az1 + Az2)*r2*(cos(theta2)-cos(theta1));
        real L30 =  0.5*(Ax3 + Ax0)*r1*cos(phi)*(sin(theta1)-sin(theta2)) + 0.5*(Ay3 + Ay0)*r1*sin(phi)*(sin(theta1)-sin(theta2))
                    +0.5*(Az3 + Az0)*r1*(cos(theta1)-cos(theta2));
        real dS  = 0.5*(r2*r2 - r1*r1)*dtheta;
        if (dS<0.0) dS = -dS; // Just in case ...
        negf=true;
        eneg = invmu0*(L01 + L12 + L23 + L30)/dS;
    }
    //break;
    //default :
    //    cerr << "Error (FATAL): jcomponent called with invalid direction." <<endl; exit(-1);
    //}
    intp = (posf and negf) ? 0.5 : 1.0;
    return intp*(epos+eneg);
}

//! (SPHERICAL) Spherical version of "calc_j"
void Tgrid::Tnode::sph_calc_j()
{
    nodedata[NODEDATA_J][0] = sph_jcomponent(0);
    nodedata[NODEDATA_J][1] = sph_jcomponent(1);
    nodedata[NODEDATA_J][2] = sph_jcomponent(2);
}

//! (SPHERICAL) New Spherical version of "calc_j". Based on Cartesian approach.
void Tgrid::Tnode::sph_calc_j_2()
{
    nodedata[NODEDATA_J][0] = sph_jcomponent_2(0);
    nodedata[NODEDATA_J][1] = sph_jcomponent_2(1);
    nodedata[NODEDATA_J][2] = sph_jcomponent_2(2);
}

//! (SPHERICAL) Spherical version of "calc_node_j_recursive"
void Tgrid::Tcell::sph_calc_node_j_recursive()
{
    int f,f2,ch,d=0;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_node_j_recursive();
    else {
        // NODE LOOP
        if (isrefined_face(d,1))
            for(f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_j();
        else
            face[d][1]->node[2]->sph_calc_j();
        for (d=1; d<3; d++)
            if (isrefined_face(d,1))
                for (f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_j();
    }
}

//! (SPHERICAL) New Spherical version of "calc_node_j_recursive". Based on Cartesian approach.
void Tgrid::Tcell::sph_calc_node_j_2_recursive()
{
    int f,f2,ch,d=0;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_node_j_2_recursive();
    else {
        // NODE LOOP
        if (isrefined_face(d,1))
            for(f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_j_2();
        else
            face[d][1]->node[2]->sph_calc_j_2();
        for (d=1; d<3; d++)
            if (isrefined_face(d,1))
                for (f=0; f<4; f++) for(f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_j_2();
    }
}

//! (SPHERICAL) Spherical version of "calc_node_j"
void Tgrid::sph_calc_node_j()
{
    int i,j,k;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++)
                cells[flatindex(i,j,k)]->sph_calc_node_j_recursive();
}

//! (SPHERICAL) New Spherical version of "calc_node_j". Based on Cartesian approach.
void Tgrid::sph_calc_node_j_2()
{
    int i,j,k;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++)
                cells[flatindex(i,j,k)]->sph_calc_node_j_2_recursive();
}


//! (SPHERICAL) Spherical version of "CN1_ne": interpolate ne to node for Ue calculation
void Tgrid::Tnode::sph_CN1_ne()
{
    real ne=0.;
    register Tcell *c;
    gridreal weightsum = 0.0;
    gridreal invweightsum;
    for (int a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        register const gridreal weight = c->sph_dV;
        ne+=weight*c->rho_q;
        weightsum+= weight;
    }
    invweightsum = (weightsum==0.) ? 0. : 1./weightsum;
    nodedata[NODEDATA_ne][0]=nodedata[NODEDATA_ne][1]=nodedata[NODEDATA_ne][2]=ne*invweightsum;
}

//! (SPHERICAL) Spherical version of "CN1_ne" + CN boundary calculations.
void Tgrid::Tnode::sph_CNb1_ne(int sph_CN_boundary_flag)
{
    int a;
    gridreal result = 0.0;
    int d;
    register Tcell *c;
    for (a=0; a<8; a++) {
        c = cell[0][0][a];
        if (c == 0) continue;
        // 1 internal region
        if (sph_CN_boundary_flag == 1)  {
            result += 1*c->rho_q;
        }
        // 6 faces
        if (sph_CN_boundary_flag == 6)  {
            result += 2*c->rho_q;
        }
        // 12 edges
        if (sph_CN_boundary_flag == 12) {
            result += 4*c->rho_q;
        }
        // 8 nodes
        if (sph_CN_boundary_flag == 8)  {
            result += 8*c->rho_q;
        }
    }
    nodedata[NODEDATA_ne][0] = nodedata[NODEDATA_ne][1] = nodedata[NODEDATA_ne][2] = result/8;
}

//! (SPHERICAL) Spherical version of "CN_ne_recursive"
void Tgrid::Tcell::sph_CN_ne_recursive()
{
    int ch,d,f,f2;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->sph_CN_ne_recursive();
    else {
        // NODE LOOP
        d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1))
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_ne();
        else
            face[d][1]->node[2]->sph_CN1_ne(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        for (d=1; d<3; d++) if (isrefined_face(d,1))
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_CN1_ne();
    }
}

//! (SPHERICAL) Spherical version of "CN_recursive"
void Tgrid::Tcell::sph_CNb_ne_recursive(int sph_CN_boundary_flag)
{
// NODE LOOP
// Idea: To pass through all nodes, pass through all right-pointing faces,
// and all upper-right corner points thereof
// IF YOU MODIFY THIS, MODIFY ALSO Tcell::calc_node_E_recursive below!!
    int d = 0;          // choose to pass in X-direction
    face[d][1]->node[2]->sph_CNb1_ne(sph_CN_boundary_flag);
}

//! (SPHERICAL) Spherical version of "CN_ne"
void Tgrid::sph_CN_ne()
{
    int i,j,k;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++)
                cells[flatindex(i,j,k)]->sph_CN_ne_recursive();
}

//! (SPHERICAL) Spherical version of "CN_ne" + CN boundary calcultion
void Tgrid::sph_CNb_ne()
{
    int i,j,k,c;
    // interpolation for the internal nodes which touch 8 internal cells
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_CNb_ne_recursive(1);
            }
    // interpolation for the internal side nodes which touch 4 internal cells
    // 6 faces
    i = 0;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    i = nx-2;
    for (j=1; j<ny-2; j++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    j = 0;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    j = ny-2;
    for (i=1; i<nx-2; i++) for (k=1; k<nz-2; k++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    k = 0;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    k = nz-2;
    for (i=1; i<nx-2; i++) for (j=1; j<ny-2; j++) {
            c = flatindex(i,j,k);
            cells[c]->sph_CNb_ne_recursive(6);
        }
    // interpolation for the internal side nodes which touch 2 internal cells
    // 12 different edges
    i = 0;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = 0;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = 0;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = 0;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = nx-2;
    j = 0;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = nx-2;
    j = ny-2;
    for (k=1; k<nz-2; k++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = nx-2;
    k = 0;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    i = nx-2;
    k = nz-2;
    for (j=1; j<ny-2; j++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    j = 0;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    j = 0;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    j = ny-2;
    k = 0;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    j = ny-2;
    k = nz-2;
    for (i=1; i<nx-2; i++) {
        c = flatindex(i,j,k);
        cells[c]->sph_CNb_ne_recursive(12);
    }
    // interpolation for the internal side nodes which touch 1 internal cell
    // 8 nodes
    i = 0;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = 0;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = 0;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = 0;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = nx-2;
    j = 0;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = nx-2;
    j = 0;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = nx-2;
    j = ny-2;
    k = 0;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
    i = nx-2;
    j = ny-2;
    k = nz-2;
    c = flatindex(i,j,k);
    cells[c]->sph_CNb_ne_recursive(8);
}

//! (SPHERICAL) Spherical version of "calc_ue1": Calculate Ue directly in node.
void Tgrid::Tnode::sph_calc_ue1()
{
    real rhoq =  nodedata[NODEDATA_ne][0];
    real invrhoq = (rhoq > 0) ? 1./rhoq : 0.;
    if (r0 < Params::R_zeroFields) {
        for(int d=0; d<3; d++) nodedata[NODEDATA_UE][d] =0.;
        return;
    }
    // Ji comes in Cartesian coordinates and J in spherical, so we need to transform J to Cartesian
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    gridreal J[3] = {nodedata[NODEDATA_J][0], nodedata[NODEDATA_J][1], nodedata[NODEDATA_J][2]};
    gridreal ue[3];
    sph_transf_S2C_A1(r, J);
    real ue2=0.;
    for(int d=0; d<3; d++) {
#ifndef IGNORE_ELECTRIC_FIELD_HALL_TERM
        ue[d] = (nodedata[NODEDATA_Ji][d] - J[d])*invrhoq;
#else
        ue[d] = nodedata[NODEDATA_Ji][d]*invrhoq;
#endif
        ue2+= sqr(ue[d]);
    }
    // Check maximum electron fluid velocity (CONSTRAINT)
    if (Params::Ue_max > 0 && ue2 > Params::Ue_max2) {
        const real norm = Params::Ue_max/sqrt(ue2);
        for(int d=0; d<3; d++) ue[d] *= norm;
#ifndef NO_DIAGNOSTICS
        Tgrid::fieldCounter.cutRateUe += 1.0;
#endif
    }
    sph_transf_C2S_A1(r, ue);
    for (int d=0; d<3; d++) nodedata[NODEDATA_UE][d] = ue[d];
}

//! (SPHERICAL) Spherical version of "calc_node_ue_recursive"
void Tgrid::Tcell::sph_calc_node_ue_recursive()
{
    int ch,d,f,f2;
    if (haschildren)
        for(ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_node_ue_recursive();
    else {
        // NODE LOOP
        d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1))
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_ue1();
        else
            face[d][1]->node[2]->sph_calc_ue1(); // assume nodes are numbered 0=(x,y), 1=(x+dx,y), 2=(x+dx,y+dy), 3=(x,y+dy)
        for (d=1; d<3; d++) if (isrefined_face(d,1))
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_ue1();
    }
}

//! (SPHERICAL) Spherical version of "calc_node_ue"
void Tgrid::sph_calc_node_ue()
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_calc_node_ue_recursive();
            }
}

//! (SPHERICAL) Calculate node E (recursive)
void Tgrid::Tcell::sph_calc_node_E_app_recursive(void)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->sph_calc_node_E_app_recursive();
    } else {
        // NODE LOOP
        int d = 0;          // choose to pass in X-direction
        if (isrefined_face(d,1)) {
            int f,f2;
            for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_E1_app();
        } else {
            // without refinement we are interesting only this line
            face[d][1]->node[2]->sph_calc_E1_app();
        }
        for (d=1; d<3; d++) if (isrefined_face(d,1)) {
                int f,f2;
                for (f=0; f<4; f++) for (f2=0; f2<4; f2++) refintf[d][1]->face[f]->node[f2]->sph_calc_E1_app();
            }
    }
}

//! (SPHERICAL)  Copy facedata to celldata
void Tgrid::Tcell::FC_copy_recursive(TFaceDataSelect fs, TCellDataSelect cs)
{
    if (haschildren) {
        int ch;
        for (ch=0; ch<8; ch++) child[0][0][ch]->FC_recursive(fs,cs);
    } else {
        int d;
        for (d=0; d<3; d++)
            celldata[cs][d] = faceave(d,1,fs);
    }
}

//! (SPHERICAL) Copy node value to cell
void Tgrid::Tcell::NC_copy_recursive(TNodeDataSelect ns,TCellDataSelect cs)
{
    int d;
    int dir = 1; //right face; face 1
    for (d=0; d<3; d++) {
        celldata[cs][d] = face[0][dir]->node[2]->nodedata[ns][d];
    }
}

//! (SPHERICAL) Spherical copy node value to cell
void Tgrid::Tcell::sph_NC_copy_recursive(TNodeDataSelect ns,TCellDataSelect cs)
{
    int d;
    int dir = 1; //right face; face 1
    gridreal r_node[3];
    gridreal A_node[3];
    gridreal r_cell[3];
    gridreal A_cell[3] = {0, 0, 0};
    for (d=0; d<3; d++) r_node[d] = face[0][dir]->node[2]->sph_centroid[d];
    for (d=0; d<3; d++) r_cell[d] = sph_centroid[d];
    for (d=0; d<3; d++) A_node[d] = face[0][dir]->node[2]->nodedata[ns][d];
    sph_transf_S2C_A1(r_node, A_node);
    for (d=0; d<3; d++) A_cell[d] = A_node[d];
    sph_transf_C2S_A1(r_cell, A_cell);
    for (d=0; d<3; d++) celldata[cs][d] = A_cell[d];
}

//! (SPHERICAL) Calculate node E
void Tgrid::Tnode::sph_calc_E1_app(void)
{
    real t = Params::t;
    real V = abs(Params::V);
    real x_max = Params::box_xmax;
    real z_max = Params::box_zmax;
    real S = V*t;
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    //gridreal x = r[0]*sin(r[1])*cos(r[2]);
    gridreal z = r[0]*cos(r[1]);
    real E[3] = {0.0, 0.0, 0.0};
    //if (x >= x_max - S) // (front to back. Cartesian)
    //if (z >= x_max - S) // (front to back. Cartesian)
    //if (r[0] >= x_max - S) // (front to back. Spherical)
    //if (r[0] <= S) // (front to back. Spherical)
    if (z >= 0.0 && z >= x_max - S) {
        E[0] = -4.0e-4;
        E[1] =  0.0e-4;
        E[2] =  0.0e-4;
    } else {
        E[0] = 0.0;
        E[1] = 0.0;
        E[2] = 0.0;
    }
    sph_transf_C2S_A(r, E);
    nodedata[NODEDATA_E][0] = E[0];
    nodedata[NODEDATA_E][1] = E[1];
    nodedata[NODEDATA_E][2] = E[2];
}

//! (SPHERICAL) Set vector in the cell center
void Tgrid::Tcell::sph_Cell_BC_recursive(TCellDataSelect cs)
{
    real t = Params::t;
    real V = abs(Params::V);
    real x_min = Params::box_xmin;
    real x_max = Params::box_xmax;
    real S = V*t;
    real A[3] = {0.0, 0.0, 0.0};
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    gridreal x    = r[0]*sin(r[1])*cos(r[2]);
    gridreal z    = r[0]*cos(r[1]);
    // Radial propagation
    if (Params::sph_propagation_type == 0) {
        A[0] = Params::boundary_Bx;
        A[1] = Params::boundary_By;
        A[2] = Params::boundary_Bz;
        celldata[cs][0] = A[0];
        celldata[cs][1] = A[1];
        celldata[cs][2] = A[2];
    }
    // Flat front propagation
    if (Params::sph_propagation_type == 1) {
        //Propagation along x-axis
        if (Params::sph_propagation_dir == 0) {
            if (x >= 0.0 && x >= x_max - S) {
                A[0] = Params::boundary_Bx;
                A[1] = Params::boundary_By;
                A[2] = Params::boundary_Bz;
                sph_transf_C2S_A(r, A);
                celldata[cs][0] = A[0];
                celldata[cs][1] = A[1];
                celldata[cs][2] = A[2];
            }
        }
        //Propagation along z-axis
        if (Params::sph_propagation_dir == 1) {
            if (z >= 0.0 && z >= x_max - S) {
                A[0] = Params::boundary_Bx;
                A[1] = Params::boundary_By;
                A[2] = Params::boundary_Bz;
                sph_transf_C2S_A(r, A);
                celldata[cs][0] = A[0];
                celldata[cs][1] = A[1];
                celldata[cs][2] = A[2];
            }
        } else {
            return;
        }
    }
}

//! (SPHERICAL) Set vector in the cell center
void Tgrid::Tcell::C_2_recursive(TCellDataSelect cs)
{
    real t = Params::t;
    real V = abs(Params::V);
    real x_min = Params::box_xmin;
    real x_max = Params::box_xmax;
    real S = V*t;
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    //gridreal x    = r[0]*sin(r[1])*cos(r[2]);
    gridreal z    = r[0]*cos(r[1]);
    real A[3]     = {0.0, 0.0, 0.0};
    //if (x >= x_max - S) // (front to back. Cartesian)
    if (z >= x_max - S) // (front to back. Cartesian)
        //if (r[0] >= x_max - S) // (front to back. Spherical)
    {
        A[0] =  0.0e-9;
        A[1] =  1.0e-9;
        A[2] =  0.0e-9;
    } else {
        A[0] = 0.0;
        A[1] = 0.0;
        A[2] = 0.0;
    }
    sph_transf_C2S_A(r, A);
    celldata[cs][0] = A[0];
    celldata[cs][1] = A[1];
    celldata[cs][2] = A[2];
}

//! (SPHERICAL) Calculate cell Ue
void Tgrid::Tcell::sph_calc_ue_app_recursive(void)
{
    real t = Params::t;
    real V = abs(Params::V);
    real x_max = Params::box_xmax;
    real S = V*t;
    gridreal r[3] = {sph_centroid[0], sph_centroid[1], sph_centroid[2]};
    //gridreal x = r[0]*sin(r[1])*cos(r[2]);
    gridreal z = r[0]*cos(r[1]);
    real Ui[3] = {0.0, 0.0, 0.0};
    real ue[3] = {0.0, 0.0, 0.0};
    // Charge density in thel cell
    real chargedensity = rho_q;
    //if (r[0] >= x_max - S) // (front to back. Spherical)
    //if (x >= x_max - S) // (front to back. Cartesian)
    if (z >= x_max - S) { // (front to back. Cartesian)
        Ui[0] =  0.0;
        Ui[1] =  0.0;
        Ui[2] = -4e5;
    } else {
        Ui[0] = 0.0;
        Ui[1] = 0.0;
        Ui[2] = 0.0;
    }
    //for (int d=0; d<3; d++) ue[d] = Ui[d] - (celldata[CELLDATA_J][d])/rho_q;
    for (int d=0; d<3; d++) ue[d] = Ui[d];
    // Ue transformation from Cartesian to hybrid coordinates
    //sph_transf_C2S_V(r, Ue);
    //sph_transf_S2H_V(Ue);
    sph_transf_C2S_A(r, ue);
    for (int d=0; d<3; d++) celldata[CELLDATA_UE][d] = ue[d];
}

//! (SPHERICAL) Copy node value to cell
void Tgrid::NC_copy(TNodeDataSelect ns,TCellDataSelect cs)
{
    int i,j,k,c;
    //ForInterior(i,j,k) {
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->NC_copy_recursive(ns,cs);
            }
}

//! (SPHERICAL) Copy node value to cell. Spherical version.
void Tgrid::sph_NC_copy(TNodeDataSelect ns,TCellDataSelect cs)
{
    int i,j,k,c;
    //ForInterior(i,j,k) {
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_NC_copy_recursive(ns,cs);
            }
}

//! (SPHERICAL) Set vector in the cell center
void Tgrid::C_2(TCellDataSelect cs)
{
    int i,j,k,c;
    ForAll(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->C_2_recursive(cs);
    }
}

//! (SPHERICAL) Set vector in the cell center
void Tgrid::sph_Cell_BC(TCellDataSelect cs)
{
    int i,j,k,c;
    // Radial propagation
    if (Params::sph_propagation_type == 0) {
        i = 1;
        for (j=0; j<ny; j++) for (k=0; k<nz; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_Cell_BC_recursive(cs);
            }
    }
    // Flat front propagation
    if (Params::sph_propagation_type == 1) {
        i = nx-2;
        for (j=0; j<ny; j++) for (k=0; k<nz; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_Cell_BC_recursive(cs);
            }
    }
}

//! (SPHERICAL) Copy facedata to celldata
void Tgrid::FC_copy(TFaceDataSelect fs, TCellDataSelect cs)
{
    int i,j,k,c;
    ForInterior(i,j,k) {
        c = flatindex(i,j,k);
        cells[c]->FC_copy_recursive(fs,cs);
    }
}

//! (SPHERICAL) Copy celldata
void Tgrid::copy_celldata_temp(int cTo, int cFrom, TCellDataSelect cs)
{
    int d = 2;
    //for (d=0; d<3; d++)
    cells[cTo]->celldata[cs][d] = cells[cFrom]->celldata[cs][d];
}

//! (SPHERICAL) Calculate node E
void Tgrid::sph_calc_node_E_app(void)
{
    int i,j,k,c;
    for (i=0; i<nx-1; i++) for (j=0; j<ny-1; j++) for (k=0; k<nz-1; k++) {
                c = flatindex(i,j,k);
                cells[c]->sph_calc_node_E_app_recursive();
            }
}

#endif

