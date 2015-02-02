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
#include <iomanip>
#include <vector> 
#include <stdio.h> 
#include <string.h> 
#include <math.h> 
#include <stdlib.h> 
#include "Ptracer.h"
#include "Config.h"
#include "Track.h"
#include "PointReader.h" 
#include "TrackWriter.h" 
#include "variables.H" 
#include "fileheader.H" 

using namespace std;
extern double mp;


Ptracer::Ptracer(Config *conf)
{
    cfg = conf;

    /* read simulation box boundaries and other header info */
    Theader head(cfg->hcfiles[0].c_str());
    simBox[0]=head.getreal("xmin0");
    simBox[1]=head.getreal("xmax0");
    simBox[2]=head.getreal("xmin1");
    simBox[3]=head.getreal("xmax1");
    simBox[4]=head.getreal("xmin2");
    simBox[5]=head.getreal("xmax2");
    gridDims[0] = head.getint("n1")-2;   // -2 takes away ghost cells  
    gridDims[1] = head.getint("n2")-2;
    gridDims[2] = head.getint("n3")-2;
    cellWidth = head.getreal("dx");

    /* Open grids. */
    unsigned int i=0;
    for(; i<cfg->hcfiles.size();i++)
    {
        grids[i] = gridcache.open(cfg->hcfiles[i].c_str(),Gamma,Invmu0,Mass,Pseudobackground); 
        if(!grids[i]) cerr << "Error: cannot open HC file " << cfg->hcfiles[i] << endl, exit(1);
        if(cfg->verbose) cout << "Opened grid "<< i << " for " << cfg->hcfiles[i] << endl;
    }
    grids[i] = NULL;

    pgrid = NULL;
    if(cfg->pressuretermfile.size()>0)
    {
        pgrid = gridcache.open(cfg->pressuretermfile.c_str(),Gamma,Invmu0,Mass,Pseudobackground); 
        if(!pgrid) cerr << "Error: cannot open HC file " << cfg->pressuretermfile << endl, exit(1);
        if(cfg->verbose) cout << "Opened grid for pressure term file " << cfg->pressuretermfile << endl;
    }

    /* Setup boundaries */
    if( cfg->xmin > cfg->xmax or cfg->ymin > cfg->ymax or cfg->zmin > cfg->zmax)
        cerr << "Error: minimum is bigger than maximum. Check boundary limits in confguration file." << endl
             << cfg->xmin << " " << cfg->xmax << endl
             << cfg->ymin << " " << cfg->ymax << endl
             << cfg->zmin << " " << cfg->zmax << endl,
        exit(-1);

    /* are given values in the box? */
    if(cfg->xmin < simBox[0] or cfg->xmin > simBox[1]) cfg->xmin = simBox[0];
    if(cfg->ymin < simBox[2] or cfg->ymin > simBox[3]) cfg->ymin = simBox[2];
    if(cfg->zmin < simBox[4] or cfg->zmin > simBox[5]) cfg->zmin = simBox[4];
    if(cfg->xmax < simBox[0] or cfg->xmax > simBox[1]) cfg->xmax = simBox[1];
    if(cfg->ymax < simBox[2] or cfg->ymax > simBox[3]) cfg->ymax = simBox[3];
    if(cfg->zmax < simBox[4] or cfg->zmax > simBox[5]) cfg->zmax = simBox[5];
}


inline bool Ptracer::boundaryCheck(Tdimvec &X)
{
    double x=X(0); double y=X(1); double z=X(2);
    if( x <= cfg->xmin or x >= cfg->xmax or  
        y <= cfg->ymin or y >= cfg->ymax or  
        z <= cfg->zmin or z >= cfg->zmax)
        return false;
    if(cfg->planetary_boundary and (x*x+y*y+z*z) <= cfg->planetary_boundary)
        return false;
    if( isnan(x) or isnan(y) or isnan(z))
        return false;
    return true;
}


void Ptracer::readInitialPoints( const char* fname, bool iscfg)
{
    PointReader pr(cfg,simBox);
    pr.readPoints(fname, iscfg, plist, plist_bw);

    if(plist.empty() and plist_bw.empty())
        cerr << "Error: found no initial values." << endl, exit(-1);
}


//NOTE: dont use with a=c, a=b or b=c...
inline void crossProd(double *a, double *b, double *c) 
{
    c[0] =  a[1]*b[2] - b[1]*a[2];
    c[1] = -a[0]*b[2] + b[0]*a[2];
    c[2] =  a[0]*b[1] - b[0]*a[1];
}


inline bool Ptracer::addPressureTerm(Tdimvec &X, double *E)
{
    Tvariable var;

    if(!pgrid) // no pressure term file
        return true;

    if(!pgrid->intpol(X,cfg->intpolorder,true))
        return false; 
    var.select("Bx0",Gamma,Invmu0,Mass);
    E[0] += var.get(*pgrid,X);
    var.select("By0",Gamma,Invmu0,Mass);
    E[1] += var.get(*pgrid,X);
    var.select("Bz0",Gamma,Invmu0,Mass);
    E[2] += var.get(*pgrid,X);
    return true;
}


inline void  Ptracer::setDirection(bool bw, vector<Track>::iterator &it, 
                                    vector<Track>::iterator &it_end, double &ds)
{
    if(bw)
        it = plist_bw.begin(), it_end = plist_bw.end(), ds=-ds;
    else
        it = plist.begin(), it_end = plist.end();
}


inline double Ptracer::getOtherVar(int id, double *v, vector<Track>::iterator &it)
{
    switch(id)
    {
        case BUN_EK : 
            return  0.5*it->mass*(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
        case BUN_VX: 
            return v[0]; 
        case BUN_VY: 
            return v[1]; 
        case BUN_VZ: 
            return v[2]; 
        case PARID: 
            return  it->id;
        default :
            cerr << "Error: getOtherVar called with invalid variable id."<< endl, exit(-1);
    }
}


inline void Ptracer::saveVars( int n, Tdimvec &X, vector<Track>::iterator &it, Tmetagrid& g, double *v)
{
    Tvariable var;
    int k;

    if(cfg->endpts and n>0) n = 1; 
    
    it->x[n] = X(0); 
    it->y[n] = X(1); 
    it->z[n] = X(2);
    it->used = n+1;

    for(k=0; k < cfg->tvi_len; k++)
    {
        if(!var.select(cfg->tvars_intpol[k],Gamma,Invmu0,Mass))
            cerr << "Error: invalid trace variable: " << cfg->tvars_intpol[k] << endl
                 << "       Change the value of TRACEVARS to a correct value.", exit(-1);
        it->data[k][n] = var.get(g,X);
    }

    for(int l=0; l < cfg->tvo_len; l++)
        it->data[k+l][n] = getOtherVar(cfg->tvo_flags[l],v,it);
}


inline bool Ptracer::getUe(Tdimvec &X, double *ue)
{
    double dnq=0, apu_dnq=0, j[3], dnqu[3]={0.,0.,0.};
    Tvariable var;

    var.select("jx",Gamma,Invmu0,Mass);
    j[0] = var.get(*grids[0],X);
    var.select("jy",Gamma,Invmu0,Mass);
    j[1] = var.get(*grids[0],X);
    var.select("jz",Gamma,Invmu0,Mass);
    j[2] = var.get(*grids[0],X);

    for(int l=0; grids[l]; l++)
    {
        if(!grids[l]->intpol(X,cfg->intpolorder,true))
            return false; 

        var.select("n",Gamma,Invmu0,Mass);
        apu_dnq = var.get(*grids[l],X)*cfg->hcf_charge[l]/cfg->hcf_mass[l];
        dnq += apu_dnq;

        var.select("vx",Gamma,Invmu0,Mass);
        dnqu[0] += var.get(*grids[l],X)*apu_dnq;

        var.select("vy",Gamma,Invmu0,Mass);
        dnqu[1] += var.get(*grids[l],X)*apu_dnq;

        var.select("vz",Gamma,Invmu0,Mass);
        dnqu[2] += var.get(*grids[l],X)*apu_dnq;
    }
    ue[0] = 1./dnq * (dnqu[0]-j[0]);
    ue[1] = 1./dnq * (dnqu[1]-j[1]);
    ue[2] = 1./dnq * (dnqu[2]-j[2]);
    return true;
}


// Particle trace using E and Buneman scheme. (Used when E is not -UexB). 
// Main part copied (more or less) straight from Hybrid simulation code: HybridSimulation::PropagateV
// NOTE: Use this when the electron pressure is included in the simulation.
void Ptracer::trace_BxUe_gradpe(bool bw)
{
    Tvariable var;
    Tmetagrid& g = *grids[0];
    Tdimvec X(0,0,0);
    double ds = cfg->stepsize_p;
    vector<Track>::iterator it, it_end;

    setDirection(bw,it,it_end,ds);

    for(;it!=it_end;it++)    
    {
        double Efield[3], B[3], Ue[3]={0.,0.,0.}, t[3], s[3], dv[3], vm[3], v0[3], vp[3], t2, b2;
        double v[3] = {it->v0[0], it->v0[1], it->v0[2]};
        double qmideltT2 = 0.5*ds*it->charge/it->mass;
        X[0]=it->x[0], X[1]=it->y[0], X[2]=it->z[0];

        if(cfg->verbose) cout << "Tracing particle: mass=" << it->mass << " charge=" << it->charge << ".";

        int n;
        for(n=0;n < cfg->maxsteps-1; n++)
        {
            if(!g.intpol(X,cfg->intpolorder,true)) break;
            saveVars( n,X,it,g,v);
            if(!getUe(X,Ue)) break;

            var.select("Bx",Gamma,Invmu0,Mass);
            B[0] = var.get(g,X);
            var.select("By",Gamma,Invmu0,Mass);
            B[1] = var.get(g,X);
            var.select("Bz",Gamma,Invmu0,Mass);
            B[2] = var.get(g,X);

            crossProd(B,Ue,Efield);
            addPressureTerm(X,Efield);

            dv[0] = qmideltT2 * Efield[0];
            dv[1] = qmideltT2 * Efield[1];
            dv[2] = qmideltT2 * Efield[2];

            t[0] = qmideltT2 * B[0];
            t[1] = qmideltT2 * B[1];
            t[2] = qmideltT2 * B[2];

            t2 = t[0]*t[0] + t[1]*t[1] + t[2]*t[2];
            b2 = 2./(1.+t2);
            
            s[0] = b2 * t[0];
            s[1] = b2 * t[1];
            s[2] = b2 * t[2];

            vm[0] = v[0] + dv[0];
            vm[1] = v[1] + dv[1];
            vm[2] = v[2] + dv[2];

            crossProd(vm,t,v0);
            v0[0] += vm[0];
            v0[1] += vm[1];
            v0[2] += vm[2];
            
            crossProd(v0,s,vp);
            vp[0] += vm[0];
            vp[1] += vm[1];
            vp[2] += vm[2];

            v[0] = vp[0] + dv[0];
            v[1] = vp[1] + dv[1];
            v[2] = vp[2] + dv[2];

            X[0] = X(0)+ds*v[0];
            X[1] = X(1)+ds*v[1];
            X[2] = X(2)+ds*v[2];

            if(!boundaryCheck(X)) break;
        }
        if(cfg->verbose) 
            cerr << "  Stopped. Step = " << n << ", point: (" << it->x[it->used-1] 
                 << "," << it->y[it->used-1] << "," << it->z[it->used-1] << ")" <<endl;
    }
}


void Ptracer::trace_ExB_drift(bool bw)
{
    Tvariable var;
    Tmetagrid& g = *grids[0];
    Tdimvec X(0,0,0);
    double ds = cfg->stepsize_p;
    vector<Track>::iterator it, it_end;

    setDirection(bw,it,it_end,ds);

    for(;it!=it_end;it++)    
    {
        double E[3], B[3], Ue[3]={0.,0.,0.}; // t[3], s[3], dv[3], vm[3], v0[3], vp[3], t2, b2;
        double v[3] = {it->v0[0], it->v0[1], it->v0[2]};
        double B2 = 0.;
        X[0]=it->x[0], X[1]=it->y[0], X[2]=it->z[0];

        if(cfg->verbose) cout << "Tracing particle: mass=" << it->mass << " charge=" << it->charge << ".";

        int n;
        for(n=0;n < cfg->maxsteps-1; n++)
        {
            if(!g.intpol(X,cfg->intpolorder,true)) break;
            saveVars( n,X,it,g,v);
            if(!getUe(X,Ue)) break;

            var.select("Bx",Gamma,Invmu0,Mass);
            B[0] = var.get(g,X);
            var.select("By",Gamma,Invmu0,Mass);
            B[1] = var.get(g,X);
            var.select("Bz",Gamma,Invmu0,Mass);
            B[2] = var.get(g,X);

            crossProd(B,Ue,E);
            addPressureTerm(X,E);
            crossProd(E,B,v);

            B2 = B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
            v[0] /= B2;
            v[1] /= B2;
            v[2] /= B2;

            X[0] = X(0)+ds*v[0];
            X[1] = X(1)+ds*v[1];
            X[2] = X(2)+ds*v[2];

            if(!boundaryCheck(X)) break;
        }
        if(cfg->verbose) 
            cerr << "  Stopped. Step = " << n << ", point: (" << it->x[it->used-1] 
                 << "," << it->y[it->used-1] << "," << it->z[it->used-1] << ")" <<endl;
    }
}


// Particle trace using U_e and Buneman scheme. 
// U_e is calculated from several hc files...
void Ptracer::trace_BxUe(bool bw)
{
    Tvariable var;
    Tmetagrid& g = *grids[0];
    Tdimvec X(0,0,0);
    double ds = cfg->stepsize_p;
    vector<Track>::iterator it, it_end;

    setDirection(bw,it,it_end,ds);

    for(;it!=it_end;it++)    
    {
        double v1[3] = {it->v0[0], it->v0[1], it->v0[2]};
        double om[3],w[3],wXom[3],wXomXom[3],ue[3]={0.,0.,0.};
        double o2;
        double a = ds*it->charge/it->mass;
        X[0]=it->x[0], X[1]=it->y[0], X[2]=it->z[0];

        if(cfg->verbose) cout << "Tracing particle: mass=" << it->mass << " charge=" << it->charge << ".";

        int n;
        for(n=0;n < cfg->maxsteps-1; n++)
        {
            if(!g.intpol(X,cfg->intpolorder,true)) break;

            var.select("Bx",Gamma,Invmu0,Mass);
            om[0] = 0.5*a*var.get(*grids[0],X);
            var.select("By",Gamma,Invmu0,Mass);
            om[1] = 0.5*a*var.get(*grids[0],X);
            var.select("Bz",Gamma,Invmu0,Mass);
            om[2] = 0.5*a*var.get(*grids[0],X);

            if(!getUe(X,ue)) break;

            /* Save orbit data */
            saveVars(n,X,it,g,v1);

            /* Propagate velocity */
            o2 = om[0]*om[0] + om[1]*om[1] + om[2]*om[2];

            /* W = v_n-u_e */
            w[0] = v1[0]-ue[0];
            w[1] = v1[1]-ue[1];
            w[2] = v1[2]-ue[2];

            crossProd(w,om,wXom);
            crossProd(wXom,om,wXomXom);
            
            w[0] = wXom[0] + wXomXom[0];
            w[1] = wXom[1] + wXomXom[1];
            w[2] = wXom[2] + wXomXom[2];

            v1[0] = v1[0] + (2./(1.+o2)) * w[0];
            v1[1] = v1[1] + (2./(1.+o2)) * w[1];
            v1[2] = v1[2] + (2./(1.+o2)) * w[2];

            /* propagate particle. */ 
            X[0] = X(0)+ds*v1[0];
            X[1] = X(1)+ds*v1[1];
            X[2] = X(2)+ds*v1[2];

            if(!boundaryCheck(X)) break;
        }
        if(cfg->verbose) 
            cerr << "  Stopped. Step = " << n << ", point: (" << it->x[it->used-1] 
                 << "," << it->y[it->used-1] << "," << it->z[it->used-1] << ")" <<endl;
    }
}


void Ptracer::trace()
{
    bool bw = not (cfg->direction == "forward");
    bool fw = not (cfg->direction == "backward");

    if(cfg->bunemanversion == "U")
    {        
        if(bw) trace_BxUe(true);
        if(fw) trace_BxUe(false);
    }
    else if(cfg->bunemanversion == "E") //version E.
    {
        if(bw) trace_BxUe_gradpe(true);
        if(fw) trace_BxUe_gradpe(false);
    }
    else if(cfg->bunemanversion == "ExB") // ExB drift.
    {
        if(bw) trace_ExB_drift(true);
        if(fw) trace_ExB_drift(false);
    }
}


void Ptracer::write()
{
    TrackWriter tw(cfg,gridDims,simBox,cellWidth);

    for(unsigned int i=0; i < cfg->file_formats.size(); i++)
    {
        if(cfg->file_formats[i] == "vtk")
            tw.writeVTK(plist,plist_bw);

        if(cfg->file_formats[i] == "3D")
            tw.write3D(plist,plist_bw);

        if(cfg->file_formats[i] == "matlab")
            tw.writeMatlab(plist,plist_bw);
    }
    tw.writeCfg(plist,plist_bw);
}



