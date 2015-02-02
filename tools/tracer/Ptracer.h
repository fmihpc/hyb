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

#ifndef PTRACER_H

#include <vector>
#include <stdlib.h>
#include "Track.h"
#include "variables.H" 
#include "Config.h"
#include "gridcache.H" 
using namespace std;


class Ptracer{
    private:
        vector<Track> plist;       // particle orbits forward or stream(field)line forward.
        vector<Track> plist_bw;    // particle orbits backward or stream(field)line backward.
        Config *cfg;

        // grid info
        Tmetagrid *grids[MAX_HCF];
        Tmetagrid *pgrid;  // grid for pressure file
        TGridCache gridcache;
        double simBox[6];  // order: x min max, y min max, z min max.
        int gridDims[3];   // not including ghost cells.
        double cellWidth;

        // variables.C needs this
        double Gamma, Invmu0, Mass;
        bool Pseudobackground;

        inline bool boundaryCheck(Tdimvec &X);
        inline void setDirection( bool bw, vector<Track>::iterator &it, 
                                vector<Track>::iterator &it_end, double &ds);
        inline bool getUe(Tdimvec &X, double *ue);
        inline void saveVars( int n, Tdimvec &X, vector<Track>::iterator &it, Tmetagrid& g, double *v);
        inline bool addPressureTerm(Tdimvec &X, double *E);
        inline double getOtherVar(int id, double *v, vector<Track>::iterator &it);

    public:
        Ptracer(Config *);
        void readInitialPoints(const char* fname, bool iscfg);
        void write();
        void trace();
        void trace_BxUe_gradpe(bool bw); // tracer for simulation with pressure term.
        void trace_BxUe(bool bw);        // tracer E(lorentz) = -UexB.
        void trace_ExB_drift(bool bw);
};

#define PTRACER_H
#endif


