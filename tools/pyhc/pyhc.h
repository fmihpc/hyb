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

#ifndef PYHC_H 
#include "metagrid.H"
//#include "GridOpener.h"
#include "gridcache.H"


class PGrid   
{
    Tmetagrid *g;
    TGridCache gc;
    double Gamma, Invmu0, Mass;
    bool Pseudobackground;

    public:

    PGrid() {g = NULL;}

    void open(const char *fn)
    {
//        GridOpener go;
//        g = go.open(fn, Gamma, Invmu0, Mass, Pseudobackground);
        g = gc.open(fn, Gamma, Invmu0, Mass, Pseudobackground);
    }

    bool isok();
    double xmin0();
    double xmin1();
    double xmin2();
    double xmax0();
    double xmax1();
    double xmax2();
    double intpol(double x, double y, double z, const char *vnam);
    double zintpol(double x, double y, double z, const char *vnam);
};

#define PYHC_H
#endif
