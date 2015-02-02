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

#include "pyhc.h"
#include "metagrid.H"
#include "variables.H"
#include <iostream>

using namespace std;


double PGrid::intpol(double x, double y, double z, const char *vnam)
{
    Tvariable var;
    Tdimvec Xc(x,y,z);
    double ret;

    g->intpol(Xc,1,true);
    var.select(vnam,Gamma,Invmu0,Mass);
    ret = var.get(*g,Xc);
    return ret;
}

double PGrid::zintpol(double x, double y, double z, const char *vnam)
{
    Tvariable var;
    Tdimvec Xc(x,y,z);
    double ret;

    g->intpol(Xc,0,true);
    var.select(vnam,Gamma,Invmu0,Mass);
    ret = var.get(*g,Xc);
    return ret;
}

bool PGrid::isok()
{
    if(!g)
        return false;
    return g->good();
}

double PGrid::xmin0()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmin[0];
}
double PGrid::xmin1()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmin[1];
}
double PGrid::xmin2()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmin[2];
}
double PGrid::xmax0()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmax[0];
}
double PGrid::xmax1()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmax[1];
}
double PGrid::xmax2()
{
    real xmin[3], xmax[3];
    g->getbox(xmin,xmax);
    return xmax[2];
}
