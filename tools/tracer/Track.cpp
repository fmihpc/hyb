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

#include "Track.h"
#include <iostream>
#include <vector>
#include <string>
#include <string.h>
using namespace std;


Track::Track(double x0, double y0, double z0, int l, int Id)
{
    size = l;
    x = new double[l];
    y = new double[l];
    z = new double[l];
                 x[0] = x0;
    y[0] = y0;
    z[0] = z0;

    for(int i=0;i<ORBIT_DATA_MAX;i++)
    {
        data[i] = new double[l];
        data[i][0] = 0.;
    }

    v0[0] = 0.;
    v0[1] = 0.;
    v0[2] = 0.;
    mass = charge = 0.;
    used = 1;
    id = Id;
    cerr << "Error: simple Track constructor is for testing" << endl;
    exit(-1);
}


Track::Track(double x0, double y0, double z0, double vx, double vy, double vz, double m, double c, int l,
         char** tvars_intpol, int tvi_len, char** tvars_other, int tvo_len, int Id)
{
    size = l;
    x = new double[l];
    y = new double[l];
    z = new double[l];
    x[0] = x0;
    y[0] = y0;
    z[0] = z0;

    for(int i=0;i<ORBIT_DATA_MAX;i++)
    {
        data[i] = new double[l];
        data[i][0] = 0.;
    }

    string help;
    for(int i=0;i<tvi_len;i++)
    {
        help = tvars_intpol[i];
        data_names.push_back(help);
    }
    for(int i=0;i<tvo_len;i++)
    {
        help = tvars_other[i];
        data_names.push_back(help);
    }

    v0[0] = vx;
    v0[1] = vy;
    v0[2] = vz;
    mass = m;
    charge = c;
    used = 1;
    id = Id;
}

// Copy constructor
Track::Track(const Track &other)
{
    size = other.size;
    used = other.used;
    x = new double[size];
    y = new double[size];
    z = new double[size];
    memcpy(x, other.x, used*sizeof(double));
    memcpy(y, other.y, used*sizeof(double));
    memcpy(z, other.z, used*sizeof(double));

    for(int i=0;i<ORBIT_DATA_MAX;i++)
    {
        data[i] = new double[size];
        memcpy(data[i], other.data[i], used*sizeof(double));
    }
    data_names = other.data_names;

//    vector<string>::iterator it = data_names.begin();
//    for(;it<data_names.end();it++)
//        cout << "debug: " << (*it) << " " << endl;

    v0[0] = other.v0[0];
    v0[1] = other.v0[1];
    v0[2] = other.v0[2];
    mass = other.mass;
    charge = other.charge;
    id = other.id;
}


Track& Track::operator=(const Track& other) 
{
    if(this != &other)
    {
        if(used > 0)
        {
            delete [] x;
            delete [] y;
            delete [] z;
            for(int i=0;i<ORBIT_DATA_MAX;i++)
                delete [] data[i];
        }

        size = other.size;
        used = other.used;
        x = new double[size];
        y = new double[size];
        z = new double[size];
        memcpy(x, other.x, used*sizeof(double));
        memcpy(y, other.y, used*sizeof(double));
        memcpy(z, other.z, used*sizeof(double));

        for(int i=0;i<ORBIT_DATA_MAX;i++)
        {
            data[i] = new double[size];
            memcpy(data[i], other.data[i], used*sizeof(double));
        }

        data_names = other.data_names;

        v0[0] = other.v0[0];
        v0[1] = other.v0[1];
        v0[2] = other.v0[2];
        mass = other.mass;
        charge = other.charge;
        id = other.id;
    }
    return *this;
}


Track::~Track()
{
    size = 0;
    used = 0;
    delete [] x;
    delete [] y;
    delete [] z;
    for(int i=0;i<ORBIT_DATA_MAX;i++)
        delete [] data[i];
    mass = charge = v0[0] = v0[1] = v0[2] = 0.;
}



