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

#ifndef TRACK 

#include <Config.h>
#include <vector>
#include <string>
#include <stdlib.h>
using namespace std;


class Track{
    private:
        friend class Ptracer;
        friend class Ftracer;
        friend class TrackWriter;
        double *x,*y,*z;                // orbit coordinates.
        double *data[ORBIT_DATA_MAX];   // orbit data.
        vector<string>  data_names;     // names for variables in data.
        int used;                       // how many used of size.
        int size;                       // maxlength.
        int id;                         // id number for particle.

        /* Used only by particle tracing. */
        double mass;        
        double charge;      
        double v0[3]; // initial velocity.
    
    public:
        Track();
        Track(double x0, double y0, double z0, int l, int Id);
        Track(double x0, double y0, double z0, double vx, double vy, double vz, double m, double c, int l, 
              char** tvars_intpol, int tvi_len, char** tvars_other, int tvo_len, int Id);
        Track(const Track &other);
        Track& operator=(const Track& other);
        ~Track();
};

#define TRACK 
#endif
