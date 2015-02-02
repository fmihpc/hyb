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

#ifndef POINTREADER_H

#include <string>
#include "Track.h"
#include "Config.h"
#include <string.h>
#include <stdlib.h>


class PointReader{
        bool fromCfg;
        string fname;
        Config *cfg;
        double mass;
        double charge;
        double velocity[3];
        double length_scale;
        double boxLimits[6];   // xmin,xmax,ymin,ymax,..
        int track_len;
        int npts;

    public:
        PointReader( Config*, double *);
        void readPoints( const char* fn, bool iscfg, vector<Track> &plist, vector<Track> &plist_bw);
};

#define POINTREADER_H
#endif
