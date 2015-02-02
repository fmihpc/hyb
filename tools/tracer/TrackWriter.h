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

#ifndef TRACKWRITER_H

#include <iostream>
#include "Track.h"
#include "Config.h"
#include <stdlib.h>


class TrackWriter{
    private:
        Config *cfg;
        double cellWidth;
        double simBox[6];  // order: x min max, y min max, z min max.
        int gridDims[3];   // not including ghost cells.
        int fileNum;       // if overwrite==false, then this will tell what number to put in file extension.

        void outFilename( string &fn, const char* ext);
        int currentFileNum();
        void removeLastTrace();
        void writeVTKtrace( ofstream &of, vector<Track> &plist, vector<Track> &plist_bw);
        void writeASCII( vector<Track> &plist, vector<Track> &plist_bw, string ext);

    public:
        TrackWriter(Config *c, int *gd, double *sb, double cw);
        void writeVTK( vector<Track> &plist, vector<Track> &plist_bw);
        void write3D( vector<Track> &plist, vector<Track> &plist_bw)
        { string s="3D"; writeASCII( plist, plist_bw, s);}
        void writeMatlab( vector<Track> &plist, vector<Track> &plist_bw)
        { string s="m"; writeASCII( plist, plist_bw, s);}
        void writeCfg(vector<Track> &plist, vector<Track> &plist_bw);
};

#define TRACKWRITER_H
#endif
