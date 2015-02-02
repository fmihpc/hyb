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

#ifndef MAPS_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "gengrid.H"

extern void RampMap2D(should_be_const real st[MAXDIM][VECLEN], real xy[MAXDIM][VECLEN], smallnat vlen);
extern void RampMap3D(should_be_const real st[MAXDIM][VECLEN], real xy[MAXDIM][VECLEN], smallnat vlen);

extern real thetacoeff, phicoeff, theta1, phi1;	// used by ArcMap3D
extern void ArcMap2D(should_be_const real st[MAXDIM][VECLEN], real xy[MAXDIM][VECLEN], smallnat vlen);
extern void ArcMap3D(should_be_const real st[MAXDIM][VECLEN], real xy[MAXDIM][VECLEN], smallnat vlen);

#define MAPS_H

#endif

