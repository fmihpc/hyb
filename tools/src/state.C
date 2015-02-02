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

#ifdef __GNUC__
#  pragma implementation "state.H"
#endif

#include "state.H"
#include <math.h>

real params::gam = 5.0/3.0;
real params::gamminus1 = 5.0/3.0-1;
real params::inv_gamminus1 = 1/(5.0/3.0-1);
real params::mu0 = 4*M_PI*1e-7;
real params::invmu0 = 1/(4*M_PI*1e-7);
real params::sqrtmu0 = 0.001120998243279586;
real params::invsqrtmu0 = 892.0620580763856;
real params::c = 3e8;
real params::invc2 = 1.0/(3e8*3e8);
real params::eps0 = 1.0/(4*M_PI*1e-7*3e8*3e8);

