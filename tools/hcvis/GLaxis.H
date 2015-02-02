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

#ifndef GLAXIS_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <GL/gl.h>

extern void DrawAxis
    (float x1, float y1, float z1,
	 float x2, float y2, float z2,
	 float dx, float dy, float dz,
	 double u1, double u2,
	 bool include_label, GLuint font);

#define GLAXIS_H

#endif

