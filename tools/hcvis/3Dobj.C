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
#  pragma implementation "3Dobj.H"
#endif

#include "3Dobj.H"

T3DFieldLineBunchObject::T3DFieldLineBunchObject
    (int n1, const double point1[3],const double point2[3],
	 TVectorField vecfield, TTraceDirection dir, TFieldLineDistribution distr,
	 TLoopThresholdType thresholdtype, double threshold)
{
	n = n1;
	int d;
	for (d=0; d<3; d++) {
		r1[d] = point1[d];
		r2[d] = point2[d];
	}
	VectorField = vecfield;
	TraceDirection = dir;
	Distribution = distr;
	LoopThresholdType = thresholdtype;
	LoopThreshold = threshold;
}
