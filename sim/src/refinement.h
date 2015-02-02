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

#ifndef REFINEMENT_H
#define REFINEMENT_H

#include <cmath>
#include <cstdio>
#include "definitions.h"
#include "atmosphere.h"

/** \brief Hierarchical grid refinement profile
 *
 * input is the cell coordinates
 * output is the cellsize wanted [m]
 * Params::dx (the toplevel size can be used)
 * The condition for recursive refinement:
 * size(of the cell) > output  ->  refine the cell
 * output is the output of these functions for the center of the cell
 * The grid is automatically constructed such that the grid spacing
 * is everywhere at least as small as the value returned by this
 * function at the corresponding point, subject only to the condition
 * that the maximum refinement level (maxGridRefinementLevel) is not exceeded.
 */
class GridRefinementProfile : public ScalarField
{
public:
    GridRefinementProfile();
    GridRefinementProfile(std::string funcName,std::vector<real> args);
    ~GridRefinementProfile();
    real getValue(const gridreal r[3]);
private:
    gridreal (GridRefinementProfile::*ptr)(const gridreal[]);
    gridreal defaultFunction(const gridreal r[3]);
    // GRID REFINEMENT PROFILES
    gridreal refineSpherical(const gridreal r[]);
    gridreal refineCartesian(const gridreal r[]);
    gridreal refineSphericalAndCartesian(const gridreal r[]);
    gridreal refineTitan(const gridreal r[]);
};

void initializeGridRefinement();

#endif

