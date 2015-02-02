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

#ifndef FORBIDSPLITJOIN_H
#define FORBIDSPLITJOIN_H

#include <cmath>
#include <cstdio>
#include "definitions.h"
#include "atmosphere.h"

//! Forbid split&join profiles
class ForbidSplitAndJoinProfile : public ScalarField
{
public:
    ForbidSplitAndJoinProfile();
    ForbidSplitAndJoinProfile(std::string funcName,std::vector<real> args);
    ~ForbidSplitAndJoinProfile();
    real getValue(const gridreal r[]);
private:
    bool (ForbidSplitAndJoinProfile::*ptr)(const gridreal[]);
    bool defaultFunction(const gridreal r[]);
    // FORBID SPLITANDJOIN PROFILES
    bool forbidSplitAndJoinDefault(const gridreal r[]);
    bool forbidSplitAndJoinInsideSphere(const gridreal r[]);
    bool forbidSplitAndJoinOutsideSphere(const gridreal r[]);
    bool forbidSplitAndJoinInsideCuboid(const gridreal r[]);
    bool forbidSplitAndJoinOutsideCuboid(const gridreal r[]);
    bool forbidSplitAndJoinTitan(const gridreal r[]);
};

void initializeForbidSplitJoin();

#endif

