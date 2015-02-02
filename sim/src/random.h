/** This file is part of the HYB simulation platform.
 *
 *  Copyright 2014- Finnish Meteorological Institute
 *
 *  Implementation of the algorithm used in Numerical Recipe's ran2
 *  generator copied from GNU Scientific Library version 1.9. ran2
 *  implementation is copyrighted undel GPL v2 by:
 *
 *  Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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

#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include "definitions.h"

#define N_SHUFFLE 32

//! State of the random generator
typedef struct {
    unsigned long int x;
    unsigned long int y;
    unsigned long int n;
    unsigned long int shuffle[N_SHUFFLE];
}
Tstate;

//! Portable random generator
class Tportrand
{
private:
    Tstate state;
    unsigned long int next_int();
public:
    void init(unsigned long int seed);
    Tportrand() {
        init(1);
    }
    Tportrand(unsigned long int seed) {
        init(seed);
    }
    double next();
    bool save(std::ostream& os);
    bool save(const char *fn);
    bool load(std::istream& is);
    bool load(const char *fn);
};

extern Tportrand portrand;
#define uniformrnd() portrand.next()
extern fastreal gaussrnd();
extern fastreal derivgaussrnd(fastreal x0);

/** \brief Probabilistic real2int rounding
 * 
 * probround(x) (x >= 0) gives either floor(x) or ceil(x), with probability
 * depending on which one is closer. For example, probround(2.3) gives 2
 * with 70% probability and 3 with 30% probability.
 */
inline int probround(real x)
{
    if(x <= 0) {
        return 0;
    }
    const int f = int(floor(x));
    return (uniformrnd() < x-f) ? f+1 : f;
}

#endif

