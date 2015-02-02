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

#include <fstream>
#include <cstring>
#include "random.h"
#include "simulation.h"

using namespace std;

Tportrand portrand;

const long int m1 = 2147483563, a1 = 40014, q1 = 53668, r1 = 12211;
const long int m2 = 2147483399, a2 = 40692, q2 = 52774, r2 = 3791;
#define N_DIV (1 + 2147483562/N_SHUFFLE)

//! Initialize random generator
void Tportrand::init(unsigned long int s)
{
    if (s == 0) {
        s = 1;
    }
    state.y = s;
    for(int i=0; i<8; i++) {
        long int h = s / q1;
        long int t = a1 * (s - h * q1) - h * r1;
        if (t < 0) {
            t += m1;
        }
        s = t;
    }
    for(int i=N_SHUFFLE-1; i>=0; i--) {
        long int h = s / q1;
        long int t = a1 * (s - h * q1) - h * r1;
        if (t < 0) {
            t += m1;
        }
        s = t;
        state.shuffle[i] = s;
    }
    state.x = s;
    state.n = s;
    return;
}

//! Next random integer
unsigned long int Tportrand::next_int()
{
    const unsigned long int x = state.x;
    const unsigned long int y = state.y;
    long int h1 = x / q1;
    long int t1 = a1 * (x - h1 * q1) - h1 * r1;
    long int h2 = y / q2;
    long int t2 = a2 * (y - h2 * q2) - h2 * r2;
    if (t1 < 0) {
        t1 += m1;
    }
    if (t2 < 0) {
        t2 += m2;
    }
    state.x = t1;
    state.y = t2;
    {
        unsigned long int j = state.n / N_DIV;
        long int delta = state.shuffle[j] - t2;
        if (delta < 1) {
            delta += m1 - 1;
        }
        state.n = delta;
        state.shuffle[j] = t1;
    }
    return state.n;
}

//! Next random double
double Tportrand::next()
{
    float x_max = 1 - 1.2e-7f;
    float x = next_int() / 2147483563.0f;
    if (x > x_max) {
        return x_max;
    }
    return x;
}

//! Save state of random generator in stream
bool Tportrand::save(ostream& os)
{
    os << "# State of a portrand generator\n";
    os << state.x << ' ' << state.y << ' ' << state.n << '\n';
    for (int i=0; i<N_SHUFFLE; i++) {
        os << state.shuffle[i] << '\n';
    }
    return os.good();
}

//! Load state of random generator from stream
bool Tportrand::load(istream& is)
{
    char buf[80];
    is.getline(buf,78);
    if (strcmp(buf,"# State of a portrand generator")) {
        ERRORMSG("file is not new portrand state file");
        doabort();
        return false;
    }
    is >> state.x >> state.y >> state.n;
    for (int i=0; i<N_SHUFFLE; i++) {
        is >> state.shuffle[i];
    }
    is.getline(buf,78);    // eat newline character
    return is.good();
}

//! Save state of random generator in file
bool Tportrand::save(const char *fn)
{
    ofstream os(fn);
    if (!os.good()) {
        return false;
    }
    save(os);
    return os.good();
}

//! Load state of random generator from file
bool Tportrand::load(const char *fn)
{
    ifstream is(fn);
    if (!is.good()) {
        return false;
    }
    load(is);
    return is.good();
}

/** \brief Gaussian randomness
 *
 *  Generate a Gaussian deviate with zero mean and unit
 *  standard deviation: f(x) = (1/(2*pi))*exp(-0.5*x^2), x real.
 *  Algorithm: Generate random pairs (x,y) from unit square
 *  -1 <= x <= 1, -1 <= y <= 1 until (x,y) is within
 *  the unit circle. Compute fac = sqrt(-2.0*log(r2)/r2),
 *  where r2 = x^2 + y^2. Then, x*fac and y*fac are two Gaussian
 *  random numbers.
 */
fastreal gaussrnd()
{
    static fastreal saved;
    static bool is_saved = false;
    fastreal x,y,r2,fac,result;
    if (is_saved) {
        result = saved;
        is_saved = false;
    } else {
        do {
            x = 2*uniformrnd() - 1;
            y = 2*uniformrnd() - 1;
            r2 = x*x + y*y;
        } while (r2 >= 1.0);
        // On average, this do loop is executed 4/pi = 1.27324 times
        fac = sqrt(-2.0*log(r2)/r2);
        result = x*fac;
        saved = y*fac;
        is_saved = true;
    }
    return result;
}

/** \brief Deriv Gaussian randomness
 *
 * Return a random number distributed according to
 * f(x) = c*max(0,x)*exp(-0.5*(x-x0)^2) where the normalization constant c
 * is chosen so that the integrate(f(x),x,-inf,inf)=1 (notice that f(x)=0 for x<=0).
 *
 * Method: F(x)=c*xm*exp(-0.5*(x-xm)^2-0.5*(x0-xm)^2), where xm=0.5*(x0+sqrt(x0^2+4)),
 * is a majorant, i.e. F(x) >= f(x) for all x>=0 and x0. The majorant is Gaussian
 * with unit standard deviation and mean equal to xm. (Note that xm is the abscissa
 * of the maximum of f(x), i.e. f'(xm)=0.) Generate random numbers x from
 * the majorant Gaussian and accept it with probability f(x)/F(x).
 * The area under the majorant curve F(x) is close to unity for x0>=0 so that
 * only a few trials are needed. For x0<0 it is asymptotically proportional
 * to (-x0) so that more and more trials are needed. Therefore, avoid calling
 * the function with x0 < -10.
 */
fastreal derivgaussrnd(fastreal x0)
{
    fastreal x,majorant,pdf;
    const fastreal invsqrt2 = 1.0/sqrt(2.0);
    const fastreal sqrt_halfpi = sqrt(0.5*pi);
    const fastreal xm = 0.5*(x0 + sqrt(sqr(x0) + 4.0));
    const fastreal c = 1.0/(exp(-0.5*x0*x0) + x0*sqrt_halfpi*erfc(-x0*invsqrt2));
    const fastreal d = sqr(x0-xm);
    const fastreal cxm = c*xm;
restart:
    x = xm + gaussrnd();
    if (x < 0) goto restart;
    majorant = cxm*exp(-0.5*(sqr(x-xm) + d));
    pdf = c*x*exp(-0.5*sqr(x-x0));
    if (uniformrnd()*majorant > pdf) goto restart;
    return x;
}

