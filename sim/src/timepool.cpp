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

#if defined(__linux) || defined(__sgi) || defined(__DECCXX)
#define HAVE_GETRUSAGE 1
#else
#define HAVE_GETRUSAGE 0
#endif

#include <cstring>
#include <iostream>
#include <iomanip>
#include <ctime>
#if HAVE_GETRUSAGE
#include <sys/time.h>
#include <sys/resource.h>
#endif
#include "timepool.h"
#include "simulation.h"

using namespace std;

#if HAVE_GETRUSAGE
//! Get CPU seconds
inline double GetCPUSeconds()
{
    struct rusage ru;
    getrusage(RUSAGE_SELF,&ru);
    return ru.ru_utime.tv_sec + ru.ru_stime.tv_sec + 1e-6*(ru.ru_utime.tv_usec + ru.ru_stime.tv_usec);
}
#else
//! Get CPU seconds
inline double GetCPUSeconds()
{
    return clock()/double(CLOCKS_PER_SEC);
}
#endif

//! Get CPU seconds
double Ttimepool::cputime() const
{
    return GetCPUSeconds();
}

//! Constructor
Ttimepool::Ttimepool() : n(0), attached_index(-1)
{
    for (int i=0; i<MAX_TIMEPOOLS; i++) {
        t[i] = 0.0;
        str[i] = 0;
    }
    cputime0 = GetCPUSeconds();
    cputimeLast = cputime0;
}

//! Attach a timepool
void Ttimepool::attach(const char *s)
{
    int j;
    bool found = false;
    for (j=0; j<n; j++) if (!strcmp(s,str[j])) {
            found = true;
            break;
        }
    if (!found) {
        if (n >= MAX_TIMEPOOLS) {
            static bool FirstTime = true;
            if (FirstTime) {
                errorlog << "*** Ttimepool::attach(\"" << s << "\"): too many tags, ignored\n";
                FirstTime = false;
            }
            return;
        }
        j = n++;
        str[j] = strdup(s);
    }
    const double c = GetCPUSeconds();
    if (attached_index >= 0) {
        t[attached_index]+= c-cputimeLast;
    }
    cputimeLast = c;
    attached_index = j;
}

//! Destructor
Ttimepool::~Ttimepool()
{
    MSGFUNCTIONCALL("Ttimepool::~Ttimepool");
    int maxlen = 0;
    double ttot = 0;
    int i,L;
    for (i=0; i<n; i++) {
        L = strlen(str[i]);
        if (L > maxlen) maxlen = L;
        ttot+= t[i];
    }
    if (ttot == 0) ttot = 1;
    mainlog.precision(1);
    mainlog.flags(ios::fixed | ios::showpoint);
    mainlog << "|--------------- TIME USAGE ---------------|\n";
    for (i=0; i<n; i++) {
        mainlog << "| ";
        mainlog.setf(ios::left);
        mainlog.width(maxlen+2);
        mainlog << str[i] << ':';
        mainlog.unsetf(ios::left);
        mainlog.setf(ios::right);
        mainlog.width(10);
        mainlog << t[i] << " s (" << 100.0*t[i]/ttot << " %) ";
        mainlog.unsetf(ios::right);
        mainlog	<< "\n";
    }
    const double cpu = GetCPUSeconds()-cputime0;
    mainlog << "| ";
    mainlog.setf(ios::left);
    mainlog.width(maxlen+2);
    mainlog << "Other" << ':';
    mainlog.unsetf(ios::left);
    mainlog.setf(ios::right);
    mainlog.width(10);
    mainlog << cpu-ttot << " s (" << 100.0*(cpu-ttot)/ttot << " %) ";
    mainlog.unsetf(ios::right);
    mainlog << "\n";
    mainlog << "| ";
    mainlog.setf(ios::left);
    mainlog.width(maxlen+2);
    mainlog << "Total" << ':';
    mainlog.unsetf(ios::left);
    mainlog.setf(ios::right);
    mainlog.width(10);
    mainlog << cpu << " s (" << 100.0*cpu/ttot << " %) ";
    mainlog.unsetf(ios::right);
    mainlog << "\n";
    mainlog << "|------------------------------------------|\n";
    MSGFUNCTIONEND("Ttimepool::~Ttimepool");
}

