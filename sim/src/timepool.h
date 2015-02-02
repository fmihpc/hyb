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

#ifndef TIMEPOOL_H
#define TIMEPOOL_H

//! Profiling
class Ttimepool
{
public:
    enum {MAX_TIMEPOOLS=30};  //!< Increase this if necessary (but probably 30 different time pools is quite enough)
private:
    double t[MAX_TIMEPOOLS];  //!< accumulated CPU time in each timepool
    char *str[MAX_TIMEPOOLS]; //!< name of each timepool
    int n;                    //!< number of timepools
    int attached_index;       //!< the index of currently attached timepool
    double cputime0;          //!< cputime() when this Ttimepool was constructed
    double cputimeLast;       //!< cputime() at previous attach() call, or cputime0 if no attach yet done
public:
    Ttimepool();
    double cputime() const;
    void attach(const char *s); //!< Call this with any string tag to start spending time in a new timepool
    void operator()(const char *s) {
        attach(s); //!< you can just say timepool("mytag") instead of timepool.attach("mytag")
    }
    ~Ttimepool();  //!< breakdown of time usage will be automatically output to mainlog when the destructor is called
};

#endif

