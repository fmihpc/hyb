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

#ifndef CACHE_H

#ifdef __GNUC__
#  pragma interface
#endif

// cache_size is in TCacheWord (=double) units
typedef double TCacheWord;

class Tcache {
private:
	TCacheWord *cache;
	int words_allocated;
	int cache_size;
	bool verbose_flag;
public:
	void init();
	Tcache() {init();}
	void global_reset() {words_allocated=0;}
	int bytes_free() const {return (cache_size - words_allocated)*sizeof(TCacheWord);}
	int size() const {return cache_size;}
	int size_in_bytes() const {return cache_size*sizeof(TCacheWord);}
	void verbose(bool flag=true) {verbose_flag=flag;}
	bool verboseflag() const {return verbose_flag;}
	TCacheWord *alloc(int sz, const char *msg=0);	// allocate sz bytes from the cache, use TCacheWord(=double) alignment
};

extern Tcache Dcache;		/* GLOBAL VARIABLE, init() called in shmemintf.C:startpes() */

#define CACHE_H

#endif
