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

#ifndef GRIDCACHE_H

#ifdef __GNUC__
#  pragma interface
#endif

#include "metagrid.H"

#define MAX_VARS 1000
#define MAX_HEADER_LINES 1000

class TGridCache {
public:
	static bool verbose;
private:
	struct TGridInfo {
		Tmetagrid *gptr;
		TGridInfo *next;
		char *filename;
		int refcount;
		int timestamp;
		double gam,mu0,mass;
		bool isSpectraFile;
		bool pseudobackground;
		TGridInfo() {gptr=0; next=0; filename=0; refcount=0; pseudobackground=false;}
		~TGridInfo();
	};
	TGridInfo *list;
	static int current_timestamp;
	int Ngrids() const;
	void purge();		// purge away the oldest grid
	void DeleteEntry(TGridInfo *p, bool allow_disposing_of_last);
public:
	void init() {list=0;}
	TGridCache() {init();}
	Tmetagrid *open(const char *fn, double& gamma, double& invmu0, double& mass, bool& pseudobackground);
	void close(Tmetagrid *ptr);
	bool isSpectraFile() { return list->isSpectraFile; }
	~TGridCache();
};

extern TGridCache theGridCache;		/* GLOBAL VARIABLE, init() called in toglwin.C:TToglWindow::SetFilename */
extern real GridDimensionScaling;
extern void (*MapFunctionPtr)(const real[MAXDIM][VECLEN], real[MAXDIM][VECLEN], smallnat);

#define GRIDCACHE_H

#endif

