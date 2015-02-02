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
#  pragma implementation "gridcache.H"
#endif

#include "gridcache.H"
#include "maps.H"
#include "constants.H"
#include <stdio.h>

bool TGridCache::verbose = false;
int TGridCache::current_timestamp = 0;
TGridCache theGridCache;
void (*MapFunctionPtr)(const real[MAXDIM][VECLEN], real[MAXDIM][VECLEN], smallnat) = 0;
real GridDimensionScaling = 1.0;

TGridCache::TGridInfo::~TGridInfo()
{
	if (TGridCache::verbose) cout << "TGridCache::close: disposing \"" << filename << "\"\n" << flush;
	if (filename) free(filename);
	if (gptr) delete gptr;
	filename = 0;
	gptr = 0;
}

static void FindOutMappingType(const char *fn)
{
	Theader h(fn);
	if (h.exists("maptype") && h.gettype("maptype") == Theader::STRING) {
		char *maptype = h.getstr("maptype");
		if (!strcmp(maptype,"ramp2")) {
			MapFunctionPtr = RampMap2D;
			cout << "Map type ramp2 in use\n" << flush;
		} else if (!strcmp(maptype,"ramp3")) {
			MapFunctionPtr = RampMap3D;
			cout << "Map type ramp3 in use\n" << flush;
		} else if (!strcmp(maptype,"arc3")) {
			MapFunctionPtr = ArcMap3D;
			cout << "Map type arc3 in use\n" << flush;

			{
				const real pi = M_PI;
				const real minlat = h.getreal("minlat");
				const real maxlat = h.getreal("maxlat");
				const real minlon = h.getreal("minlon");
				const real maxlon = h.getreal("maxlon");
				/* thetacoeff, phicoeff, theta1, phi1 are globals (declared in maps.H, defined in maps.C) */
				theta1 = (pi/180)*(90 - maxlat);
				const real theta2 = (pi/180)*(90 - minlat);
				phi1 = (pi/180)*minlon;
				const real phi2 = (pi/180)*maxlon;
				const real smin = h.getreal("smin");
				const real smax = h.getreal("smax");
				const real VHcoeff = h.getreal("vhcoeff");	// we have VHcoeff more grid points in vertical than in horiz
				const real C = (smax-smin)/VHcoeff;
				thetacoeff = C/(theta2-theta1);
				phicoeff = C/(phi2-phi1);
			}

		} else if (!strcmp(maptype,"arc2")) {
			MapFunctionPtr = ArcMap2D;
			cout << "Map type arc2 in use\n" << flush;

			{
				const real pi = M_PI;
				const real minlat = h.getreal("minlat");
				const real maxlat = h.getreal("maxlat");
				theta1 = (pi/180)*(90 - maxlat);
				const real theta2 = (pi/180)*(90 - minlat);
				const real smin = h.getreal("smin");
				const real smax = h.getreal("smax");
				const real VHcoeff = h.getreal("vhcoeff");	// we have VHcoeff more grid points in vertical than in horiz
				const real C = (smax-smin)/VHcoeff;
				thetacoeff = C/(theta2-theta1);
			}

		} else {
			cerr << "%%% gridcache.C:FindOutMappingType(): Unknown maptype " << maptype << " in file " << fn << "\n";
		}
	}
}

static double FindGammaFromFile(const char *fn)
{
	FILE *fp;
	fp = fopen(fn,"r");
	int ret;
	double gamma;
	int cnt = 0;
	const int maxsize = 1024;
	char s[maxsize+1];
	do {
		fgets(s,maxsize,fp);
		ret = sscanf(s,"# gamma = %lf",&gamma);
		cnt++;
	} while (ret != 1 && cnt < MAX_HEADER_LINES);
	fclose(fp);
	if (cnt >= MAX_HEADER_LINES)
		return -9999.0;
	else
		return gamma;
}

static double FindMu0FromFile(const char *fn)
{
	FILE *fp;
	fp = fopen(fn,"r");
	int ret;
	double mu0;
	int cnt = 0;
	const int maxsize = 1024;
	char s[maxsize+1];
	do {
		fgets(s,maxsize,fp);
		ret = sscanf(s,"# mu0 = %lf",&mu0);
		cnt++;
	} while (ret != 1 && cnt < MAX_HEADER_LINES);
	fclose(fp);
	if (cnt >= MAX_HEADER_LINES)
		return -9999.0;
	else
		return mu0;
}

static double FindMassFromFile(const char *fn)
{
	FILE *fp;
	fp = fopen(fn,"r");
	int ret;
	double mass;
	int cnt = 0;
	const int maxsize = 1024;
	char s[maxsize+1];
	do {
		fgets(s,maxsize,fp);
		ret = sscanf(s,"# m = %lf",&mass);
		cnt++;
	} while (ret != 1 && cnt < MAX_HEADER_LINES);
	fclose(fp);
	if (cnt >= MAX_HEADER_LINES)
		return -9999.0;
	else
		return mass;
}

static bool FindSpectraFlagFromFile(const char *fn)
{
	FILE *fp;
	fp = fopen(fn,"r");
	int ret;
        double flag = -1;
	int cnt = 0;
	const int maxsize = 1024;
	char s[maxsize+1];
	do {
		fgets(s,maxsize,fp);
		ret = sscanf(s,"# spectra = %lf",&flag);
		cnt++;
	} while (ret != 1 && cnt < MAX_HEADER_LINES);
	fclose(fp);
	if (cnt >= MAX_HEADER_LINES || flag <= 0)
		return false;
	else
		return true;
}

static bool FindPseudoBackgroundFlagFromFile(const char *fn)
{
	FILE *fp;
	fp = fopen(fn,"r");
	int ret;
	int cnt = 0;
	const int maxsize = 1024;
	char s[maxsize+1];
	char truefalse[12] = "";
	do {
		fgets(s,maxsize,fp);
		if (!s || *s != '#') break;
		ret = sscanf(s,"# pseudobackground = %10s",truefalse);
		cnt++;
	} while (ret != 1 && cnt < MAX_HEADER_LINES);
	fclose(fp);
	if (cnt >= MAX_HEADER_LINES)
		return false;
	else
		return !strcmp(truefalse,"true");
}

void TGridCache::purge()
{
	// Purge the oldest grid from the cache.
	// Find the minimum timestap in the cache.
	TGridInfo *p;
	if (!list) return;		// do nothing if the cache is empty
	int min_stamp = list->timestamp;
	TGridInfo *purged_grid = list;
	for (p=list; p; p=p->next)
		if (p->timestamp < min_stamp) {
			min_stamp = p->timestamp;
			purged_grid = p;
		}
	// Purge the purged_grid
	if (verbose) cout << "TGridCache::purge: purging grid with timestamp=" << min_stamp << "\n";
	DeleteEntry(purged_grid,false);
}

int TGridCache::Ngrids() const
{
	TGridInfo *p;
	int n = 0;
	for (p=list; p; p=p->next) n++;
	return n;
}

Tmetagrid *TGridCache::open(const char *fn, double& gamma, double& invmu0, double& mass, bool& pseudobackground)
{
	TGridInfo *p;
	// Try to find the grid first
	for (p=list; p; p=p->next)
		if (!strcmp(fn,p->filename)) {
			// Found. Increase reference counter, and return the pointer
			p->refcount++;
			gamma = p->gam;
			invmu0 = 1.0/p->mu0;
		        mass = p->mass;
			return p->gptr;
		}
	// Not found. Purge the cache until at most 3 grids remain, then add a new entry in the list.
	while (Ngrids() > 3) purge();
	if (verbose) cout << "TGridCache::open: loading \"" << fn << "\"\n" << flush;
	Tmetagrid *const thegrid = new Tmetagrid(fn);
	if (!thegrid || !thegrid->good()) {
		cerr << "*** Could not open grid file \"" << fn << "\"\n";
		return 0;
	}
	thegrid->scale(GridDimensionScaling);		// global variable, default 1
	p = new TGridInfo;
	p->gptr = thegrid;
	p->next = list;
	p->filename = strdup(fn);
	p->refcount = 1;
	p->timestamp = ++current_timestamp;
	gamma = FindGammaFromFile(fn);
	if (gamma == -9999.0) {
		static bool FirstTime = true;
		if (FirstTime) {
			cerr << "note: gamma not stored in file, using 5/3\n";
			FirstTime = false;
		}
		gamma = 5.0/3.0;
	}
	p->gam = gamma;
	p->mu0 = FindMu0FromFile(fn);
	pseudobackground = p->pseudobackground = FindPseudoBackgroundFlagFromFile(fn);
	if (p->mu0 == -9999.0) {
		static bool FirstTime = true;
		if (FirstTime) {
			cerr << "note: mu0 not stored in file, using 4*pi*1e-7\n";
			FirstTime = false;
		}
		p->mu0 = 4*3.141592653589793*1e-7;
	}
        mass = FindMassFromFile(fn);
     	if (mass == -9999.0) {
		static bool FirstTime = true;
		if (FirstTime) {
			cerr << "note: particle mass not stored in file, using m = mp\n";
			FirstTime = false;
		}
	        mass = cnst::mp;
	}
        p->mass=mass;
        p->isSpectraFile = FindSpectraFlagFromFile(fn);
	if (verbose) cout << "TGridCache::open: using gamma=" << gamma << ", mu0=" << p->mu0 << ", m=" << mass << "\n" << flush;
	invmu0 = 1.0/p->mu0;
	list = p;
	static bool FirstTime = true;
	if (FirstTime) {
		FindOutMappingType(fn);
		FirstTime = false;
	}
	return p->gptr;
}

void TGridCache::DeleteEntry(TGridInfo *p, bool allow_disposing_of_last)
{
	if (p->refcount > 0) p->refcount--;	// do not allow negative reference counts
	if (p->refcount == 0 && (list->next || allow_disposing_of_last)) {
		// list->next is here so that we do not delete our last grid
		if (p == list) {
			list = p->next;
		} else {
			TGridInfo *prev;
			for (prev=list; prev; prev=prev->next)
				if (prev->next == p) break;
			if (prev->next != p) {
				cerr << "*** TGridCache::close: internal error\n";
				return;
			}
			prev->next = p->next;
		}
		delete p;
	}
}

void TGridCache::close(Tmetagrid *ptr)
{
	const bool allow_disposing_of_last = true /*false*/;
	TGridInfo *p, *foundptr=0;
	for (p=list; p; p=p->next)
		if (ptr == p->gptr) {
			foundptr = p;
			break;
		}
	if (!foundptr) {
		cerr << "*** TGridCache::close: could not find ptr in grid cache\n";
		return;
	}
	DeleteEntry(foundptr, allow_disposing_of_last);
}

TGridCache::~TGridCache() {
	TGridInfo *p;
	while (list) {
		p = list;
		list = p->next;
		delete p;
	}
}

