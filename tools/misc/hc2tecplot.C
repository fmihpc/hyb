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

/*
 * See function usage() below for what this program does.
 */

#include <stdio.h>
#include <iostream>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "metagrid.H"

#define c8 8

static char *inputfn = 0;
static char *title = 0;
const static char *zone = 0;
static int digits = 6;
static char *outputfn = 0;
static bool IncludeGhosts = false;
static bool verbose = false;
static bool usehash = true;
static TGridIndex Np = 0;			// Number of corner points
static TGridIndex Ne = 0;			// Number of (leaf and interior) elements
static int ncd = 0;
static TGridIndex *elemtab = 0;
static real gam = 5.0/3.0;
static real mu0 = 4*3.141592653589793*1e-7;
static char format[80];

#define elemtabref(e,c) elemtab[(e)*c8+(c)]

enum {VAR_RHO,
	  VAR_RHOUX,VAR_RHOUY,VAR_RHOUZ, VAR_RHOVX,VAR_RHOVY,VAR_RHOVZ,
	  VAR_UX,VAR_UY,VAR_UZ, VAR_VX,VAR_VY,VAR_VZ,
	  VAR_V, VAR_V2,
	  VAR_U, VAR_P,
	  VAR_BX,VAR_BY,VAR_BZ,
	  VAR_BX0,VAR_BY0,VAR_BZ0,
	  VAR_BX1,VAR_BY1,VAR_BZ1,
	  VAR_B,VAR_B0,VAR_B1, VAR_B2,VAR_B02,VAR_B12};
#define VAR_LAST (VAR_B12+1)

static bool enabled[VAR_LAST] = {false};
const static char *varnames[VAR_LAST];



// ------------------------ Txyzlist ---------------------------

class Txyzlist {
private:
	struct Txyzlistitem {
		real x,y,z;
		Txyzlistitem *next;
	};
	Txyzlistitem *ptr,*lastptr,*rower;
	TGridIndex nitems;
public:
	void init() {ptr=lastptr=rower=0; nitems=0;}
	Txyzlist() {init();}
	TGridIndex add(real x1, real y1, real z1);
	void reset() {rower=ptr;}
	bool isend() const {return rower==0;}
	void next() {rower=rower->next;}
	real X() const {return rower->x;}
	real Y() const {return rower->y;}
	real Z() const {return rower->z;}
	TGridIndex length() const {return nitems;}
	friend ostream& operator<<(ostream& o, const Txyzlist& L);
	~Txyzlist();
};

ostream& operator<<(ostream& o, const Txyzlist& L)
{
	Txyzlist::Txyzlistitem *p;
	o << "[";
	for (p = L.ptr; p; p=p->next) o << " (" << p->x << "," << p->y << "," << p->z << ")";
	o << "]";
	return o;
}

TGridIndex Txyzlist::add(real x1, real y1, real z1)
{
	const TGridIndex result = nitems;
	Txyzlistitem *n;
	n = new Txyzlistitem;
	nitems++;
	n->x = x1;
	n->y = y1;
	n->z = z1;
	n->next = 0;
	if (ptr == 0) {
		ptr = lastptr = rower = n;
	} else {
		lastptr->next = n;
		lastptr = n;
	}
	return result;
}

Txyzlist::~Txyzlist()
{
	Txyzlistitem *p;
	while (ptr) {
		p = ptr;
		ptr = p->next;
		delete p;
	}
}

static Txyzlist xyzlist;

// ---------------------- Txyzhash -----------------------------

class Txyzhash {
private:
	struct THashEntry {
		long xq,yq,zq;
		TGridIndex p;
		THashEntry *next;
	};
	TGridIndex N;			// length of tab
	TGridIndex nitems;		// number of entries stored
	real Dx,x0,y0,z0;		// parameters needed in Quantize
	int nbits;				// needed in HashFunction
	THashEntry **tab;
	TGridIndex HashFunction(long,long,long) const;
	void Quantize(real,real,real, long&,long&,long&) const;
public:
	Txyzhash(TGridIndex n, real Dx1, real x01, real y01, real z01);
	TGridIndex search(real x, real y, real z) const;		// Return -1 if not found, otherwise index which is >= 0
	void store(real x, real y, real z, TGridIndex p);		// Use store only if search failed
	void stats(ostream& o) const;
	~Txyzhash();
};

Txyzhash::Txyzhash(TGridIndex n, real Dx1, real x01, real y01, real z01)
{
	typedef THashEntry *THashEntryPtr;
	N = n;
	nitems = 0;
	Dx = Dx1;
	x0 = x01;
	y0 = y01;
	z0 = z01;
	nbits = 16;
	tab = new THashEntryPtr [N];
	memset(tab,0,sizeof(THashEntryPtr)*N);
}

TGridIndex Txyzhash::HashFunction(long x, long y, long z) const
{
	return (((unsigned long)x >> 5) + ((unsigned long)y) + ((unsigned long)z << 5)) % N;
}

void Txyzhash::Quantize(real x, real y, real z, long& xq, long& yq, long& zq) const
{
	xq = long((1<<nbits)*((x-x0)/Dx) + 0.5);
	yq = long((1<<nbits)*((y-y0)/Dx) + 0.5);
	zq = long((1<<nbits)*((z-z0)/Dx) + 0.5);
}

TGridIndex Txyzhash::search(real x, real y, real z) const
{
	long xq,yq,zq;
	Quantize(x,y,z, xq,yq,zq);
	const TGridIndex i = HashFunction(xq,yq,zq);
	THashEntry *ptr;
	for (ptr=tab[i]; ptr; ptr=ptr->next)
		if (ptr->xq == xq && ptr->yq == yq && ptr->zq == zq) return ptr->p;
	return -1;
}

void Txyzhash::store(real x, real y, real z, TGridIndex p)
{
	long xq,yq,zq;
	Quantize(x,y,z, xq,yq,zq);
	const TGridIndex i = HashFunction(xq,yq,zq);
	THashEntry *ptr;
	ptr = new THashEntry;
	ptr->xq = xq;
	ptr->yq = yq;
	ptr->zq = zq;
	ptr->p = p;
	ptr->next = tab[i];
	tab[i] = ptr;
	nitems++;
}

void Txyzhash::stats(ostream& o) const
{
	o << "Hash table length is " << N << " and it contains " << nitems << " entries.\n";
	TGridIndex i,maxlen=0,totlen=0,nnull=0;
	THashEntry *ptr;
	for (i=0; i<N; i++) {
		TGridIndex len = 0;
		for (ptr = tab[i]; ptr; ptr=ptr->next) len++;
		if (len > maxlen) maxlen = len;
		if (len == 0) nnull++;
		totlen+= len;
	}
	if (totlen != nitems) o << "*** totlen=" << totlen << "!=nitems !\n";
	const double avelen = nitems/double(N);
	o << "Average chain length is " << avelen << ", maximum length is " << maxlen << ".\n";
	o << "Number of null chains is " << nnull << " (" << 100.0*nnull/double(N) << " %)\n";
#if 0
	// Output statistics to hashstats.dat
	FILE *fp;
	fp = fopen("hashstats.dat","w");
	for (i=0; i<N; i++) {
		TGridIndex len = 0;
		for (ptr = tab[i]; ptr; ptr=ptr->next) len++;
		fprintf(fp,"%d\n",int(len));
	}
	fclose(fp);
#endif
}

Txyzhash::~Txyzhash()
{
	TGridIndex i;
	THashEntry *ptr,*p;
	for (i=0; i<N; i++) {
		ptr = tab[i];
		while (ptr) {
			p = ptr;
			ptr = p->next;
			delete p;
		}
	}
	delete [] tab;
}

static TGridIndex SelectHashTableLength(TGridIndex x)
/*
  Theory: We select primes which are close to powers of 2.
  In Mathematica:
  nextPrime[x_Integer]:= Module[{m},m=x; While[!PrimeQ[m],m++];m];
  Table[nextPrime[2^n],{n,2,32}] -->
  {5,11,17,37,67,131,257,521,1031,2053,4099,8209,16411,32771,65537,131101,
  262147,524309,1048583,2097169,4194319,8388617,16777259,33554467,67108879,
  134217757,268435459,536870923,1073741827,2147483659,4294967311}
 */
{
	static TGridIndex possibleLengths[] = {
		2053,4099,8209,16411,32771,65537,131101,262147,
		524309,1048583,2097169,4194319,8388617,16777259,33554467,67108879,134217757,
		268435459,536870923,1073741827,2147483659 /*,4294967311*/,0};
	int i;
	for (i=0; possibleLengths[i] > 0; i++) if (possibleLengths[i] > x) return possibleLengths[i];
	return possibleLengths[8];
}

// -------------------------------------------------------------

static void usage()
{
	clog << "usage: hc2tecplot [-t title] [-z zone] [-d digits] [-g] [-v] [-F]\n";
	clog << "                  [-l varlist] [-o outfile] inputfile.hc\n";
	clog << "Translation of a 3D HC file into Tecplot ASCII file. Options:\n";
	clog << "-t title   The TITLE appearing in the file. The default is composed\n";
	clog << "           of the input filename and the #t = .. setting in file header\n";
	clog << "-z zone    The ZONE string. Default is \"Whole grid\"\n";
	clog << "-d digits  Number of significant ASCII digits. Default is 6.\n";
	clog << "-g         Include ghosts cells also\n";
	clog << "-v         Verbose output\n";
	clog << "-F         Produce Full corner vertex lists, do not attempt to remove duplicates\n";
	clog << "           by hashing. Try -F if you have problems otherwise. With -F,\n";
	clog << "           the output file will be many times larger because it contains\n";
	clog << "           exactly 8 corners for every cell.\n";
	clog << "-o outfile Name of the output file. Default is to send to standard output.\n";
	clog << "-l varlist Comma-separated list of variables to write to Tecplot file.\n";
	clog << "           The default is rho,rhoUx,rhoUy,rhoUz,U,Bx,By,Bz\n";
	clog << "           The variables can be selected from (including aliases):\n";
	clog << "             rho rhoUx rhoUy rhoUz rhovx rhovy rhovz vx vy vz Ux Uy Uz\n";
	clog << "             v v2 U P\n";
	clog << "             Bx By Bz Bx0 Bxy0 Bz0 Bx1 By1 Bz1 B B0 B1 B2 B02 B12\n";
	clog << "           If the varlist starts with a comma, the variables are appended\n";
	clog << "           to the default varlist.\n";
	exit(0);
}

static void SetVar(const char *var)
{
	int a;
	if (!var) {cerr << "*** hc2tecplot warning: Syntax error in varlist\n"; return;}
	bool found = false;
	for (a=0; a<VAR_LAST; a++)
		if (!strcmp(varnames[a],var)) {
			found = true;
			break;
		}
	if (found)
		enabled[a] = true;
	else
		cerr << "*** hc2tecplot warning: Variable name '" << var << "' not recognized\n";
}

static void ParseVarlist(char *varlist)
{
	int a;
	bool reset_them = true;
	if (*varlist == ',') {varlist++; reset_them=false;}
	if (reset_them) {
		for (a=0; a<VAR_LAST; a++) enabled[a] = false;
	}
	const char *delim = ",";
	SetVar(strtok(varlist,delim));
	while (1) {
		char *const var = strtok(0,delim);
		if (!var) break;
		SetVar(var);
	} 
}

static void ParseArgs(int argc, char *argv[])
{
	int a;
	for (a=1; a<argc; a++)
		if (argv[a][0] == '-') {
			switch (argv[a][1]) {
			case 't': title = argv[++a]; break;
			case 'z': zone = argv[++a]; break;
			case 'o': outputfn = argv[++a]; break;
			case 'd': digits = atoi(argv[++a]); break;
			case 'g': IncludeGhosts = true; break;
			case 'v': verbose = true; break;
			case 'l': ParseVarlist(argv[++a]); break;
			case 'F': usehash = false; break;
			default: usage();
			}
		} else {
			if (a != argc-1) usage();
			inputfn = argv[a];
		}
	if (!inputfn) usage();
}

static void CheckInputFileAndTakeTime(const char *fn, double& t)
{
	FILE *fp;
	fp = fopen(fn,"r");
	if (!fp) {cerr << "*** hc2tecplot: Could not open input file \"" << fn << "\"\n"; exit(1);}
	fscanf(fp,"# t = %lf",&t);
	fclose(fp);
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
	} while (ret != 1 && cnt < 10);
	fclose(fp);
	if (cnt >= 10)
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
	} while (ret != 1 && cnt < 10);
	fclose(fp);
	if (cnt >= 10)
		return -9999.0;
	else
		return mu0;
}

static void MakeTitle(const char *inputfn, double t, char *& title)
{
	const int L = strlen(inputfn);
	const int n = L + 80;
	const char *ifn1;
	for (ifn1=inputfn+L-1; ifn1>inputfn; ifn1--) if (*ifn1 == '/') {ifn1++; break;}
	char *tit = new char [n];
	sprintf(tit,"%s, t = %g",ifn1,t);
	title = strdup(tit);
}

inline bool is_noteworthy_cell(Tmetagrid& g, TGridIndex e)
{
	const TCellType ct = g.celltype(e);
	return IncludeGhosts ? (ct == GHOST_CELL || ct == INTERIOR_CELL) : (ct == INTERIOR_CELL);
}

inline void FindCorner(const Tdimvec& Xc, unsigned c, real halfdx, real& x, real& y, real& z)
{
	x = Xc(0) + (GetBit(0,c) ? halfdx : -halfdx);
	y = Xc(1) + (GetBit(1,c) ? halfdx : -halfdx);
	z = Xc(2) + (GetBit(2,c) ? halfdx : -halfdx);
}

static void Phase1(Tmetagrid& g)
// Fill in xyzlist, elemtab, Ne and Np.
{
	TGridIndex i,e,p;
	unsigned c;
	real x,y,z;
	Tdimvec Xc;
	// Compute Ne
	Ne = 0;
	for (i=g.first(); !g.isover(i); i=g.next(i)) {
		if (!g.isleaf(i) || !is_noteworthy_cell(g,i)) continue;
		Ne++;
	}
	if (verbose) clog << "Phase1: Ne=" << Ne << "\n";
	// Allocate elemtab
	elemtab = new TGridIndex [c8*Ne];
	// Fill in xyzlist and elemtab, computing Np at the same time
	xyzlist.init();
	real xmin[3],xmax[3];
	g.getbox(xmin,xmax);
	Txyzhash hashtable(usehash ? SelectHashTableLength(Ne/2) : 10,g.BasegridCellSize(),xmin[0],xmin[1],xmin[2]);
	for (e=0,i=g.first(); !g.isover(i); i=g.next(i)) {
		if (!g.isleaf(i) || !is_noteworthy_cell(g,i)) continue;
		const real dx = g.cellsize(i);
		const real halfdx = 0.5*dx;
		g.centroid(i,Xc);
		for (c=0; c<c8; c++) {
			FindCorner(Xc,c,halfdx,x,y,z);
			if (usehash) {
				p = hashtable.search(x,y,z);
				if (p < 0) {
					p = xyzlist.add(x,y,z);
					hashtable.store(x,y,z,p);
				}
			} else
				p = xyzlist.add(x,y,z);
			elemtabref(e,c) = p;
		}
		e++;
	}
	if (usehash && verbose) hashtable.stats(clog);
	Np = xyzlist.length();
	if (verbose) clog << "Phase1: Np=" << Np << "\n";
}

static void WriteHeader(FILE *fp)
{
	fprintf(fp,"TITLE = \"%s\"\n",title);
	fprintf(fp,"VARIABLES = \"X\",\"Y\",\"Z\"");
	int v;
	for (v=0; v<VAR_LAST; v++) if (enabled[v]) fprintf(fp,",\"%s\"",varnames[v]);
	fprintf(fp,"\n");
	fprintf(fp,"ZONE T=\"%s\", N=%d, E=%d, ET=brick, F=FEPOINT\n\n",zone,Np,Ne);
}

static void WriteValues(FILE *fp, Tmetagrid& g)
{
	const real rho = g.CT(0);
	const real rhoUx = g.CT(1);
	const real rhoUy = g.CT(2);
	const real rhoUz = g.CT(3);
	const real U = g.CT(4);
	real Bx1,By1,Bz1, Bx0,By0,Bz0;
	if (ncd >= 8) {
		Bx1 = g.CT(5);
		By1 = g.CT(6);
		Bz1 = g.CT(7);
	} else {
		Bx1 = By1 = Bz1 = 0;
	}
	if (ncd >= 11) {
		Bx0 = g.CT(8);
		By0 = g.CT(9);
		Bz0 = g.CT(10);
	} else {
		Bx0 = By0 = Bz0 = 0;
	}
	const real Bx = Bx0 + Bx1, By = By0 + By1, Bz = Bz0 + Bz1;

	real vars[VAR_LAST];
	vars[VAR_RHO] = rho;
	vars[VAR_RHOUX] = vars[VAR_RHOVX] = rhoUx;
	vars[VAR_RHOUY] = vars[VAR_RHOVY] = rhoUy;
	vars[VAR_RHOUZ] = vars[VAR_RHOVZ] = rhoUz;
	const real invrho = (rho != 0) ? 1.0/rho : 1.0;
	const real vx=rhoUx*invrho, vy=rhoUy*invrho, vz=rhoUz*invrho;
	vars[VAR_UX] = vars[VAR_VX] = vx;
	vars[VAR_UY] = vars[VAR_VY] = vy;
	vars[VAR_UZ] = vars[VAR_VZ] = vz;
	vars[VAR_U] = U;
	vars[VAR_V2] = sqr(vx) + sqr(vy) + sqr(vz);
	vars[VAR_V] = sqrt(vars[VAR_V2]);
	// 16
	vars[VAR_BX] = Bx;
	vars[VAR_BY] = By;
	vars[VAR_BZ] = Bz;
	vars[VAR_BX1] = Bx1;
	vars[VAR_BY1] = By1;
	vars[VAR_BZ1] = Bz1;
	vars[VAR_BX0] = Bx0;
	vars[VAR_BY0] = By0;
	vars[VAR_BZ0] = Bz0;
	vars[VAR_B2] = sqr(Bx) + sqr(By) + sqr(Bz);
	vars[VAR_B] = sqrt(vars[VAR_B2]);
	vars[VAR_B02] = sqr(Bx0) + sqr(By0) + sqr(Bz0);
	vars[VAR_B0] = sqrt(vars[VAR_B02]);
	vars[VAR_B12] = sqr(Bx1) + sqr(By1) + sqr(Bz1);
	vars[VAR_B1] = sqrt(vars[VAR_B12]);
	const real invmu0=1.0/mu0;
	vars[VAR_P] = (gam-1)*(U - 0.5*rho*vars[VAR_V2] - 0.5*invmu0*vars[VAR_B12]);
	int a;
	for (a=0; a<VAR_LAST; a++) if (enabled[a]) fprintf(fp,format,double(vars[a]));
	fprintf(fp,"\n");
}

static void Phase2(Tmetagrid& g, FILE *fp)
{
	if (verbose) clog << "Phase2. Writing corner coordinate and variable values ...\n";
	Tdimvec Xcorner;
	WriteHeader(fp);
	xyzlist.reset();
	sprintf(format," %%.%dg",digits);
	xyzlist.reset();
	for (; !xyzlist.isend(); xyzlist.next()) {
		Xcorner[0] = xyzlist.X();
		Xcorner[1] = xyzlist.Y();
		Xcorner[2] = xyzlist.Z();
		g.intpol(Xcorner,1,true);
		fprintf(fp,format,double(Xcorner[0]));
		fprintf(fp,format,double(Xcorner[1]));
		fprintf(fp,format,double(Xcorner[2]));
		WriteValues(fp,g);
	}
}

static void Phase3(FILE *fp)
{
	if (verbose) clog << "Phase3. Writing connectivity information\n";
	TGridIndex e;
	fprintf(fp,"\n");
	for (e=0; e<Ne; e++) {
		fprintf(fp,"%d %d %d %d %d %d %d %d\n",
				elemtabref(e,4)+1,elemtabref(e,5)+1,elemtabref(e,7)+1,elemtabref(e,6)+1,
				elemtabref(e,0)+1,elemtabref(e,1)+1,elemtabref(e,3)+1,elemtabref(e,2)+1);
	}
}

int main(int argc, char *argv[])
{
	double t;
	if (argc <= 1) usage();
	varnames[VAR_RHO] = "rho";
	varnames[VAR_RHOUX] = "rhoUx";
	varnames[VAR_RHOUY] = "rhoUy";
	varnames[VAR_RHOUZ] = "rhoUz";
	varnames[VAR_RHOVX] = "rhovx";
	varnames[VAR_RHOVY] = "rhovy";
	varnames[VAR_RHOVZ] = "rhovz";
	varnames[VAR_UX] = "Ux";
	varnames[VAR_UY] = "Uy";
	varnames[VAR_UZ] = "Uz";
	varnames[VAR_VX] = "vx";
	varnames[VAR_VY] = "vy";
	varnames[VAR_VZ] = "vz";
	varnames[VAR_V2] = "v2";
	varnames[VAR_V] = "v";
	varnames[VAR_U] = "U";
	varnames[VAR_P] = "P";
	varnames[VAR_BX] = "Bx";
	varnames[VAR_BY] = "By";
	varnames[VAR_BZ] = "Bz";
	varnames[VAR_BX0] = "Bx0";
	varnames[VAR_BY0] = "By0";
	varnames[VAR_BZ0] = "Bz0";
	varnames[VAR_BX1] = "Bx1";
	varnames[VAR_BY1] = "By1";
	varnames[VAR_BZ1] = "Bz1";
	varnames[VAR_B] = "B";
	varnames[VAR_B0] = "B0";
	varnames[VAR_B1] = "B1";
	varnames[VAR_B2] = "B2";
	varnames[VAR_B02] = "B02";
	varnames[VAR_B12] = "B12";
	enabled[VAR_RHO] = true;
	enabled[VAR_RHOUX] = true;
	enabled[VAR_RHOUY] = true;
	enabled[VAR_RHOUZ] = true;
	enabled[VAR_U] = true;
	enabled[VAR_BX] = true;
	enabled[VAR_BY] = true;
	enabled[VAR_BZ] = true;
	ParseArgs(argc,argv);
	CheckInputFileAndTakeTime(inputfn,t);
	if (!title) MakeTitle(inputfn,t,title);
	if (!zone) zone = "Whole grid";
	FILE *outfp;
	bool outfp_is_file = true;
	if (outputfn) {
		outfp = fopen(outputfn,"w");
		if (!outfp) {
			cerr << "*** hc2tecplot: Could not open output file \"" << outputfn << "\"\n";
			exit(2);
		}
	} else {
		outfp = stdout;
		outfp_is_file = false;
	}
	Theader h(inputfn);
	const int dim = int(h.getint("dim"));
	const double filegamma = FindGammaFromFile(inputfn);
	if (filegamma == -9999.0) {
		if (verbose)
			cerr << "hc2tecplot: gamma value not found in file, assuming gamma=" << gam << "\n";
	} else {
		gam = filegamma;
		if (verbose) cerr << "hc2tecplot: using gamma=" << gam << " from file\n";
	}
	const double filemu0 = FindMu0FromFile(inputfn);
	if (filemu0 > 0) {
		mu0 = filemu0;
		if (verbose)
			cerr << "hc2tecplot: using mu0=" << mu0 << " from file\n";
	}
	if (verbose) cerr << "            (This affects P computation only.)\n";
	if (dim != 3) {
		cerr << "*** hc2tecplot: file \"" << inputfn << "\" dim=" << dim << ", only 3D grids are supported\n";
		exit(3);
	}
	ncd = int(h.getint("ncd"));
	if (ncd != 5 && ncd != 8 && ncd != 11) {
		cerr << "*** hc2tecplot: file \"" << inputfn << "\" ncd=" << ncd << ", only 5,8,11 are supported\n";
		exit(4);
	}
	Tmetagrid *const gridptr = new Tmetagrid(inputfn);
	Phase1(*gridptr);
	Phase2(*gridptr,outfp);
	Phase3(outfp);
	if (outfp_is_file) fclose(outfp);
	return 0;
}
