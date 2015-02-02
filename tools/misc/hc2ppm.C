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

#include <cstdio>
#include <cstring>
#include <ctime>
#include <cstdlib>
#include <cmath>
#include "metagrid.H"
using namespace std;

int main(int argc, char *argv[])
{
	int verbose=1;			// Verbose flag (0:silent, 1:normal, 2:verbose)
	int arg,start=1;
	int zflag = 0;			// Gzip flag
	int gflag = 0;			// GIF flag
	int linflag = 0;		// Linear interpolation flag (default is zeroth order "interpolation")
	int ignore_boundary = 1;
//	int zero_dead = 1;		// Put dead cell states to all zero
	int multiplex = 1;		// Pixel replication factor, 1 gives original size
	int a = 0;
	int slicing_dim = -1;		// slicing_dim=0,1,2 for x,y,z; -1 represents for default (see below)
	double slicing_xval = 0;
//	Tcell *c;
//	FILE *fp;
//	Tmodel::SetFatalAbortFlag(0);
//	Tmodel::SetErrorReportFlag(0);
	if (argc <= 1 || argc >= 2 && !strcmp(argv[1],"-help")) {
		cerr << "usage: hc2ppm [-zg01iV] [-v n] [-s scale] [-X xval] [-Y yval] [-Z zval]\n";
		cerr << "              file1.hc ...\n";
		cerr << "  creates file1.ppm, file1.ppm ... from HC grid file(s) file1.hc ...\n";
		cerr << "  The resulting files can be displayed e.g. by xv.\n";
		cerr << "  -z compresses the files with gzip. They are named *.ppm.gz in that case.\n";
		cerr << "  -g write GIF files instead of PPMs. They are named *.gif. ppmtogif is needed.\n";
		cerr << "  -0 uses zeroth order interpolation (the default)\n";
		cerr << "  -1 uses linear interpolation\n";
		cerr << "  -i include boundary and dead cells also\n";
		cerr << "  -V Verbose mode\n";
		cerr << "  -v n     select nth variable (MHD: n=1..8)\n";
		cerr << "  -s n     use n (n=1,2,3...) pixels per smallest cell (default n=1)\n";
		cerr << "  -X xval  use x=xval slice in case of 3D HC files\n";
		cerr << "  -Y yval  use y=yval slice\n";
		cerr << "  -Z zval  use z=zval slice\n";
		cerr << "The default is to slice in Z=0.5*(zmin+zmax) plane.\n";
		return 0;
	}
	while (argv[start][0] == '-') {
		if (strchr("zg01iV",argv[start][1])) {
			int i;
			for (i=1; argv[start][i]; i++) {
				if (argv[start][i] == 'z') {
					zflag = 1;
					gflag = 0;
				} else if (argv[start][i] == 'g') {
					gflag = 1;
					zflag = 0;
				} else if (argv[start][i] == '0') {
					linflag = 0;
				} else if (argv[start][i] == '1') {
					linflag = 1;
				} else if (argv[start][i] == 'i') {
					ignore_boundary = 0;
				} else if (argv[start][i] == 'V') {
					verbose = 2;
				} else {
					cerr << "### h hc2ppm: Unrecognized option flag -" << argv[start][i] << " ignored\n";
				}
			}
		} else if (argv[start][1] == 'v') {
			a = atoi(argv[start+1]) - 1;
			start++;
		} else if (argv[start][1] == 's') {
			multiplex = atoi(argv[start+1]);
			if (multiplex <= 0 || multiplex > 50) {
				cerr << "### Bad -s argument " << multiplex << " - multiplexing reset to unity\n";
				multiplex = 1;
			}
			start++;
		} else if (argv[start][1] == 'X') {
			slicing_dim = 0;
			slicing_xval = atof(argv[++start]);
		} else if (argv[start][1] == 'Y') {
			slicing_dim = 1;
			slicing_xval = atof(argv[++start]);
		} else if (argv[start][1] == 'Z') {
			slicing_dim = 2;
			slicing_xval = atof(argv[++start]);
		} else {
			cerr << "### hc2ppm: Unrecognized option \"" << argv[start] << "\" ignored\n";
		}
		start++;
	}
	for (arg=start; arg<argc; arg++) {
		const int L = strlen(argv[arg]);
		char *ppmfile = new char [L+5];
		if (argv[arg][L-1] == 'c' && argv[arg][L-2] == 'h' && argv[arg][L-3] == '.') {
			strcpy(ppmfile,argv[arg]);
			ppmfile[L-2] = 'p';
			ppmfile[L-1] = 'p';
			ppmfile[L] = 'm';
			if (zflag) {
				ppmfile[L+1] = '.';
				ppmfile[L+2] = 'g';
				ppmfile[L+3] = 'z';
				ppmfile[L+4] = '\0';
			} else if (gflag) {
				ppmfile[L-2] = 'g';
				ppmfile[L-1] = 'i';
				ppmfile[L] = 'f';
				ppmfile[L+1] = '\0';
			} else
				ppmfile[L+1] = '\0';
		} else {
			cerr << "*** hc2ppm: file name \"" << argv[arg] << "\" does not have suffix .hc - ignored\n";
			delete [] ppmfile;
			continue;
		}
		if (verbose >= 2) cout << "Loading \"" << argv[arg] << "\" ..." << flush;
//		if (verbose != 1) Tgrid_HC::silent = 1;
//		Tgrid_HC::fastload = 1;
		clock_t t1,t2;
		double cputime;
		t1 = clock();
		Tmetagrid g(argv[arg]);
		t2 = clock();
		cputime = (t2-t1)/double(CLOCKS_PER_SEC);
		if (verbose >= 2) cout << " done, " << cputime << " seconds\n" << flush;
		if (!g.good()) {
			cerr << "*** hc2ppm: File \"" << argv[arg] << "\" ignored\n";
			continue;
		}
		if ((g.dimension() != 2 && g.dimension() != 3) || g.dimension() > MAXDIM ) {
			if (g.dimension() > MAXDIM)
				cerr << "*** hc2ppm: MAXDIM=" << MAXDIM << " is too small to read this HC file - recompile!\n";
			else
				cerr << "*** hc2ppm: file \"" << argv[arg] << "\" is not 2- or 3-dimensional - ignored\n";
			continue;
		}
#if 0
		if (zero_dead) {
			Tcell *c;
			for (c=g.FirstCell(); c; c=c->Next())
				if (c->CellType() == DEAD_CELL) c->u = 0;
		}
#endif
		real xmin,xmax,ymin,ymax,zmin=0,zmax=0;
//		if (ignore_boundary) {
			real Xmin[3],Xmax[3];
			g.getbox(Xmin,Xmax);
			xmin = Xmin[0]; xmax = Xmax[0];
			ymin = Xmin[1]; ymax = Xmax[1];
			zmin = Xmin[2]; zmax = Xmax[2];
//			g.get_interior_domain_minmax(0,xmin,xmax);
//			g.get_interior_domain_minmax(1,ymin,ymax);
//			if (g.dimension()==3) g.get_interior_domain_minmax(2,zmin,zmax);
//		} else {
//			g.get_domain_minmax(0,xmin,xmax);
//			g.get_domain_minmax(1,ymin,ymax);
//			if (g.dime*nsion()==3) g.get_domain_minmax(2,zmin,zmax);
//		}
		const real mindx = g.MinimumGridSpacing()/multiplex;
		const int nx=int((xmax-xmin)/mindx+0.5);
		const int ny=int((ymax-ymin)/mindx+0.5);
//		const int nz=int((zmax-zmin)/mindx+0.5);
		real *Mr = new real [nx*ny];
		int i,j;
		const int NV = g.Ncelldata();
//		TNVvec u;
		real *u = new real [NV];
		Tdimvec X;
		real maxMr=0,minMr=0;
		if (a < 0) a = 0;
		if (a >= NV) a = NV-1;
		if (g.dimension()==3) {
			if (slicing_dim == -1) {
				slicing_dim = 2;			// by default, slice in Z=0.5*(zmin+zmax) plane
				slicing_xval = 0.5*(zmin + zmax);
			}
			if (verbose >= 1) {
				cout << "Slicing 3D dataset in ";
				switch (slicing_dim) {
				case 0: cout << "X";break;
				case 1: cout << "Y";break;
				case 2: cout << "Z";break;
				}
				cout << "=" << slicing_xval << " plane\n";
			}
		}
		if (verbose >= 2) cout << "Interpolating ..." << flush;
		t1 = clock();
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
			const real x = xmin + (i+0.5)*mindx;
			const real y = ymin + (j+0.5)*mindx;
			if (g.dimension()==3) {
				switch (slicing_dim) {
				case 0:
					X[0] = slicing_xval;
					X[1] = x;
					X[2] = y;
					break;
				case 1:
					X[0] = x;
					X[1] = slicing_xval;
					X[2] = y;
					break;
				case 2:
					X[0] = x;
					X[1] = y;
					X[2] = slicing_xval;
					break;
				}
			} else {
				X[0] = x;
				X[1] = y;
			}
			if (linflag)
				g.intpol(X,1);
			else
				g.intpol(X,0);
			smallnat a1;
			for (a1=0; a1<NV; a1++) u[a1] = g.CT(a1);
			Mr[i*ny+j] = u[a];
			if (i==0 && j==0)
				maxMr = minMr = u[a];
			else {
				maxMr = max(maxMr,u[a]);
				minMr = min(minMr,u[a]);
			}
		}
		cout << flush;
		t2 = clock();
		cputime = (t2-t1)/double(CLOCKS_PER_SEC);
		if (verbose >= 2) cout << " done, " << cputime << " seconds\n" << flush;
		if (maxMr <= minMr) maxMr = minMr + 1;
		unsigned char *M = new unsigned char [nx*ny];
		const real coeff = 256.0/(maxMr - minMr);
		for (i=0; i<nx; i++) for (j=0; j<ny; j++) {
			real Mscaled = (Mr[i*ny+j] - minMr)*coeff;
			if (Mscaled > 255.5) Mscaled = 255.5;
			if (Mscaled < 0.5) Mscaled = 0.5;
			M[i*ny+j] = int(floor(Mscaled));
		}
		FILE *fp;
		if (zflag) {
			char pipename[1024];
			sprintf(pipename,"gzip >%s",ppmfile);
			fp = popen(pipename,"w");
		} else if (gflag) {
			char pipename[1024];
			sprintf(pipename,"ppmtogif >%s",ppmfile);
			fp = popen(pipename,"w");
		} else
			fp = fopen(ppmfile,"w");
		if (!fp) {
			cerr << "*** Could not open output\n";
			exit(1);
		}
		if (verbose >= 2) cout << "Writing output \"" << ppmfile << "\" ..." << flush;
		fprintf(fp,"P5\n");
		fprintf(fp,"# 8-bit PGM file created from \"%s\". Created by hc2ppm.\n",argv[arg]);
		fprintf(fp,"%d %d\n255\n",nx,ny);
		for (j=0; j<ny; j++) for (i=0; i<nx; i++)
			putc(M[i*ny+j],fp);
		if (zflag || gflag)
			pclose(fp);
		else
			fclose(fp);
		if (verbose >= 2)
			cout << " done.\n" << flush;
		else if (verbose == 1)
			cout << "hc2ppm: Wrote \"" << ppmfile << "\"\n";
		delete [] M;
		delete [] Mr;
		delete [] ppmfile;
		delete [] u;
	}
	return 0;
}
