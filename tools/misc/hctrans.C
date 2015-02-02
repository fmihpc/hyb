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

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "metagrid.H"
using namespace std;

static void usage()
{
	clog <<
 "Usage: hctrans [-ascii | -float | -double ] [-remove_gaps] [-chop] [-precision prec]\n"
 "               input.hc output.hc\n"
 "reads a HC file and writes another HC file in the given format.\n"
 "The default output format is ASCII. Default ASCII precision is 10 digits.\n";
	exit(0);
}

static int outputformat = 0;		// 0:ascii, 1:float, 2: double
static bool remove_gaps = false;
static bool chopping = false;
static int precision = 10;

#if defined(_UNICOS)
#  define finite(x) (-DBL_MAX < (x) && (x) < DBL_MAX)
#endif

int main(int argc, char *argv[])
{
	if (argc < 3) usage();
	int a;
	
	for (a=1; a<argc && argv[a][0] == '-'; a++) {
		if (!strcmp(argv[a],"-ascii")) {
			outputformat = 0;
		} else if (!strcmp(argv[a],"-float")) {
			outputformat = 1;
		} else if (!strcmp(argv[a],"-double")) {
			outputformat = 2;
		} else if (!strcmp(argv[a],"-remove_gaps")) {
			remove_gaps = true;
		} else if (!strcmp(argv[a],"-chop")) {
			chopping = true;
		} else if (!strcmp(argv[a],"-precision") || !strcmp(argv[a],"-prec")) {
			precision = atoi(argv[++a]);
		} else if (!strcmp(argv[a],"-help") || !strcmp(argv[a],"--help")) {
			usage();
		}
	}
	if (argc != a+2) usage();
	char *input = argv[a];
	char *output = argv[a+1];
	Tmetagrid g(input);
	TRealFormat rf=REALFORMAT_ASCII;
	switch (outputformat) {
	case 0: rf = REALFORMAT_ASCII; break;
	case 1: rf = REALFORMAT_FLOAT; break;
	case 2: rf = REALFORMAT_DOUBLE; break;
	}
	g.realformat(rf);
	g.set_remove_gaps(remove_gaps);
	if (chopping) {
		const real chop_epsilon = pow(10.0,real(-precision));
		TGridIndex i;
		const smallnat ncd = g.Ncelldata();
		smallnat a;
		for (i=g.first(); !g.isover(i); i=g.next(i)) {
			g > i;
			for (a=0; a<ncd; a++)
				if (!finite(g.CT(a)) || fabs(g.CT(a)) < chop_epsilon) g.CT(a) = 0;
			g << i;
		}
	}
	// Copy comment lines from the beginning, Need hchead on search path
	char *s = new char [80 + strlen(input) + strlen(output)];
	sprintf(s,"hchead %s | fgrep '#' >%s",input,output);
	system(s);
	g.save(output,precision,ios::app);
	return 0;
}
