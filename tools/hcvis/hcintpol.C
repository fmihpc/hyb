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
#include <iomanip>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include "variables.H"
#include "gridcache.H"

FILE *input = 0, *output = 0;
static char *hcfile = 0;

static void usage()
{
	cerr << "usage: hcintpol [-n] [-v varlist] hcfile [<]pointfile >output\n";
	cerr << "  interpolates at the given points from the given HC file.\n";
	cerr << "The pointfile must have 3 coordinates (x y z) per line.\n";
	cerr << "The coordinates must be given in meters (unless the simulation\n";
	cerr << "  follows untypical conventions).\n";
	cerr << "If y and/or z are not given, they are assumed zero.\n";
	cerr << "If the grid is 1D or 2D only the 1 or 2 first coordinates are used.\n";
	cerr << "The ASCII output by default consists of all the variables.\n";
	cerr << "A subset of variables can be selected with -v var1,var2,...,varN\n";
	cerr << "By default the first line is a '#' comment line listing";
	cerr << " the variable names.\n";
	cerr << "Options: -n            Supress the comment line\n";
	cerr << "         -v varlist    Define variables to extract\n";
	cerr << "         -i            Ignore ghost cells (zero data for them)\n";
	cerr << "         -z            use Zeroth order interpolation (default first order)\n";
	cerr << "NOTE: The variables are not necessarily in the asked order!\n";
	cerr << "      The '#' comment line tells the ordering.\n";
	cerr << "Use hcintpol -fullhelp to get list of variables and their descriptions.\n";
}

int main(int argc, char *argv[])
{
	int a = 1;
	int intpol_order = 1;		// 1=linear, 0=zeroth order
	bool SupressComment = false;
	bool IgnoreGhost = false;
	char *VariableList = 0;
	if (argc < 2) {usage(); exit(1);}
	while (!strcmp(argv[a],"-n") || !strcmp(argv[a],"-v") || !strcmp(argv[a],"-i") || !strcmp(argv[a],"-z") || !strcmp(argv[a],"-fullhelp")) {
		if (!strcmp(argv[a],"-n")) {
			SupressComment = true;
			a++;
		} else if (!strcmp(argv[a],"-i")) {
			IgnoreGhost = true;
			a++;
		} else if (!strcmp(argv[a],"-v")) {
			VariableList = argv[a+1];
			a+= 2;
		} else if (!strcmp(argv[a],"-z")) {
			intpol_order = 0;
			a++;
		} else if (!strcmp(argv[a],"-fullhelp")) {
			usage();
			Tvariable var;
			const int n = var.Nvars();
			int maxwidth = strlen("Variable")-1;
			for (int i=0; i<n; i++) {
				var.select(i);
				const int L = strlen(var.selected());
				if (L > maxwidth) maxwidth = L;
			}
			cerr << '\n';
			cerr << setw(maxwidth+2) << left << "Variable" << "Description\n";
			cerr << setw(maxwidth+2) << left << "--------" << "-----------\n";
			for (int i=0; i<n; i++) {
				var.select(i);
				cerr << setw(maxwidth+2) << left << var.selected();
				char *s = new char [strlen(var.description()) + 1];
				strcpy(s,var.description());
				for (int j=0; s[j]; j++) if (s[j] == '\n') s[j] = ' ';
				cerr << s << '\n';
				delete [] s;
			}
			return 0;
		}
	}
	hcfile = argv[a];
	input = stdin;
	if (argc-a == 2) {
		input = fopen(argv[a+1],"r");
		if (!input) {cerr << "*** hcintpol: cannot open input file \"" << argv[a+1] << "\"\n"; exit(2);}
	} else if (argc-a != 1) {usage(); exit(3);}
	output = stdout;
	char s[1026];
	TGridCache gridcache;
	double Gamma, Invmu0, Mass; bool Pseudobackground;
	Tmetagrid* gp = gridcache.open(hcfile,Gamma,Invmu0,Mass,Pseudobackground);
	if (!gp) {cerr << "*** hcintpol: cannot open HC file \"" << hcfile << "\"\n"; exit(4);}
        bool isSpectraFile = gridcache.isSpectraFile();
        if(isSpectraFile == true) { intpol_order = 0; }
	Tmetagrid& g = *gp;
	int warn = 0, ok = 0;
	Tvariable var;
        if(var.Nvars() > MAX_VARS) { cerr << "*** hcintpol: too many variables in variables.C (\"" << hcfile << "\")\n"; exit(5); }
        if(g.Ncelldata() > MAX_VARS) {
	   if(isSpectraFile == false) {
	      cerr << "*** hcintpol: too many data columns in \"" << hcfile << "\"\n";
	   } else {
	      cerr << "*** hcintpol: too many energy bins in \"" << hcfile << "\"\n";
	   }
	   exit(5);
	}
	// Find out vars
	int i;
	char *varnames[MAX_VARS];
	double varvalues[MAX_VARS];
	bool varflags[MAX_VARS];
	for (i=0; i<var.Nvars(); i++) {var.select(i); varnames[i] = strdup(var.selected()); varflags[i] = true;}
	if (VariableList && isSpectraFile == false) {
		for (i=0; i<var.Nvars(); i++) varflags[i] = false;
		int ntokens = 1;
		for (i=0; VariableList[i]; i++) if (VariableList[i] == ',') ntokens++;
		if (ntokens >= MAX_VARS) {cerr << "*** Too many variables selected\n"; exit(8);}
		int commapos[MAX_VARS];
		int j = 0;
		for (i=0; VariableList[i]; i++) if (VariableList[i] == ',') commapos[j++] = i;
		int tok_start[MAX_VARS],tok_end[MAX_VARS];
		tok_start[0] = 0;
		for (j=0; j<ntokens-1; j++) {
			tok_start[j+1] = commapos[j]+1;
			tok_end[j] = commapos[j]-1;
		}
		tok_end[ntokens-1] = strlen(VariableList)-1;
		char ss[MAX_VARS];
		bool errors = false;
		for (j=0; j<ntokens; j++) {
			int k = 0;
			for (i=tok_start[j]; i<=tok_end[j]; i++) ss[k++] = VariableList[i];
			ss[k] = '\0';
			bool found = false;
			for (i=0; i<var.Nvars(); i++) if (!strcmp(varnames[i],ss)) {varflags[i] = true; found = true; break;}
			if (!found) {cerr << "*** Unrecognized variable name \"" << ss << "\", ignored\n"; errors=true;}
		}
		if (errors) {
			cerr << "Valid variables are: ";
			for (i=0; i<var.Nvars()-1; i++) cerr << varnames[i] << ',';
			cerr << varnames[var.Nvars()-1] << "\n";
		}
	}
	if (!SupressComment) {
		fprintf(output,"# x y z");
	        if(isSpectraFile == false) {
		   for (i=0; i<var.Nvars(); i++) if (varflags[i]) fprintf(output," %s",varnames[i]);
		} else {
		   for (i=0; i<g.Ncelldata(); i++) fprintf(output," Ebin%d",i);
		}
 	        fprintf(output,"\n");
	}
	while (!feof(input)) {
	        fgets(s,1025,input);
		if (feof(input)) break;
		const int L = strlen(s);
		if (s[L-1] != '\n') {cerr << "*** hcintpol: too long input line (>1024)\n"; exit(5);}
		s[L-1] = '\0';		// removed newline
		double x=0,y=0,z=0;
		const int scanret = sscanf(s,"%lf%lf%lf",&x,&y,&z);
		if (scanret < 1) {cerr << "*** hcintpol: syntax error in input file\n"; exit(6);}
		Tdimvec X(x,y,z);
	        if(isSpectraFile == false) {
			if (!g.intpol(X,intpol_order,true))
		    {
		      warn++;
		      for (i=0; i<var.Nvars(); i++) varvalues[i] = -999;
		    }
		    else
		    {
				if (IgnoreGhost) 
				{
					const TGridIndex ii = g.find(X);
					if (g.celltype(ii) != INTERIOR_CELL)
					{
						for (i=0; i<var.Nvars(); i++) varvalues[i] = -999;
						goto writefile;
					}
				}

		      ok++;
		      for (i=0; i<var.Nvars(); i++)
			  {
					if (!varflags[i]) {varvalues[i] = -999; continue;}
					var.select(varnames[i],Gamma,Invmu0,Mass);
					varvalues[i] = var.get(g,X);
		      }
		   }
writefile:
		   fprintf(output,"%g %g %g",x,y,z);
		   for (i=0; i<var.Nvars(); i++) if (varflags[i]) fprintf(output," %g",varvalues[i]);
		   fprintf(output,"\n");
		}
	        else { // spectra file
		   fprintf(output,"%g %g %g",x,y,z);
		   if (!g.intpol(X,intpol_order,true)) {
		      warn++;
		      for (i=0; i<g.Ncelldata(); i++) { fprintf(output," -999"); }
		      fprintf(output,"\n");
		   }
		   else {
		      vector<double> spectra = var.getSpectra(g,X);
		      for (unsigned int k=0;k<spectra.size();k++) { fprintf(output," %g",spectra[k]); }
		      fprintf(output,"\n");
		   }
		}
	}
	if (warn > 0) cerr << "- hcintpol: " << warn << "/" << (ok+warn) << " points out of domain (set as -999 in output)\n";
	return 0;
}
