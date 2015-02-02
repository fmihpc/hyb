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
#include <sstream>
#include <stdexcept>
#include <map>
#include <vector>
#include <string>
#include <string.h>
#include <unistd.h>
#include <variables.H>
#include "gridcache.H"

using namespace std;

double Gamma, Invmu0, Mass; bool Pseudobackground;
char *hcfile = NULL;
char *varstr = NULL;
vector<string> vars;
const char *default_vars[] = {"rho","rhovx","rhovy","rhovz","U1","Bx1","By1","Bz1","Bx0","By0","Bz0",NULL};

void usage()
{
    cerr << "usage: hc2gridxyz -f hcfile > outputfile" << endl;
    cerr << "       -f hcfile  input file" << endl;
    cerr << "       -h         help" << endl;
    cerr << "       -v         comma separated list of variables (see: hcintpol -fullhelp)" << endl;
    cerr 
     << "hc2gridxyz can be used to create an ascii grid file from a hc" << endl
     << "file. The grid file has the columns:" << endl << endl
     << "xc yc zc reflevel var1 var2 var3..." << endl << endl
     << "where xc, yc and zc are centroind coodinates of the cell," << endl
     << "reflevel the refinement level of the cell (0=base cell," << endl
     << "1=first refined level etc), varX are the variables used asked" << endl
     << "to get from the hc file." << endl << endl
     << "Only leaf cells (i.e. cells actually having data in them)" << endl
     << "are listed, no parent cells for refined children cells are" << endl
     << "stored." << endl;
     
    exit(0);
}

// parse variables given in cmd line
void parseVars(char* varlist)
{       
    char *var = strtok(varlist,",");
    while(var)
    {
        vars.push_back(var);
        var = strtok(NULL,",");
    }
}

// parse cmd line options
void options(int argc, char**argv)
{
    char c;
    opterr = 0;
    while((c=getopt(argc,argv,"hf:v:m:")) != -1)
    {
        switch (c)
        {
            case 'f':
                hcfile = optarg;
                break;
            case 'h':
                usage();
	        break;
            case 'v':
                parseVars(optarg);
                break;
            default:
	        cerr << "Error: unknown option" << endl;
	        exit(-1);
        }
    }
    if(!hcfile)
        usage();
    if(vars.size()==0)
        for(int i=0;default_vars[i];i++)
            vars.push_back(default_vars[i]);
}

// check is a cell is "extra", that is, a ghost, removed or dead cell
inline bool isExtra(Tmetagrid &g, TGridIndex &c)
{
    const TCellType ct = g.celltype(c);
    if(ct == DEAD_CELL or ct == GHOST_CELL or ct == REMOVED_CELL) 
        return true;
    return false;
}

int main(int argc, char **argv)
{
    // parse cmd line options and open a hc file
    TGridCache gridcache;
    options(argc,argv);
    Tmetagrid* gp = gridcache.open(hcfile,Gamma,Invmu0,Mass,Pseudobackground);
    if(!gp) {cerr << "*** hcintpol: cannot open HC file \"" << hcfile << "\"\n"; exit(4);}
    Tmetagrid& g = *gp;

    cout << scientific;
    cout.precision(5);

    // write coordinates of cells, refinement level and asked variables in the file
    TGridIndex c;
    Tdimvec Xc;
    Tvariable var;
    int totcells=0;
    int reflevel=0;
    for (c=g.first(); !g.isover(c); c=g.next(c)) 
    {
       //if(isExtra(g,c)) { continue; }       
       if(!g.isleaf(c)) { reflevel++; continue; } // check if the cell has children
       g.centroid(c,Xc); // get centroid coordinates of the cell
       g.intpol(Xc,0,true); // get variables in the cell
       cout << Xc[0] << " " << Xc[1] << " " << Xc[2] << " " << reflevel << " ";
       for(unsigned int i=0; i<vars.size();i++)
	 {
	    var.select(vars[i].c_str(),Gamma,Invmu0,Mass);
	    cout << var.get(g,Xc) << " ";
	 }
       cout << endl;
       totcells++;
       reflevel=0;
    }
   return 1;
}

