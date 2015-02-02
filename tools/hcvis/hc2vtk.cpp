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
map<string,int> table;
vector<string> points;
vector<string> vars;
const char *default_vars[] = {"rho","rhovx","rhovy","rhovz","U1","Bx1","By1","Bz1","Bx0","By0","Bz0",NULL};

void usage()
{
    cerr << "usage: hc2vtk -f hcfile > outputfile" << endl;
    cerr << "       -f hcfile  input file" << endl;
    cerr << "       -h         help" << endl;
    cerr << "       -v         comma separated list of variables (see: hcintpol -fullhelp)" << endl;
    exit(0);
}

void terminate(const char *err)
{
    cerr << "Error: " << err << endl;
    exit(-1);
}

void parseVars(char* varlist)
{       
    char *var = strtok(varlist,",");
    while(var)
    {
        vars.push_back(var);
        var = strtok(NULL,",");
    }
}

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
                terminate("unknown option");
        }
    }
    if(!hcfile)
        usage();
    if(vars.size()==0)
        for(int i=0;default_vars[i];i++)
            vars.push_back(default_vars[i]);
}

inline string coord2str(const double &X, const double &Y, const double &Z, int mod = 10)
{
    std::ostringstream o;
    o.precision(5);
//    if (!(o << x-(x%mod) <<' '<< y-(y%mod) <<' '<< z-(z%mod)))
    if (!(o << X <<' '<< Y <<' '<< Z))
        throw runtime_error("cord2str: could not convert coordinate to string");
    return o.str();
}

inline void cellIndex(map<string,int> &table, string &key)
{
    map<string,int>::iterator it = table.find(key);
    if(it == table.end())
        terminate("coordinate not in table");
    cout << it->second << ' ';
}

inline void putToTable(string &key)
{
    static int n = 0;
    if(table.find(key) == table.end())
    {
        table[key] = n++; 
        points.push_back(key);
    }
}

inline bool isExtra(Tmetagrid &g, TGridIndex &c)
{
    const TCellType ct = g.celltype(c);
    if(ct == DEAD_CELL or ct == GHOST_CELL or ct == REMOVED_CELL) 
        return true;
    return false;
}

void writeScalarData(Tmetagrid &g, int totcells)
{
    Tvariable var;
    Tdimvec Xc;
    TGridIndex c;
    cout << "CELL_DATA " << totcells << endl;
    for(unsigned int i=0; i<vars.size();i++)
    {
        cout << "SCALARS " << vars[i]<< " float 1" << endl;
        cout << "LOOKUP_TABLE default" << endl;
        for (c=g.first(); !g.isover(c); c=g.next(c)) 
        {
            if(!g.isleaf(c) or isExtra(g,c)) 
                continue;
            g.centroid(c,Xc);
            g.intpol(Xc,0,true);
            var.select(vars[i].c_str(),Gamma,Invmu0,Mass);
            cout << var.get(g,Xc) << endl;
        }
    }
}

int main(int argc, char **argv)
{
    string key;
    TGridCache gridcache;
    cout.precision(5);

    options(argc,argv);
    Tmetagrid* gp = gridcache.open(hcfile,Gamma,Invmu0,Mass,Pseudobackground);
    if(!gp) {cerr << "*** hcintpol: cannot open HC file \"" << hcfile << "\"\n"; exit(4);}
    Tmetagrid& g = *gp;

    TGridIndex c;
    Tdimvec Xc;
    int totcells=0;
    for (c=g.first(); !g.isover(c); c=g.next(c)) 
    {
        if(!g.isleaf(c) or isExtra(g,c)) 
            continue;

        totcells++;
	g.centroid(c,Xc);
        double halfdx = 0.5*g.cellsize(c);
        key = coord2str( Xc[0]-halfdx,Xc[1]-halfdx,Xc[2]-halfdx); putToTable(key);
        key = coord2str( Xc[0]+halfdx,Xc[1]-halfdx,Xc[2]-halfdx); putToTable(key);
        key = coord2str( Xc[0]-halfdx,Xc[1]+halfdx,Xc[2]-halfdx); putToTable(key);
        key = coord2str( Xc[0]+halfdx,Xc[1]+halfdx,Xc[2]-halfdx); putToTable(key);
        key = coord2str( Xc[0]-halfdx,Xc[1]-halfdx,Xc[2]+halfdx); putToTable(key);
        key = coord2str( Xc[0]+halfdx,Xc[1]-halfdx,Xc[2]+halfdx); putToTable(key);
        key = coord2str( Xc[0]-halfdx,Xc[1]+halfdx,Xc[2]+halfdx); putToTable(key);
        key = coord2str( Xc[0]+halfdx,Xc[1]+halfdx,Xc[2]+halfdx); putToTable(key);
    }

    cout << "# vtk DataFile Version 2.0" << endl;
    cout << "HC file grid" << endl;
    cout << "ASCII" << endl;
    cout << "DATASET UNSTRUCTURED_GRID" << endl;

    // POINT_DATA 
    int totpoints = table.size();
    cout << "POINTS " << totpoints << " float " << endl;
    for(unsigned int i=0; i<points.size();i++)
        cout << points[i] <<endl;

    // CELL_DATA
    cout << "CELLS " << totcells << " " << 9*totcells << endl;
    map<string,int>::iterator it;
    for (c=g.first(); !g.isover(c); c=g.next(c)) 
    {
        if(!g.isleaf(c) or isExtra(g,c)) 
            continue;

        cout << "8 ";
	g.centroid(c,Xc);
        double halfdx = 0.5*g.cellsize(c);
        key = coord2str( Xc[0]-halfdx,Xc[1]-halfdx,Xc[2]-halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]+halfdx,Xc[1]-halfdx,Xc[2]-halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]-halfdx,Xc[1]+halfdx,Xc[2]-halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]+halfdx,Xc[1]+halfdx,Xc[2]-halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]-halfdx,Xc[1]-halfdx,Xc[2]+halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]+halfdx,Xc[1]-halfdx,Xc[2]+halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]-halfdx,Xc[1]+halfdx,Xc[2]+halfdx); cellIndex(table,key);
        key = coord2str( Xc[0]+halfdx,Xc[1]+halfdx,Xc[2]+halfdx); cellIndex(table,key);
        cout << endl;
    }
    cout << "CELL_TYPES " << totcells << endl;
    for(int i=0; i<totcells;i++)
        cout << "11" << endl;
    
    writeScalarData(g,totcells);
    return 1;
}

