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

#include "Config.h"
#include "Token.h"
#include "variables.H"
#include "constants.H"
#include <ctype.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <fstream>

using namespace std;


// Trace variables that are not interpolated:
const char* other_vars[] = {"bunEk","bunVx","bunVy","bunVz","parID",NULL};
const char* formats[] = {"vtk","3D","matlab",NULL};


Config::Config()
{
    version = VERSION;

    //Default values:
    maxsteps = 2000;
    stepsize_p = 0.3;
    tracetype = "particle";
    bunemanversion = "U";
    direction = "forward";
    verbose = false;
    overwrite = true;
    endpts = false;
    intpolorder = 1;    
    xmin = -numeric_limits<double>::max();
    ymin = -numeric_limits<double>::max();
    zmin = -numeric_limits<double>::max();
    xmax =  numeric_limits<double>::max();
    ymax =  numeric_limits<double>::max();
    zmax =  numeric_limits<double>::max();
    planetary_boundary = 0.;
    file_formats.push_back(formats[0]);
    out_dir = "./";
    pressuretermfile = "";
    
    // set up trace variables
    tvi_len = 0;
    tvo_len = 0;
    for(int i=0;i<ORBIT_DATA_MAX;i++)
    {
        tvars_intpol[i] = NULL;
        tvars_other[i]  = NULL;
    }
    tvo2id["bunEk"] = BUN_EK; 
    tvo2id["bunVx"] = BUN_VX; 
    tvo2id["bunVy"] = BUN_VY; 
    tvo2id["bunVz"] = BUN_VZ; 
    tvo2id["parID"] = PARID; 
}


bool Config::format_ok(string value)
{
    for(int k=0;formats[k];k++)
        if(value == formats[k])
            return true;
    return false;
}


bool Config::check_other_var(string value)
{
    for(int k=0;other_vars[k];k++)
        if(value == other_vars[k])
            return true;
    return false;
}


void Config::fill_trace_vars2(string line)
{
    Tvariable dummy;
    vector<string> tokens;    
    tokenize(line,tokens);
    tvo_len = 0; tvi_len = 0;

    for(unsigned int i=0; i<tokens.size();i++)
    {
        if(check_other_var(tokens[i]))
        {
            tvars_other[tvo_len] = new char[tokens[i].size()];
            memcpy(tvars_other[tvo_len], tokens[i].c_str(), tokens[i].size()+1);
            tvo_flags[tvo_len] = tvo2id[tokens[i]];
            tvo_len++;             
        }
        else if( dummy.check(tokens[i].c_str())) // check if variable is a hcintpol variable.
        {
            tvars_intpol[tvi_len] = new char[tokens[i].size()];
            memcpy(tvars_intpol[tvi_len], tokens[i].c_str(), tokens[i].size()+1);
            tvi_len++;
        }
        else
        {
            cerr << "Error: unrecognized trace variable." << endl;
            exit(-1);
        }
    }

    if(tvi_len+tvo_len<1)
        cerr << "Error: you must specify a trace variable" << endl, exit(-1);            

    if(tvi_len+tvo_len >= ORBIT_DATA_MAX)
        cerr << "Error: maximum value of trace variables is " << ORBIT_DATA_MAX << endl, exit(-1);
}


void Config::fill_file_formats(string line)
{
    bool found = false;
    vector<string> tokens;    
    tokenize(line,tokens);

    for(unsigned int i=0; i<tokens.size();i++)
    {
        if(not format_ok(tokens[i]))
        {
            cerr << "Error: invalid file format. Check config file." << endl;
            exit(-1);
        }
        if(not found)
        {
            file_formats.clear(); // remove default
            found = true;
        }
        file_formats.push_back(tokens[i]);
    }
}


void Config::fill_hcfile(string line)
{
    vector<string> tokens;    
    tokenize(line,tokens);
    int l = tokens.size();

    if(l==3)
    {
        if(hcfiles.size() > MAX_HCF-1)
            cerr << "Error: to many hcfiles " << endl, exit(-1);

        hcfiles.push_back(tokens[0]);
        double m;
        getDouble(tokens[1],m);
        hcf_mass.push_back(m);
        double q;
        getDouble(tokens[2],q);
        hcf_charge.push_back(q*cnst::qe);

        if(hcfiles.size()!=hcf_mass.size() or hcfiles.size()!=hcf_charge.size()) // sanity check
            cerr << "Error: mass or charge missing in HCF directive" << endl, exit(-1);
        return;
    }
    cerr << "Error: wrong number of arguments to HCF." << endl;
    exit(-1);
}


void Config::readCfg(const char* fname)
{
    ifstream of(fname);
    if(of.fail())
        cerr << "Error: cannot open configuration file" << fname <<"'." << endl, exit(-1);
    of.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );
    cfg_fname = fname;

    char l[310];
    int i; 
    do{  //loop lines.
        try
        {
            of.getline(l,sizeof(l));
        }
        catch(ifstream::failure e)
        {
            if(of.eof())
                break;
            cerr << "Error: IO error reading configuration file! " << e.what() <<  endl;
            exit(-1);
        }

        string var = "";
        string value = "";
        i=0;

        while(isspace(l[i])) i++;             // eat whitespace
        if(l[i]==0 or l[0] == '#') continue;  // is comment
        while(isalnum(l[i]) or l[i] == '_')   // read first token             
            var.push_back(l[i++]);     
        if(var == "EOC") break;               // at end? 
        while(isspace(l[i])) i++;             // eat whitespace

        if(l[i]==0)
            cerr << "Error: configuration option with no value " << endl << "     line:    " << l << endl, exit(-1);
        value = &l[i];

        /* check options */ 
        if(var == "HCF") 
        {
            fill_hcfile(value);
            if(hcfiles.empty())
                cerr << "Error: no hcfiles specified (HCF option)" << endl,exit(-1);
            continue;
        }
        if(var == "MAXSTEPS") 
        {
            maxsteps = atoi(value.c_str());
            if(maxsteps <= 0)
                cerr << "Error: maxsteps must be positive: maxsteps = " << maxsteps  << endl,exit(-1);
            continue;
        }
        if(var == "VERBOSE") 
        {
            verbose = atoi(value.c_str());
            if(verbose != 0 and verbose != 1)
                cerr << "Error: verbose must be 0 or 1 " << verbose << endl,exit(-1);
            continue;
        }
        if(var == "INTPOLORDER") 
        {
            intpolorder = atoi(value.c_str());
            if(intpolorder != 0 and intpolorder != 1)
                cerr << "Error: intpolorder must be 0 or 1 " << verbose << endl,exit(-1);
            continue;
        }
        if(var == "STEPSIZE") 
        {
            stepsize_p = strtod(value.c_str(),NULL);
            if(stepsize_p <= 0)
                cerr << "Error: stepsize must be positive: stepsize_p = " << stepsize_p  << endl,exit(-1);
            continue;
        }
        if(var == "ENDPOINTSONLY") 
        {
            endpts = atoi(value.c_str());
            if(endpts != 0 and endpts != 1)
                cerr << "Error: endpointsonly must be 0 or 1 " << endpts << endl, exit(-1);
            continue;
        }
        if(var == "PRESSURE_TERM_FILE") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            pressuretermfile = value;
            continue;
        }
        if(var == "TRACEVARS") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            fill_trace_vars2(value);
            continue;
        }
        if(var == "OVERWRITE") 
        {
            overwrite = atoi(value.c_str());
            if(overwrite != 0 and overwrite != 1)
                cerr << "Error: overwrite must be 0 or 1 " << overwrite << endl, exit(-1);
            continue;
        }
        if(var == "PLANETARY_BOUNDARY") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            planetary_boundary = strtod(value.c_str(),NULL);
            if(planetary_boundary < 0)
                cerr << "Error: planetary_boundary must be positive: planetary_boundary = " << planetary_boundary  << endl,exit(-1);

            planetary_boundary *= planetary_boundary; 
            continue;
        }
        if(var == "XMIN") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            xmin = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "YMIN") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            ymin = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "ZMIN") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            zmin = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "XMAX") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            xmax = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "YMAX") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            ymax = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "ZMAX") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            zmax = strtod(value.c_str(),NULL);
            continue;
        }
        if(var == "FORMATS")
        {
            fill_file_formats(value);
            continue;
        }
        if(var == "OUT_DIR")
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);
            out_dir = value; 
            continue;
        }
        if(var == "BUNEMANVERSION") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            if( value != "E" and value != "U" and value != "ExB")
                cerr << "Error: bunemanversion should be either E, U or ExB. bunemanversion = " << bunemanversion << endl,exit(-1);

            bunemanversion = value;
            continue;
        }
        if(var == "DIRECTION") 
        {
            while(isspace(*(value.end()-1)))  // remove trailing whitespace
                value.erase(value.end()-1);

            if( value != "forward" and
                value != "backward"  and
                value != "both")
            {
                cerr << "Error: trace direction  must be one of forward, backward or both: direction = " << direction << endl;
                exit(-1);
            }
            direction = value;
            continue;
        }
        cerr << "Error: unrecognized characters in option file"<<endl;
        cerr << "       option: >" << var << "<" << endl;
        cerr << "       value:  >" << value << "<" << endl;
        cerr << "       line:    " << l << endl;
        exit(-1);
    }while(!of.eof());
}


void Config::writeDefaultCfg()
{
//    ofstream cf("iontracer.cfg");
//    if(!cf.good())
//    {
//        cerr << "Error: Cannot open configuration file 'iontracer.cfg' for output." << endl;
//        exit(-1);
//    }
//    cout << "Writing default configuration file to 'iontracer.cfg'" << endl;

    cout << "#####################################################################################################" << endl
        << "#                                iontracer configuration file                                       #" << endl
        << "#                                 iontracer: version " << version << "                                          #" << endl 
        << "#####################################################################################################" << endl <<endl
        << "# HC file, relative mass, relative charge (give separate line for all particle hc files)" << endl
        << "# E.g. HCF H+_hybstate_01100000.hc 1 1  " << endl
        << "#      HCF O+_hybstate_01100000.hc 16 1 " << endl
        << "HCF " << endl << endl

        << "# Electron pressure term file contains -(nabla p) for E field." << endl
        << "# This should be the DBUG file at the same time step (check hc header of hc file if unsure look at B0)." <<endl
        << "# NOTE: comment out if there is no electron pressure term." << endl
        << "# PRESSURE_TERM_FILE " << endl << endl

        << "# File formats specifies the output file format (values: vtk, matlab, 3D)." << endl
        << "FORMATS " << file_formats[0] << endl << endl

        << "# Out put dir specifies a directory where trace files are written. " << endl
        << "OUT_DIR " << out_dir << endl << endl

        << "# Trace variables specifies the data that is saved along the particle orbit (or field line)." << endl
        << "# This should be a space seperated list of the hcintpol variables OR: '"  << endl
        << "#      bunEk                 (kinetic energy of particle)" << endl
        << "#      bunVx, bunVy, bunVz   (velocity component of particle)" << endl
        << "#      parID                 (unique id for each particle)" << endl
        << "TRACEVARS " << endl << endl

        << "# Buneman version: specifies tracing method." << endl
        << "# The alternatives are the two propagators from hyb code or ExB drift (values: E, U or ExB)." << endl
        << "# NOTE: you have to use E if elctron pressure term is included. ExB stands for tracing using ExB drift." << endl
        << "BUNEMANVERSION " << bunemanversion << endl << endl

        << "# Direction of trace (values: forward, backward, both)" << endl
        << "DIRECTION " << direction << endl << endl

        << "# Maximum steps for ODE solvers (integer)" << endl
        << "MAXSTEPS " << maxsteps << endl << endl

        << "# Step size (dt) when doing a particle trace (in seconds)." << endl
        << "STEPSIZE " << stepsize_p << endl << endl

        << "# Interpolation order (0 or 1. 0 gives zero interpolation and 1 linear interpolation.)" << endl
        << "INTPOLORDER " << intpolorder << endl << endl

        << "# Verbose output messages (0 or 1)" << endl
        << "VERBOSE " << verbose << endl << endl

        << "# Overwrite the trace file (0 or 1)" << endl
        << "OVERWRITE " << overwrite << endl << endl

        << "# Save data at end and starting points only." << endl
        << "# This is to save memory when there is alot of particles (0 or 1)" << endl
        << "ENDPOINTSONLY " << endpts << endl << endl

        << "# Planetary boundary spcifies the radie of a origin centerd sphere which acts as " << endl 
        << "# boundary for tracing (if 0 then there is no planetary boundary). (in meters)" << endl
        << "PLANETARY_BOUNDARY " << planetary_boundary << endl << endl

        << "# xmin,...,zmax specifies walls where tracing stops, other than simulation box walls." << endl
        << "# Values are ignored if they are outside the box. (in meters)" << endl
        << "XMIN " << xmin << endl 
        << "YMIN " << ymin << endl 
        << "ZMIN " << zmin << endl 
        << "XMAX  " << xmax << endl 
        << "YMAX  " << ymax << endl 
        << "ZMAX  " << zmax << endl << endl
        << "EOC" << endl << endl

        << "#####################################################################################################" << endl
        << "#                                    INITIAL POINTS SECTION                                         #" << endl
        << "#####################################################################################################" <<endl <<endl
        << "## You can give initial values for tracing bellow" << endl
        << "## (Another option is using the -p option and a seperate file.)" << endl << endl 

        << "## POINT DATA EXAMPLES:" << endl << endl

        << "## Example 1. EXPLCIT DATA POINT" << endl
        << "## Arguments: x, y, z, vx, vy, vz, relative mass and relative charge."  << endl
        << "# 6e6 0 0 0 0 0 1 1"  << endl << endl

        << "## IMPLCIT POINT DATA DIRECTIVES " << endl 
        << "## The bellow examples use directives that assume that the relative mass, relative charge and velocity"  << endl 
        << "## have been set to some value with the 'set' directive."  << endl 
        << "## Default values are mass=1, charge=1 and velocity=(0,0,0)."  << endl << endl

        << "## Example 2: Setting the mass, charge and velocity."  << endl
        << "## This sets the particle mass, charge and velocity for the line, single point and circle directives."  << endl
        << "# set mass 12"  << endl
        << "# set charge 2"  << endl 
        << "# set velocity -400e3 0 0"  << endl << endl

        << "## Example 3: Setting the length scale."  << endl
        << "## This sets the length scale for the coordinate system. The next command "  << endl
        << "## means that x becomes x * 3400e3. Default value is 1."  << endl
        << "# set length_scale 3400e3"  << endl << endl

        << "## Example 4: Data point." << endl
        << "## Arguments: x, y and z."  << endl
        << "# 6e6 0 0"  << endl << endl

        << "## Example 5: Data point in spherical coordinates." << endl
        << "## Arguments: spherical r theta phi (angles in degrees)"  << endl
        << "# spherical 8e6 45 0"  << endl << endl

        << "## Example 6: Line of points between a point a and b."  << endl
        << "## Arguments: ax, ay, az, bx, by, bz and number of points"  << endl
        << "# line 6e6 0 0 -6e6 0 0 30"  << endl << endl

        << "## Example 7: A  circle of points with midpoint (x,y,z)."  << endl
        << "## Arguments: x, y, z, radius, phi, theta, number of points."  << endl
        << "## phi and theta are the spherical coordinates of the unit normal to the plane of the circle."  << endl
        << "## The angles are given in degrees."  << endl
        << "# circle 0 0 0 6e6 0 45 30"  << endl <<endl

        << "## Example 8: randsphere generates a sphere of uniformaly distributed random points."  << endl
        << "## Arguments: originx, originy, originz, radius, radial velocity and number of points."  << endl
        << "# randsphere 0 0 0 7000e3 100e3 1000"  << endl <<endl;

//    cf.close();
}


void Config::writeCfg(ostream &of)
{
    time_t rawtime;
    struct tm *timeinfo;

    time(&rawtime);
    timeinfo = localtime( &rawtime);
    of   << "#" << endl << "# iontracer version " << version << " configuration " << endl 
         << "#   " << asctime(timeinfo) <<"#"  <<endl <<endl;

    for(unsigned int i=0; i<hcfiles.size();i++)
        of <<  "HCF  " << hcfiles[i] <<" "<< hcf_mass[i] <<" "<< (hcf_charge[i]/cnst::qe) <<endl;

    of << "TRACEVARS           "; 
    for(int n=0;n<tvi_len;n++)
        of << tvars_intpol[n] << " ";
    for(int n=0;n<tvo_len;n++)
        of << tvars_other[n] << " ";
    of << endl;

    of << "FORMATS              "; 
    for(unsigned int n=0;n<file_formats.size();n++)
        of << file_formats[n] << " ";
    of << endl;

    of   << "OUT_DIR              " << out_dir << endl
         << "MAXSTEPS             " << maxsteps << endl
         << "STEPSIZE             " << stepsize_p << endl
         << "BUNEMANVERSION       " << bunemanversion << endl
         << "DIRECTION            " << direction << endl
         << "VERBOSE              " << verbose << endl
         << "OVERWRITE            " << overwrite << endl
         << "ENDPOINTSONLY        " << endpts << endl
         << "INTPOLORDER          " << intpolorder << endl
         << "XMIN                 " << xmin << endl
         << "XMAX                 " << xmax << endl
         << "YMIN                 " << ymin << endl
         << "YMAX                 " << ymax << endl
         << "ZMIN                 " << zmin << endl
         << "ZMAX                 " << zmax << endl
         << "PLANETARY_BOUNDARY   " << sqrt(planetary_boundary) << endl; // NOTE! sqrt because
    if(pressuretermfile.size()>0)
        of   << "PRESSURE_TERM_FILE   " << pressuretermfile << endl << endl;
    else
        of   << "#PRESSURE_TERM_FILE"   << endl << endl;
    of   << "EOC" <<endl;
}


Config::~Config()
{
    for(int n=0;n<tvi_len;n++)
        delete tvars_intpol[n];
    for(int n=0;n<tvo_len;n++)
        delete tvars_other[n];
}



