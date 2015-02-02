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

#include<iostream>
#include<vector>
#include<string>
#include<unistd.h>
#include<stdlib.h>
#include"Track.h"
#include"Config.h"
#include"Ptracer.h"

using namespace std;

void usage(int r, const char *err)
{
    cout << "usage: iontracer configfile" << endl;
    cout << "       -p file    use pointfile instead of point data in config file" << endl;
    cout << "                  Format: text file with similiar directives as in end of config file" << endl;
    cout << "                          (See examples at end of config file.)" << endl;
    cout << "       -c         write default configuration file to stdout" << endl;
    cout << "       -h         help" << endl;
    cout << "       -v         print version" << endl;
    cout << err << endl;
    exit(r);
}


int main(int argc, char** argv)
{
    char c;
    Config cfg;
    string cfgfile;
    char *pointfile = NULL;

    while((c=getopt(argc, argv, "hvp:c")) != -1)
    {
        switch(c)
        {
            case 'h':
                usage(0,"");
                break;
            case 'p':
                pointfile = optarg;
                break;
            case 'v':
                cout << cfg.version << endl;
                exit(0);
                break;
            case 'c':
                cfg.writeDefaultCfg();
                exit(0);
            default:
                usage(3,"Option error");
        };
    }
    if(!argv[optind]) usage(-1,"  Error: Please specify a configuration file (to create one use the -c option).");
    cfgfile = argv[optind];

    cfg.readCfg(cfgfile.c_str());
    if(cfg.verbose) cfg.writeCfg(cout);

    Ptracer tracer(&cfg);
    if(cfg.pressuretermfile.size() == 0 and cfg.bunemanversion == "E" and cfg.verbose)
        cerr << "\nWARNING!!!: no pressure term file given.\n" <<endl;

    if(pointfile == NULL)        
        tracer.readInitialPoints(cfgfile.c_str(),true);
    else
        tracer.readInitialPoints(pointfile,false);

    tracer.trace();       
    tracer.write();       
    return 0;
}


