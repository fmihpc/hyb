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

#include "simulation.h"
#include "params.h"

using namespace std;

extern Params simuConfig;

//! Show program usage
void showUsage()
{
    cout << Params::codeVersion << "\n\n";
    cout << "Usage: hyb [-f hybrid.cfg] [-cont breakpoint.dat] [-paradyn] [-justinit] [-dump default.cfg]\n\n";
    cout << "-f: start a simulation run with a given config file\n\n";
    cout << "-cont: continue a simulation run from a breakpoint\n\n";
    cout << "-paradyn: produce a complete parameter log without actually running a simulation\n\n";
    cout << "-justinit: initialize a simulation run and exit\n\n";
    cout << "-dump: dump default config file variables (incomplete)\n\n";
}

//! Main program
int main(int argc, char *argv[])
{
    // Check command line arguments
    if (argc < 3) {
        showUsage();
        return -1;
    }
    bool configFileGiven = false;
    bool paramsDynamics = false;
    bool justInitialize = false;
    // Go through command line arguments
    for (int i = 1; i < argc; ++i) {
        // Config file
        if (strcmp(argv[i], "-f") == 0 && i+1 < argc) {
            Params::configFileName = argv[i+1];
            configFileGiven = true;
        }
        // Whole state dump file
        if (strcmp(argv[i], "-cont") == 0 && i+1 < argc) {
            Params::wsFileName = argv[i+1];
            Params::wsFileGiven = true;
        }
        // Default variable dumpping
        else if (strcmp(argv[i], "-dump") == 0 && i+1 < argc) {
            const char *dumpFileName = argv[i+1];

            // Check whether the dumpfile already exists
            fstream filetest;
            filetest.open(dumpFileName,ios::in);
            if(filetest.is_open() == true) {
                cerr << "ERROR [HybridMain]: given dumpfile already exists or write error\n";
                filetest.close();
                return -1;
            }
            // Dump variables
            cout << "Dumping default variables into configuration file " << dumpFileName  << ".. " << flush;
            simuConfig.init(true);
            simuConfig.dumpVars(dumpFileName);
            cout << "done!\n" << flush;
            return 0;
        }
        // Parameter dynamics
        else if (strcmp(argv[i], "-paradyn") == 0) {
            paramsDynamics = true;
            cout << "Creating parameter dynamics logfile.. " << flush;
        }
        // Just initialization
        else if (strcmp(argv[i], "-justinit") == 0) {
            justInitialize = true;
        }
    }
    // User should give the config file
    if (configFileGiven == false) {
        cerr << "ERROR [HybridMain]: no config file given\n";
        return -1;
    }
    // Check whether the config file exists
    fstream filetest;
    filetest.open(Params::configFileName,ios::in);
    if(filetest.is_open() == false) {
        cerr << "ERROR [HybridMain]: cannot find config file (" << Params::configFileName << ")\n";
        return -1;
    } else {
        filetest.close();
    }
    // Construct Simulation object named simu
    Simulation simu;
    // Construct parameter dynamics logfile
    if (paramsDynamics == true) {
        simu.runOnlyParameterDynamics();
        cout << "done!\n" << flush;
        return 0;
    }
    // Exit
    if(justInitialize == true) {
        cout << "Initilization done!\n";
        return 0;
    }
    // Run simulation
    simu.run();
    // Finalize simu
    simu.finalize();
    return 0;
}

