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

#ifndef CONFIG_H 
#include <vector>
#include <string>
#include <map>
#include <limits>
#include <stdlib.h>

using namespace std;

#define VERSION "0.2.1"
#define ORBIT_DATA_MAX 10
#define MAX_HCF 8

extern const char* formats[]; 
enum {BUN_EK,BUN_VX,BUN_VY,BUN_VZ,PARID}; // variable ids for variables that are hcintpol variables

class Config{
    private:

        void fill_trace_vars2(string );
        void fill_file_formats(string line);
        void fill_hcfile(string line);
        bool check_other_var(string );
        bool format_ok(string value);
    
    public:
        /* Configuration variables. */
        string pressuretermfile;
        vector<string> hcfiles;  
        vector<double> hcf_mass; // relative mass.
        vector<double> hcf_charge;
        vector<string> file_formats; 
        int maxsteps;
        double stepsize_p;       // stepsize for particle trace.
        string tracetype;

        /* Arrays for trace variables. There are two kinds those that are interpolated and those that are not. */
        char *tvars_intpol[ORBIT_DATA_MAX];  // trace vars that are interpolated. 
        char *tvars_other[ORBIT_DATA_MAX];   // other var names. E.g. kinetic energy. 
        int tvo_flags[ORBIT_DATA_MAX];       // flags other vars.
        map<string,int> tvo2id;              // name to id map for other vars.
        string bunemanversion;
        string direction;
        bool verbose;
        bool overwrite;
        bool endpts;
        int intpolorder;
        string out_dir;
        string version;

        double xmin; // these are for box boundaries other than original simulation box
        double ymin;
        double zmin;
        double xmax;
        double ymax;
        double zmax;
        double planetary_boundary; // NOTE: square of planetary boundary. 

        /* Helper variables */
        int tvi_len; 
        int tvo_len; 
        string cfg_fname;

        Config();
        void readCfg(const char*);
        void writeDefaultCfg();
        void writeCfg(ostream &of);
        ~Config();
};

#define CONFIG_H 
#endif
