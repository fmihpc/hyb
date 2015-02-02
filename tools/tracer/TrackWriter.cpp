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

#include "TrackWriter.h"
#include "Token.h"
#include "constants.H"
#include <iostream>
#include <fstream>
#include <string.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <unistd.h>

using namespace std;


TrackWriter::TrackWriter(Config *c, int *gd, double *sb, double cw)
{
    cfg = c;
    cellWidth = cw;
    string num;

    fileNum = currentFileNum();
    if(fileNum == -1 or (not cfg->overwrite)) 
        fileNum++;

    if(cfg->overwrite) // unlink all files from last trace
        removeLastTrace();

    memcpy( gridDims, gd, sizeof(int)*3);
    memcpy( simBox, sb, sizeof(double)*6);
}


void splitFilename( const char *fname, string *tokens)
{
    tokens[0] = string(), tokens[1] = string(), tokens[2] = string(); tokens[3] = string();

    int i = 0;
    for(;fname[i] and fname[i]!='.';i++) 
        tokens[0].push_back(fname[i]);

    if(fname[i]) i++;
    for(;fname[i] and fname[i]!='.';i++) 
        tokens[1].push_back(fname[i]);

    if(fname[i]) i++;
    for(;fname[i] and fname[i]!='.';i++) 
        tokens[2].push_back(fname[i]);

    if(fname[i]) i++;
    for(;fname[i] and fname[i]!='.';i++) 
        tokens[3].push_back(fname[i]);
}


void TrackWriter::writeCfg(vector<Track> &plist, vector<Track> &plist_bw)
{
    string fn;
    outFilename(fn, "cfg");
    ofstream of(fn.c_str());
    if(!of.good())
    {cerr << "Error: Cannot open save file for config info " << fn << endl;return;}
    cfg->writeCfg(of);

    if(cfg->direction == "forward")
        for( vector<Track>::iterator it = plist.begin(); it<plist.end(); it++)
            of << it->x[0] << " " << it->y[0] << " " << it->z[0] << " " << it->v0[0] << " " << it->v0[1] << " " << it->v0[2] 
               << " " << (int)(it->mass/cnst::mp) << " " << (int)(it->charge/cnst::qe) <<endl; 
    else // bw
        for( vector<Track>::iterator it = plist_bw.begin(); it<plist_bw.end(); it++)
            of << it->x[0] << " " << it->y[0] << " " << it->z[0] << " " << it->v0[0] << " " << it->v0[1] << " " << it->v0[2] 
               << " " << (int)(it->mass/cnst::mp) << " " << (int)(it->charge/cnst::qe) <<endl; 
}


void TrackWriter::removeLastTrace()
{
    string fn;
    for(int i=0; formats[i]; i++)
    {
        outFilename(fn,formats[i]);
        unlink(fn.c_str());                
    }
    outFilename(fn,"cfg");
    unlink(fn.c_str());                
}


int TrackWriter::currentFileNum()
{
    int max_fnum = -1;
    string tokens[4];
    DIR *dir;
    struct dirent *d_ent;
    string dirname;

    if(not cfg->out_dir.empty())
        dirname = cfg->out_dir;
    dir = opendir( dirname.c_str());
    if(!dir)
    {
        cerr << "Error: cannot open output dir " << dirname.c_str() << endl;
        exit(-1);
    }

    while( (d_ent=readdir(dir)) )
    {
        splitFilename(d_ent->d_name, tokens);
        if(tokens[0] == cfg->cfg_fname and tokens[1] == "trace" 
           and atoi(tokens[2].c_str()) > max_fnum)  // and tokens[2] == ext)
            max_fnum = atoi(tokens[1].c_str());
    }
    closedir(dir);
    return max_fnum;
}


void TrackWriter::outFilename( string &fn, const char* ext)
{
    /* file name format: <config filename>.trace.number.extension.  Example:  mutli.cfg.trace.0.vtk */
    string fnum;
    itos(fileNum,fnum);
    fn = cfg->out_dir + "/" + cfg->cfg_fname + "_trace";
    fn += "_" + fnum + "." + ext;
}


// Function writes trace data in Point3D format which can be read by matlab and octave. 
void TrackWriter::writeASCII( vector<Track> &plist, vector<Track> &plist_bw, string ext)
{
    int n=0;
    bool fw = not (cfg->direction == "backward");
    bool bw = not (cfg->direction == "forward");

    string fn;

    outFilename(fn, ext.c_str());

    ofstream of(fn.c_str());
    if(!of.good())
    {
        cerr << "Error: Cannot open ascii dump file " << fn << endl;
        return;
    }

    if(ext=="m")
        of << "% ";
    of << "x y z ";
    string zeros;
    for(int k=0; k < cfg->tvi_len;k++) 
    {
        of << cfg->tvars_intpol[k] << " ";
        zeros += " 0.";
    }
    for(int k=0; k < cfg->tvo_len;k++) 
    {
        of << cfg->tvars_other[k] << " ";
        zeros += " 0.";
    }
    of << endl;

    of << simBox[0] << " " << simBox[2] << " " << simBox[4] << " " <<  zeros << endl;
    of << simBox[1] << " " << simBox[2] << " " << simBox[4] << " " <<  zeros << endl;
    of << simBox[0] << " " << simBox[3] << " " << simBox[4] << " " <<  zeros << endl;
    of << simBox[1] << " " << simBox[3] << " " << simBox[4] << " " <<  zeros << endl;
    of << simBox[0] << " " << simBox[2] << " " << simBox[5] << " " <<  zeros << endl;
    of << simBox[1] << " " << simBox[2] << " " << simBox[5] << " " <<  zeros << endl;
    of << simBox[0] << " " << simBox[3] << " " << simBox[5] << " " <<  zeros << endl;
    of << simBox[1] << " " << simBox[3] << " " << simBox[5] << " " <<  zeros << endl;

    if(fw)
    {
        for( vector<Track>::iterator it = plist.begin(); it<plist.end(); n++,it++)
            for(int i=0;i<it->used;i++)
            {
                of  << it->x[i]  << " " 
                    << it->y[i]  << " " 
                    << it->z[i]  << " ";
                int k;
                for(k=0; k < cfg->tvi_len;k++)
                    of << it->data[k][i] << " ";
                for(k=0; k < cfg->tvo_len;k++)
                    of << it->data[cfg->tvi_len+k][i] << " ";
                of << endl;
            }
    }
    if(bw)
    {
        for( vector<Track>::iterator it = plist_bw.begin(); it<plist_bw.end(); n++,it++)
            for(int i=0;i<it->used;i++)
            {
                of  << it->x[i]  << " " 
                    << it->y[i]  << " " 
                    << it->z[i]  << " ";
                int k;
                for(k=0; k < cfg->tvi_len;k++)
                    of << it->data[k][i] << " ";
                for(k=0; k < cfg->tvo_len;k++)
                    of << it->data[cfg->tvi_len+k][i] << " ";
                of << endl;
            }
    }
    of.flush();
    of.close();
}


void TrackWriter::writeVTKtrace( ofstream &of, vector<Track> &plist, vector<Track> &plist_bw)
{
    bool fw = not (cfg->direction == "backward");
    bool bw = not (cfg->direction == "forward");
    vector<Track>::iterator it; 

    int pts_num = 0;
    int cells_length = 0;

    if(fw)
    {
        for(it=plist.begin();it<plist.end();it++)
        {
            pts_num += it->used;        
            cells_length += it->used+1; 
        }
    }
    if(bw)
    {
        for(it=plist_bw.begin();it<plist_bw.end();it++)
        {
            pts_num += it->used;        
            cells_length += it->used+1; 
        }
    }

    of << "POINTS " << pts_num+8 << " double" << endl;    
    /* 8 first points are the grid corners. */
    of << simBox[0] << " " << simBox[2] << " " << simBox[4] << endl;
    of << simBox[1] << " " << simBox[2] << " " << simBox[4] << endl;
    of << simBox[0] << " " << simBox[3] << " " << simBox[4] << endl;
    of << simBox[1] << " " << simBox[3] << " " << simBox[4] << endl;
    of << simBox[0] << " " << simBox[2] << " " << simBox[5] << endl;
    of << simBox[1] << " " << simBox[2] << " " << simBox[5] << endl;
    of << simBox[0] << " " << simBox[3] << " " << simBox[5] << endl;
    of << simBox[1] << " " << simBox[3] << " " << simBox[5] << endl;

    if(fw)
    {
        for(it=plist.begin();it<plist.end();it++)
            for(int i=0;i<it->used;i++)
                of << it->x[i] << " " << it->y[i] << " " << it->z[i] << endl;
    }
    if(bw)
    {
        for(it=plist_bw.begin();it<plist_bw.end();it++)
            for(int i=0;i<it->used;i++)
                of << it->x[i] << " " << it->y[i] << " " << it->z[i] << endl;
    }
        

    of << endl;
    int plen = 0;
    if(fw) plen += plist.size();
    if(bw) plen += plist_bw.size();
    of << "CELLS " << plen+8 << " " <<  cells_length+16 << endl;
    for(int i=0; i<8; i++) // grid corners
        of << 1 << " " << i << endl;
    int n = 8;
    if(fw)
    {
        for(it=plist.begin();it<plist.end();it++)
        {
            of << it->used << " ";
            for(int i=0;i<it->used;i++)
                of << n++ << " ";  
            of << endl;
        }
    }
    if(bw)
    {
        for(it=plist_bw.begin();it<plist_bw.end();it++)
        {
            of << it->used << " ";
            for(int i=0;i<it->used;i++)
                of << n++ << " ";  
            of << endl;
        }
    }
    of << endl;

    of << "CELL_TYPES " << plen+8 << endl;
    for(int i=0; i<8; i++) // grid corners
        of << 1 << endl;
    if(fw)
    {
        for(it=plist.begin();it<plist.end();it++)
            of << 4 << endl;
    }
    if(bw)
    {
        for(it=plist_bw.begin();it<plist_bw.end();it++)
            of << 4 << endl;
    }
    
    /* write data */
    of << "POINT_DATA " << pts_num+8 << endl;
    for(int j=0; j < cfg->tvi_len; j++)
    {        
        of << endl;
        of << "SCALARS  "  << cfg->tvars_intpol[j] << " double 1" << endl;
        of << "LOOKUP_TABLE default" << endl;
        for(int i=0; i<8; i++) // grid corners
            of << 0 << endl;
        if(fw)
        {
            for(it=plist.begin();it<plist.end();it++)
                for(int i=0;i<it->used;i++)
                    of << it->data[j][i] << endl; 
        }
        if(bw)
        {
            for(it=plist_bw.begin();it<plist_bw.end();it++)
                for(int i=0;i<it->used;i++)
                    of << it->data[j][i] << endl; 
        }
    }

    for(int j=0; j < cfg->tvo_len; j++)
    {        
        of << endl;
        of << "SCALARS  "  << cfg->tvars_other[j] << " double 1" << endl;
        of << "LOOKUP_TABLE default" << endl;
        for(int i=0; i<8; i++) // grid corners
            of << 0 << endl;
        if(fw)
        {
            for(it=plist.begin();it<plist.end();it++)
                for(int i=0;i<it->used;i++)
                    of << it->data[cfg->tvi_len+j][i] << endl; 
        }
        if(bw)
        {
            for(it=plist_bw.begin();it<plist_bw.end();it++)
                for(int i=0;i<it->used;i++)
                    of << it->data[cfg->tvi_len+j][i] << endl; 
        }
    }
}


void TrackWriter::writeVTK( vector<Track> &plist, vector<Track> &plist_bw)
{
    string fn;
    outFilename(fn,"vtk");

    ofstream of(fn.c_str());
    if(!of.good())
    {
        cerr << "Error: Cannot open VTK file " << fn << " for output." <<endl;
        return;
    }

    of << "# vtk DataFile Version 2.0" << endl 
       << "iontracer "  << cfg->version << " trace file." << endl
       << "ASCII" << endl
       << "DATASET UNSTRUCTURED_GRID" << endl;

    writeVTKtrace(of, plist, plist_bw);    
    of.close();
}


