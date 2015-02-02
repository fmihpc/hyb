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

#include "Track.h"
#include "Token.h"
#include "Config.h"
#include "constants.H"
#include "PointReader.h"
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>


PointReader::PointReader( Config* c, double *box)
{
    cfg = c;
    npts = 1;
    memcpy( boxLimits, box, sizeof(double)*6);
    mass = cnst::mp;
    charge = cnst::qe;
    length_scale = 1.0;
    velocity[0] = 0.;
    velocity[1] = 0.;
    velocity[2] = 0.;
    if(cfg->endpts)
        track_len = 2;
    else
        track_len = cfg->maxsteps; 
}


inline double sign(double x) // note: sign(0) = 0
{ 
    return (x>0)-(x<0);
}


//NOTE: dont use with a=c, a=b or b=c...
inline void crossProd(double *a, double *b, double *c) 
{
    c[0] =  a[1]*b[2] - b[1]*a[2];
    c[1] = -a[0]*b[2] + b[0]*a[2];
    c[2] =  a[0]*b[1] - b[0]*a[1];
}


inline void dotProd(double *a, double *b, double *c) 
{
    c[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}


void PointReader::readPoints( const char* fn, bool iscfg, vector<Track> &plist, vector<Track> &plist_bw)
{
    ifstream pf(fn);
    if(!pf.good())
    {
        cerr << "Error: cannot open point data file '" << fn << "'." << endl;
        exit(1);
    }

    string line;
    vector<string> tokens;
    if(iscfg)  // skip to end of config data
    {
        while(getline(pf,line))
        {
            tokens.clear();
            tokenize(line,tokens);
            if(!tokens.empty() and tokens[0] == "EOC")
                break;
        }
    }

    Track *tr,*trb;
    while(getline(pf,line))
    {
        if(line.size()>1024)
        {
            cerr << "Error: line too long in point data file or config file" << endl;
            exit(-1);
        }

        tokens.clear();
        tokenize(line,tokens);

        if(tokens.empty() or tokens[0][0] == '#')
            continue;

        if(tokens[0] == "line")
        {
            // syntax: line point1 point2 num
            double bx,by,bz,ex,ey,ez;
            unsigned int num;

            if(tokens.size() != 8)
            {
                cerr << "Error: wrong number of arguments to 'line' directive. Edit point data." << endl;
                exit(-1);
            }

            getDouble(tokens[1], bx);
            getDouble(tokens[2], by);
            getDouble(tokens[3], bz);
            getDouble(tokens[4], ex);
            getDouble(tokens[5], ey);
            getDouble(tokens[6], ez);
            bx *= length_scale; by *= length_scale; bz *= length_scale;
            ex *= length_scale; ey *= length_scale; ez *= length_scale;

            /* Check Boundaries */

            if( bx<boxLimits[0] or bx>boxLimits[1] or ex<boxLimits[0] or ex>boxLimits[1] or 
                by<boxLimits[2] or by>boxLimits[3] or ey<boxLimits[2] or ey>boxLimits[3] or
                bz<boxLimits[4] or bz>boxLimits[5] or ez<boxLimits[4] or ez>boxLimits[5])
            {
                cerr << "Error: initial points of line directive are outside simulation box." << endl
                     << "Check points and length_scale" << endl;
                exit(-1);
            }

            num = atoi(tokens[7].c_str());
            if(num <= 1)
            {
                cerr << "Error: the number of points for line directive must be reasonable (>1). Edit point data." << endl;
                exit(-1);
            }
            double dx = sign(ex-bx)* fabs(bx-ex)/(double)(num-1);
            double dy = sign(ey-by)* fabs(by-ey)/(double)(num-1);
            double dz = sign(ez-bz)* fabs(bz-ez)/(double)(num-1);

            for(int n=0; n<(int)num; n++)
            {
                tr = new Track( bx+n*dx, by+n*dy, bz+n*dz, 
                                velocity[0], velocity[1], velocity[2],
                                mass, charge, track_len,
                                cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                trb = new Track( bx+n*dx, by+n*dy, bz+n*dz,
                                 velocity[0], velocity[1], velocity[2],
                                 mass, charge, track_len,        
                                 cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                plist.push_back(*tr);        
                plist_bw.push_back(*trb);        
                delete tr;
                delete trb;
            }
        }
        else if(tokens[0] == "spherical")
        {
            // A point in spherical coordinates.
            // syntax: spherical r theta phi 
            double r,phi,theta;
            double x,y,z;

            if(tokens.size() != 4)
            {
                cerr << "Error: wrong number of arguments to 'spherical' directive. Edit point data." << endl;
                exit(-1);
            }
            getDouble(tokens[1], r);
            getDouble(tokens[2], theta);
            getDouble(tokens[3], phi);

            if(r <= 0)
            {
                cerr << "Error: negative or zero radius in spherical point" << endl;
                exit(-1);
            }
            
            phi   = phi*3.14159265/180.;  // convert to radians
            theta = theta*3.14159265/180.; // convert to radians 
            r *= length_scale;

            x = r*sin(theta)*cos(phi);
            y = r*sin(theta)*sin(phi);
            z = r*cos(theta);

            /* check that point is in sim box */
            if( x<boxLimits[0] or x>boxLimits[1] or 
                y<boxLimits[2] or y>boxLimits[3] or 
                z<boxLimits[4] or z>boxLimits[5]) 
            {
                cerr << "Error: spherical point goes outside of simulation box." << endl
                     << "Check spherical directive and length_scale" << endl;
                exit(-1);
            }
            tr = new Track( x, y, z,
                            velocity[0], velocity[1], velocity[2],
                            mass, charge, track_len,
                            cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
            trb = new Track( x, y, z,
                             velocity[0], velocity[1], velocity[2],
                             mass, charge, track_len,        
                             cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
            plist.push_back(*tr);        
            plist_bw.push_back(*trb);        
            delete tr;
            delete trb;
        }
        else if(tokens[0] == "randsphere") // syntax: randsphere midpoint radius radial_velocity number_of_particles
        {
            double ox,oy,oz,radius,rv;
            long int num;

            if(tokens.size() != 7)
            {
                cerr << "Error: wrong number of arguments to 'randsphere' directive. Edit point data." << endl;
                exit(-1);
            }
            getDouble(tokens[1], ox);
            getDouble(tokens[2], oy);
            getDouble(tokens[3], oz);
            getDouble(tokens[4], radius);
            getDouble(tokens[5], rv);
            getLong(tokens[6], num);
            radius *= length_scale;
            ox *= length_scale;
            oy *= length_scale;
            oz *= length_scale;

            if(radius <= 0)
                cerr << "Error: negative or zero radius for randsphere" << endl, exit(-1);

            int i=0; double dir[3]; double norm; double rp[3]; double v[3];
            while(i<num)
            {
                dir[0] = ((double) rand())/((double) RAND_MAX) - .5;
                dir[1] = ((double) rand())/((double) RAND_MAX) - .5;
                dir[2] = ((double) rand())/((double) RAND_MAX) - .5;
                norm = sqrt(dir[0]*dir[0]+dir[1]*dir[1]+dir[2]*dir[2]);
                if(norm > 1.) continue; // accept only vectors inside sphere to get even distribution.

                rp[0] = radius*dir[0]/norm + ox; // particle_position = direction_vector + sphere_origin.
                rp[1] = radius*dir[1]/norm + oy;
                rp[2] = radius*dir[2]/norm + oz;

                if( rp[0]<boxLimits[0] or rp[0]>boxLimits[1] or 
                    rp[1]<boxLimits[2] or rp[1]>boxLimits[3] or 
                    rp[2]<boxLimits[4] or rp[2]>boxLimits[5]) 
                    cerr << "Error: randsphere goes outside of simulation box." << endl
                         << "Check randsphere directive and length_scale" << endl, exit(-1);

                v[0] = rv*dir[0];   // velocity = radial_velocity*direction
                v[1] = rv*dir[1];
                v[2] = rv*dir[2];

                tr = new Track( rp[0], rp[1], rp[2], v[0], v[1], v[2],
                                mass, charge, track_len,
                                cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                trb = new Track( rp[0], rp[1], rp[2], v[0], v[1], v[2],
                                 mass, charge, track_len,        
                                 cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                plist.push_back(*tr);        
                plist_bw.push_back(*trb);        
                delete tr;
                delete trb;
                i++;
            }
        }
        else if(tokens[0] == "circle")
        {
            // syntax: circle midpoint radius normal_phi normal_theta (spherical coordinates) num
            double x,y,z,radius,phi,theta;
            long int num;

            if(tokens.size() != 8)
            {
                cerr << "Error: wrong number of arguments to 'circle' directive. Edit point data." << endl;
                exit(-1);
            }

            getDouble(tokens[1], x);
            getDouble(tokens[2], y);
            getDouble(tokens[3], z);
            getDouble(tokens[4], radius);
            getDouble(tokens[5], phi);
            getDouble(tokens[6], theta);
            getLong(tokens[7], num);

            if(radius <= 0)
            {
                cerr << "Error: negative or zero radius for circle" << endl;
                exit(-1);
            }
            
            phi   = phi*3.14159265/180.;  // convert to radians
            theta = theta*3.14159265/180.; // convert to radians 
            radius *= length_scale;
            x *= length_scale;
            y *= length_scale;
            z *= length_scale;
    
            /* new rotated basis */
            double n[3], a[3], b[3];

            n[0] = sin(theta)*cos(phi); // n normal in cartesian coordinates.
            n[1] = sin(theta)*sin(phi);
            n[2] = cos(theta);


            if(n[0]==0 and n[1]==0) // special cases
            {
                a[0]=0.; a[1]=1.; a[2]=0.;
                b[0]=1.; b[1]=0.; b[2]=0.;
            }
            else if( n[0]==0 and n[2]==0)
            {
                a[0]=1.; a[1]=0.; a[2]=0.;
                b[0]=0.; b[1]=0.; b[2]=1.;
            }
            else if(n[1]==0 and n[2]==0)
            {
                a[0]=0.; a[1]=1.; a[2]=0.;
                b[0]=0.; b[1]=0.; b[2]=1.;
            }
            else
            {
                a[0] = -n[1]*n[2];  // define a as: a dot_product n = 0 
                a[1] = -n[0]*n[2]; 
                a[2] = 2*n[1]*n[0];

                double norm = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
                a[0]/=norm; a[1]/=norm; a[2]/=norm;

                crossProd(a,n,b);   //  a x n = b
                norm = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
                b[0]/=norm; b[1]/=norm; b[2]/=norm;
            }

            double M[3][3];   // the rotation matrix
            M[0][0]=n[0]; M[1][0]=n[1]; M[2][0]=n[2];   // col 1 = n
            M[0][1]=a[0]; M[1][1]=a[1]; M[2][1]=a[2];   // col 2 = a
            M[0][2]=b[0]; M[1][2]=b[1]; M[2][2]=b[2];   // col 3 = b
            
            double dt = 2*3.14159265/((double)num);
            double ang = 0.;
            for(int i=0;i<num;i++)
            {
                double p[3];
                double rp[3];
                p[0] = 0.;
                p[1] = radius*cos(ang);
                p[2] = radius*sin(ang);

                /* Rotation */
                dotProd(M[0],p,&rp[0]);   // this is p=M*p, and it rotates the circle.
                dotProd(M[1],p,&rp[1]);
                dotProd(M[2],p,&rp[2]);                

                /* Translation */
                rp[0] += x;
                rp[1] += y;
                rp[2] += z;

                /* check that point is in sim box */
                if( rp[0]<boxLimits[0] or rp[0]>boxLimits[1] or 
                    rp[1]<boxLimits[2] or rp[1]>boxLimits[3] or 
                    rp[2]<boxLimits[4] or rp[2]>boxLimits[5]) 
                {
                    cerr << "Error: cirlce goes outside of simulation box." << endl
                         << "Check circle directive and length_scale" << endl;
                    exit(-1);
                }
                tr = new Track( rp[0], rp[1], rp[2],
                                velocity[0], velocity[1], velocity[2],
                                mass, charge, track_len,
                                cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                trb = new Track( rp[0], rp[1], rp[2],
                                 velocity[0], velocity[1], velocity[2],
                                 mass, charge, track_len,        
                                 cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
                plist.push_back(*tr);        
                plist_bw.push_back(*trb);        
                delete tr;
                delete trb;
                ang += dt;
            }
        }
        else if(tokens[0] == "multipoint")
        {
            // multipoint:  x, y, z, radial_velocity, 
            cerr << "Error: multipoint is not implemented yet." << endl,
            exit(-1);

        }
        else if(tokens[0] == "set")
        {
            if(tokens[1] == "charge")
            {
                double new_q;
                getDouble(tokens[2],new_q);
                charge = new_q*cnst::qe;
            }
            else if(tokens[1] == "mass")
            {
                double new_m;
                getDouble(tokens[2],new_m);
                mass = new_m*cnst::mp;
            }
            else if(tokens[1] == "length_scale")
            {
                double new_l;
                getDouble(tokens[2],new_l);

                if(new_l <= 0)
                {
                    cerr << "Error: negative or zero scale length" << endl;
                    exit(-1);
                }
                length_scale = new_l;
            }
            else if(tokens[1] == "velocity")
            {
                getDouble(tokens[2],velocity[0]);
                getDouble(tokens[3],velocity[1]);
                getDouble(tokens[4],velocity[2]);
            }
            else
            {
                cerr << "Error: undefined parameter for set: " << tokens[1] << endl;
                exit(-1);
            }
        }
        else if(isDouble(tokens[0]))  // try point data
        {
            // syntax distinct point data:  
            //      1. explicit: x y z vx vy vz mass charge
            //      2. or just:  x y z  (using the values set by directives)

            double x, y, z, vx, vy, vz, m, q;
            getDouble(tokens[0], x);
            getDouble(tokens[1], y);
            getDouble(tokens[2], z);
            x *= length_scale;
            y *= length_scale;
            z *= length_scale;

            if(tokens.size()==8)
            {
                getDouble(tokens[3], vx);
                getDouble(tokens[4], vy);
                getDouble(tokens[5], vz);
                getDouble(tokens[6], m);
                getDouble(tokens[7], q);
            }
            else if(tokens.size()==3)
            {
                vx = velocity[0]; 
                vy = velocity[1]; 
                vz = velocity[2]; 
                m = mass;
                q = charge;
            }
            else
            {
                cerr << "Error: wrong number of tokens in point data (correct value 3 or 8)" << endl
                     << "Check points in the point data or config file" << endl;
                exit(-1);
            }

            /* check limits */
            if( x<boxLimits[0] or x>boxLimits[1] or
                y<boxLimits[2] or y>boxLimits[3] or
                z<boxLimits[4] or z>boxLimits[5])
            {
                cerr << "Error: Some points go outside of the simulation box." << endl
                     << "Check points in the point data or config file" << endl;
                exit(-1);
            }

            if( m <= 0 or q== 0)
            {
                cerr << "Error: mass <= 0 or charge is zero" << endl;
                exit(-1);
            }

            tr = new Track(x,y,z,vx,vy,vz, cnst::mp*m, cnst::qe*q, track_len,
                           cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
            trb = new Track(x,y,z,vx,vy,vz, cnst::mp*m, cnst::qe*q, track_len,        
                            cfg->tvars_intpol, cfg->tvi_len, cfg->tvars_other,cfg->tvo_len,npts++);        
            plist.push_back(*tr);        
            plist_bw.push_back(*trb);        
            delete tr;
            delete trb;
        }
        else
        {
            cerr << "Error: syntax error in point data" << endl;
            exit(-1);
        }
    }
    pf.close();
}

