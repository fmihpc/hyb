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

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <cstdlib>
#include <fstream>

typedef float shortreal; //!< Simulation shortreal type
typedef float fastreal; //!< Simulation fastreal type
typedef double real; //!< Simulation real type
typedef shortreal gridreal; //!< Simulation gridreal type
typedef real datareal; //!< Simulation datareal type
typedef int TPDF_ID; //!< Simulation TPDF_ID type

#define ERRORMSG(msg) errorlog << "ERROR [" << __FILE__ << "/" << __LINE__ << "]: " << msg << "\n";
#define ERRORMSG2(msgA,msgB) errorlog << "ERROR [" << __FILE__ << "/" << __LINE__ << "]: " << msgA << " (" << msgB << ")\n";
#define WARNINGMSG(msg) errorlog << "WARNING [" << __FILE__ << "/" << __LINE__ << "]: " << msg << "\n";
#define WARNINGMSG2(msgA,msgB) errorlog << "WARNING [" << __FILE__ << "/" << __LINE__ << "]: " << msgA << " (" << msgB << ")\n";
#define MSGFUNCTIONCALL(name) mainlog << "\n==== FUNCTION CALL [" << name << "]\n";
#define MSGFUNCTIONEND(name) mainlog << "==== FUNCTION END  [" << name << "]\n\n";

#define sign(a) (((a)<0)? -1: 1)
#define sqr(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define quad(x) ((x)*(x)*(x)*(x))
#define min2(a,b) (((a)<(b))?(a):(b))
#define min3(a,b,c) (min2(min2(a,b),c))
#define max2(a,b) (((a)>(b))?(a):(b))
#define max3(a,b,c) (max2(max2(a,b),c))
#define vecsqr(a) (sqr(a[0])+sqr(a[1])+sqr(a[2]))
#define normvec(a) (sqrt(vecsqr(a)))
#define norm(a,b,c) (sqrt(sqr(a) + sqr(b) +sqr(c)))

#include <cmath>
#define pi M_PI

#if (defined(__i386) || defined(__alpha)) && !defined(LITTLE_ENDIAN)
#  define LITTLE_ENDIAN
#endif

#ifdef LITTLE_ENDIAN
extern void ByteConversion(int sz, unsigned char *x, int n);
#else
inline void ByteConversion(int sz, unsigned char *x, int n) {}
#endif
#define ByteConversion_input ByteConversion

//! Write float type to ostream
inline void WriteFloatsToFile(std::ostream& o, const float xf[], int n)
{
    o.write((const char*)xf,sizeof(float)*n);
}

//! Write double type to ostream
inline void WriteDoublesToFile(std::ostream& o, const double x[], int n)
{
    o.write((const char*)x,sizeof(double)*n);
}

//! Read float type from istream
inline void ReadFloatsFromFile(std::istream& o, float xf[], int n)
{
    o.read((char*)xf,sizeof(float)*n);
}

//! Read double type from istream
inline void ReadDoublesFromFile(std::istream& o, double x[], int n)
{
    o.read((char*)x,sizeof(double)*n);
}

//! 3D vector cross product for real
inline void crossProduct(const real a[3], const real b[3], real result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

//! 3D vector cross product for gridreal
inline void crossProduct(const gridreal a[3], const gridreal b[3], gridreal result[3])
{
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

//! Dot product for real
inline real dotProduct(const real a[3], const real b[3])
{
    return  a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//! Dot product for gridreal
inline gridreal dotProduct(const gridreal a[3], const gridreal b[3])
{
    return  a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

//! Normalize real vector
inline void normalize(real a[3])
{
    real length = normvec(a);
    a[0]/= length;
    a[1]/= length;
    a[2]/= length;
}

//! Normalize gridreal vector
inline void normalize(gridreal a[3])
{
    gridreal length = normvec(a);
    a[0]/= length;
    a[1]/= length;
    a[2]/= length;
}

extern void doabort();
void TermHandler(int);
std::string int2string(int nn, int zeros);
real string2double(std::string str);
real getExecutionSecs();
real getLastIntervalSecs();
std::string secsToTimeStr(double);
std::string cropPrecedingAndTrailingSpaces(std::string str);
void skipBracedContent(std::ifstream &fileStream);

//! Split string into vector of strong according to the delim string
inline std::vector<std::string> splitString(std::string str, std::string delim)
{
    std::string::size_type cutAt;
    std::vector<std::string> result;
    while( (cutAt = str.find_first_of(delim)) != str.npos ) {
        if(cutAt > 0) {
            result.push_back(str.substr(0,cutAt));
        }
        str = str.substr(cutAt+1);
    }
    if(str.length() > 0) {
        result.push_back(str);
    }
    return result;
}

//! Convert a vector of string to a vector of reals
inline std::vector<real> stringListToReals(std::string str)
{
    std::vector<std::string> splitted = splitString(str," ");
    std::vector<real> result;
    for (unsigned int ii = 0; ii < splitted.size(); ++ii) {
        result.push_back(atof(splitted[ii].c_str()));
    }
    return result;
}

//! Read an ascii file and convert it in a two dimensional vector of reals
inline std::vector< std::vector<real> > readRealsFromFile(const char* fn)
{
    std::ifstream in;
    std::vector< std::vector<real> > result;
    std::vector<std::string> strs;
    // read file lines into strs
    in.open(fn);
    if(in.good() == false) {
        return result;
    }
    bool eof = in.eof();
    while(in.eof() == false) {
        char line[1024];
        in.getline(line,1024);
        eof = in.eof();
        if(eof == false) {
            strs.push_back(std::string(line));
        }
    }
    in.close();
    // convert lines into real vectors
    for (unsigned int i = 0; i < strs.size(); ++i) {
        if (strs[i].find_first_not_of(" \n\t") == std::string::npos) {
            continue;
        }
        result.push_back(stringListToReals(strs[i]));
    }
    return result;
}

#endif

