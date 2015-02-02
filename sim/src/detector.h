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

#ifndef DETECTOR_H
#define DETECTOR_H

#include <cstdio>
#include <vector>
#include <fstream>
#include "definitions.h"
#include "particle.h"
#include "params.h"
#include "vectors.h"

#define MAX_PART_DETECTORS 30

//! All detector arguments
struct DetectorArgs {
    //TEMPLATE DETECTOR VARIABLES
    stringArg popIdStr;
    stringArg detectionFile;
    realArg2 detectionTime;
    realArg maxCounts;
    stringArg coordinateFile;
    stringArg testParticleFile;
    realArg m;
    realArg q;
    functionArg2 detectorFUNC;
    DetectorArgs();
    void clearArgs();
};

//! Field point detector
class FieldDetect
{
protected:
    std::ofstream *fs; //!< Detector output file
    gridreal r[3]; //!< Detector coordinates
public:
    FieldDetect(std::ofstream *fs1, const Tgr3v point1);
    ~FieldDetect();
    void run(void);
};

//! Particle detector base class
class PartDetect
{
protected:
    std::ofstream *fs; //!< Detector output file
    gridreal R;
    gridreal Radius2;
    std::string partDetectType; //!< Detector type
public:
    PartDetect();
    PartDetect(std::ofstream *fs1, std::vector<real> partDetectArgs);
    virtual void firstline(void);
    inline int run(const TLinkedParticle* part, const gridreal r_new[3]);
    inline void save(const TLinkedParticle* part, const gridreal r_new[3]);
    virtual bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    virtual bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
    virtual ~PartDetect();
};

//! Test particle for tracing
class TestParticle
{
public:
    gridreal x,y,z,vx,vy,vz; //!< Position and velocity
    real q,m; //!< Charge and mass
    bool propagate; //!< Propagate or not
    TestParticle(gridreal r[3], gridreal v[3], real charge, real mass);
    inline bool run(void);
private:
    inline bool boundaries(void);
    void PropagateV(void);
};

//! Particle and field detectors
class Detector
{
public:
    Detector();
    Detector(DetectorArgs args);
    std::string popIdStr;
    real maxCounts;
    std::string detectorType;
    virtual ~Detector();
    void runFieldDetects();
    void runPartDetects(const TLinkedParticle* part, const gridreal r_new[3]);
    void runTestParticles();
    std::string toString();
    std::string configDump();
    unsigned int getNumberOfDetects();
    unsigned int getNumberOfFieldDetects();
    unsigned int getNumberOfPartDetects();
    std::string detectionFile;
protected:
    std::string coordinateFile, testParticleFile;
    std::string detectionFile1st, detectionFile2nd;
    real detectionTime[2];
    real mass, charge;
    bool stillPropagating;
    std::ofstream *files;
    std::vector<FieldDetect*> fieldDetects;
    std::vector<PartDetect*> partDetects;
    std::vector<TestParticle*> testParts;
    std::vector<std::string> detectorFunctionNames;
    std::vector<std::vector<real> > detectorFuncArgs;
    std::vector<std::string> detectionFiles;
    std::vector<real> currentCounts;
};

//! Field detector set
class FieldDetectorSet : public Detector
{
public:
    FieldDetectorSet(DetectorArgs args);
    unsigned int fieldDetectsFromFile;
    unsigned int fieldDetectsFromFunc;
private:
    std::vector<Tgr3v> pointCoordinates; //point detector coordinates
    void readCoordinateFile(const std::string coordinateFile, std::vector<Tgr3v>& coordinateVector);
};

//! Particle detector set
class ParticleDetectorSet : public Detector
{
public:
    const unsigned int numberOfPartDetectTypes;
    ParticleDetectorSet(DetectorArgs args);
    real getCurrentCount(int detectId);
private:
    static std::vector<std::string> partDetectTypes;
    static bool setupPartFunc;
    static std::vector<PartDetect* (*) (std::ofstream*,std::vector<real>)> newPartDetectFuncs;
};

//! Test particle set
class TestParticleSet : public Detector
{
public:
    TestParticleSet(DetectorArgs args);
private:
    std::vector<bool> propagating;
};

//! Class to create detector objects
class DetectorFactory
{
public:
    DetectorFactory();
    ~DetectorFactory();
    static Detector* createDetector(const std::string detectorType, DetectorArgs args);
    static int getNumberOfDetectors();
    static int getNumberOfDetectors(const std::string detectorType);
private:
    static unsigned int numberOfDetectors[3]; //field, particle, and testparticle
};

//! Spherical particle detector
class PartDetect_Sphere : public PartDetect
{
protected:
    gridreal r[3];
public:
    PartDetect_Sphere(std::ofstream *fs1, std::vector<real> partDetectArgs);
    virtual void firstline(void);
    virtual bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    virtual bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Spherical particle detector
class PartDetect_SphereInto : public PartDetect_Sphere
{
public:
    PartDetect_SphereInto(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Spherical particle detector
class PartDetect_SphereInside : public PartDetect_Sphere
{
public:
    PartDetect_SphereInside(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Spherical particle detector
class PartDetect_SphereOut : public PartDetect_Sphere
{
public:
    PartDetect_SphereOut(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Planar particle detector
class PartDetect_XPlane : public PartDetect
{
protected:
    gridreal x_plane;
public:
    PartDetect_XPlane(std::ofstream *fs1, std::vector<real> partDetectArgs);
    virtual void firstline(void);
    virtual inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    virtual inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Planar particle detector
class PartDetect_XPlaneReverse : public PartDetect_XPlane
{
public:
    PartDetect_XPlaneReverse(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Line particle detector
class PartDetect_Line : public PartDetect
{
protected:
    unsigned int Npoints;
    std::vector<Tgr3v> Points;
    std::vector<Tgr3v> directions;
    gridreal mincoord[3], maxcoord[3];
public:
    PartDetect_Line(std::ofstream *fs1, std::vector<real> partDetectArgs);
    virtual void firstline(void);
    virtual inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    virtual inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Line particle detector
class PartDetect_LineInside : public PartDetect_Line
{
public:
    PartDetect_LineInside(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Ellipse particle detector
class PartDetect_Ellipse : public PartDetect
{
protected:
    gridreal a,b,e_a[3],e_b[3],e_c[3], ecc,ae,a_ecc;
    bool circle;
public:
    virtual void firstline(void);
    PartDetect_Ellipse(std::ofstream *fs1, std::vector<real> partDetectArgs);
    virtual bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    virtual inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

//! Ellipse particle detector
class PartDetect_EllipseInside : public PartDetect_Ellipse
{
public:
    PartDetect_EllipseInside(std::ofstream *fs1, std::vector<real> partDetectArgs);
    void firstline(void);
    inline bool InsideDetector(const gridreal x, const gridreal y, const gridreal z);
    inline bool InsideDetector2(const gridreal x, const gridreal y, const gridreal z);
};

#endif

