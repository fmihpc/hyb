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

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <iostream>
#include <string>
#include "detector.h"
#include "simulation.h"
#include "params.h"
#include "magneticfield.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "transformations.h"
#endif

using namespace std;

extern Tgrid g;

//! Constructor for detector argument struct
DetectorArgs::DetectorArgs()
{
    clearArgs();
}

//! Clear detector arguments inside DetectorArgs()
void DetectorArgs::clearArgs()
{
    // TEMPLATE DETECTOR VARIABLES
    popIdStr.value.clear();
    popIdStr.given = false;
    detectionFile.value.clear();
    detectionFile.given = false;
    detectionTime.value[0] = 0;
    detectionTime.value[1] = 0;
    detectionTime.given = false;
    maxCounts.value = 0;
    maxCounts.given = false;
    coordinateFile.value.clear();
    coordinateFile.given = false;
    testParticleFile.value.clear();
    testParticleFile.given = false;
    m.value = 0;
    m.given = false;
    q.value = 0;
    q.given = false;
    detectorFUNC.name.clear();
    detectorFUNC.funcArgs.clear();
    detectorFUNC.given = false;
}

//! Constructor to create detector object
Detector::Detector(DetectorArgs args)
{
    if (Params::DETECTORS >= Params::MAX_DETECTORS) {
        ERRORMSG("too many detectors");
        doabort();
    }
    //check that necessary information is given
    if (args.detectionFile.given == false) {
        ERRORMSG ("detector result file must be given");
        doabort();
    }
    if (args.detectionTime.given == false) {
        ERRORMSG ("detection begin and end times must be given");
        doabort();
    }
    detectionFile = args.detectionFile.value;
    const string::size_type pdot = detectionFile.find('.');
    const string::size_type dfl = detectionFile.size();
    detectionFile1st = detectionFile;
    detectionFile2nd.clear();
    if (pdot!=string::npos) {
        detectionFile1st.erase(pdot,dfl);
        detectionFile2nd = detectionFile.substr(pdot+1,dfl-pdot-1);
    }
    const string detFile_reassembled
        = detectionFile1st+((pdot==string::npos)?"":"."+detectionFile2nd);
    if (args.detectionFile.value.compare(detFile_reassembled) != 0) {
        ERRORMSG2 ("deviding detectionFile string unsuccessful", detFile_reassembled);
        doabort();
    }
    detectionTime[0] = args.detectionTime.value[0];
    detectionTime[1] = args.detectionTime.value[1];
    fieldDetects.clear();
    partDetects.clear();
    testParts.clear();
    detectionFiles.clear();
    currentCounts.clear();
    if (args.detectorFUNC.given == true) {
        for (unsigned int i=0; i<args.detectorFUNC.name.size(); i++) {
            detectorFunctionNames.push_back(args.detectorFUNC.name[i]);
            detectorFuncArgs.push_back(args.detectorFUNC.funcArgs[i]);
        }
    }
}

//! Dummy constructor
Detector::Detector() : detectorType("Invalid") { }

//! Dummy virtual destructor
Detector::~Detector() { }

//! String representation for a detector
string Detector::toString()
{
    stringstream ss;
    ss << "detectorType = " << detectorType << "\n";
    if (detectorType.compare("field")==0) {
        ss << "Number of point detectors for the fields: "
           << fieldDetects.size() << "\n";
        ss << "Detection time:  from " << detectionTime[0] << " s to "
           << detectionTime[1] << " s\n";
        ss << "Coordinate file: " << coordinateFile << "\n";
        ss << "Detection files: \n";
        for (unsigned int i=0; i<fieldDetects.size(); i++) {
            ss << detectionFiles[i] << "\n";
        }
    }
    if (detectorType.compare("particle")==0) {
        if (popIdStr.compare("-") == 0) {
            ss << "Collecting all populations.\n";
        } else {
            ss << "Collecting " << popIdStr << " population.\n";
        }
        ss << "Number of particle detectors: " << partDetects.size() << "\n";
        ss << "Detection time:  from " << detectionTime[0] << " s to "
           << detectionTime[1] << " s\n";
        ss << "Maximum number of counts: " << maxCounts << "\n";
        ss << "Detection files: \n";
        for (unsigned int i=0; i<partDetects.size(); i++) {
            ss << detectionFiles[i] << "\n";
        }
    }
    if (detectorType.compare("testparticle")==0) {
        ss << "Number of test particles: " << testParts.size() << "\n";
        ss << "Insertion time: " << detectionTime[0] << "\n";
        ss << "Test particle propagation ends at "<< detectionTime[1] << " s\n";
        ss << "Charge of test particles: " << charge/Params::e << " e\n";
        ss << "Mass of test particles: " << mass
           << " (= " << mass/Params::amu << " amu)";
        ss << "Detection file: " << detectionFile << "\n";
    }
    return ss.str();
}

//! Returns detector example for config file (not implemented)
string Detector::configDump()
{
    return "";
}

//! Number of particle+field detectors
unsigned int Detector::getNumberOfDetects()
{
    return fieldDetects.size() + partDetects.size();
}

//! Number of field detectors
unsigned int Detector::getNumberOfFieldDetects()
{
    return fieldDetects.size();
}

//! Number of particle detectors
unsigned int Detector::getNumberOfPartDetects()
{
    return partDetects.size();
}

//! Number of field and particle detectors and test particles
unsigned int DetectorFactory::numberOfDetectors[3] = {0, 0, 0};

// Construction functions for detector factory
Detector* newFieldDetector(DetectorArgs args)
{
    return new FieldDetectorSet(args);
}
Detector* newParticleDetector(DetectorArgs args)
{
    return new ParticleDetectorSet(args);
}
Detector* newTestParticleSet(DetectorArgs args)
{
    return new TestParticleSet(args);
}

//! Constructor
DetectorFactory::DetectorFactory() { }

//! Destructor
DetectorFactory::~DetectorFactory() { }

//! Creates new detector objects
Detector* DetectorFactory::createDetector(const string detectorType, DetectorArgs args)
{
    mainlog << "CREATING A DETECTOR [DetectorFactory::createDetector]: " << args.detectionFile.value << "\n";
    Detector* temp = NULL;
    //check if same detectionFile string is already given for another detector
    for (unsigned int i=0; i<Params::detectors.size(); i++) {
        if (args.detectionFile.value.compare(Params::detectors[i]->detectionFile)==0) {
            ERRORMSG2("detectionFile cannot be same for multiple detector sets",
                      args.detectionFile.value);
            doabort();
            return temp;
        }
    }
    if (detectorType.compare("field") == 0) {
        temp = newFieldDetector(args);
        numberOfDetectors[0]++;
    } else if (detectorType.compare("particle") == 0) {
        temp = newParticleDetector(args);
        numberOfDetectors[1]++;
    } else if (detectorType.compare("testparticle") == 0) {
        temp = newTestParticleSet(args);
        numberOfDetectors[2]++;
    } else {
        ERRORMSG2("cannot find detector type",detectorType);
        doabort();
    }
    return temp;
}

//! Return total number of detectors
int DetectorFactory::getNumberOfDetectors()
{
    return numberOfDetectors[0]+numberOfDetectors[1];
}

//! Return number of detectors of a certain type
int DetectorFactory::getNumberOfDetectors(const string detectorType)
{
    if (detectorType.compare("field") == 0) {
        return numberOfDetectors[0];
    } else if (detectorType.compare("particle") == 0) {
        return numberOfDetectors[1];
    } else if (detectorType.compare("testparticle") == 0) {
        return numberOfDetectors[1];
    } else {
        ERRORMSG2("cannot find detector type",detectorType);
        doabort();
        return 0;
    }
}

//! Constructor
FieldDetect::FieldDetect(ofstream *fs1, Tgr3v point1)
{
    for (int i=0; i<3; i++) {
        r[i] = point1[i];
    }
    fs = fs1;
    fs->precision(16);
    // output file format:
    //     t  n1 vx1 vy1 vz1 P1  n2 vx2 vy2 vz2 P2 ... uex uey uez jx jy jz Bx By Bz
    // where (n,v,P) groups are given for each population
    (*fs) << "% t ";
    for (int popi=0; popi<Params::POPULATIONS; popi++) {
        (*fs) << " n" << popi+1 << " vx" << popi+1 << " vy" << popi+1
              << " vz" << popi+1 << " P" << popi+1;
    }
    (*fs) << " uex uey uez jx jy jz Bx By Bz\n% Populations:  ";
    for (int popi=0; popi<Params::POPULATIONS; popi++) {
        string temp = Params::pops[popi]->getIdStr();
        (*fs) << temp << " ";
    }
    (*fs) << "\n% Position: x=" << r[0] << ", y=" << r[1] << ", z=" << r[2] << " [m]\n"
          "% Field values are in SI base units: time t in s, n in m-3, v in m/s, P in Pa, etc.\n"
          << flush;
}

//! Destructor
FieldDetect::~FieldDetect()
{
    fs->close();
}

//! Run field detector
void FieldDetect::run()
{
    (*fs) << Params::t << " ";
    real numberdens = 0, vx = 0, vy = 0, vz = 0, P = 0;
    vector <int> popId;
    popId.clear();
    popId.push_back(0);
    for (int popi=0; popi<Params::POPULATIONS; popi++) {
        // Output n,vx,vy,vz,P
        popId[0]=popi;
        g.cellintpol_fluid(r,numberdens,vx,vy,vz,P, popId);
        (*fs) << numberdens << " " << vx << " "
              << vy << " " << vz << " " << P << " ";
    }
    real ue[3],j[3],B[3];
    g.cellintpol(Tgrid::CELLDATA_UE,ue);
    g.cellintpol(Tgrid::CELLDATA_J,j);
    g.cellintpol(Tgrid::CELLDATA_B,B);
    (*fs) << ue[0] << " " << ue[1] << " " << ue[2] << " "
          << j[0] << " " << j[1] << " " << j[2] << " "
          << B[0] << " " << B[1] << " " << B[2] << "\n" << flush;
}

//! Dummy Constructor for a base type particle detect
PartDetect::PartDetect(void) { }

//! Constructor for a base type particle detect
PartDetect::PartDetect(ofstream *fs1, vector <real> partDetectArgs)
{
    fs = fs1;
    fs->precision(5);
    (*fs) << scientific;
    R = partDetectArgs[0];
    Radius2 = sqr(R);
}

//! Dummy firstline function
void PartDetect::firstline(void)
{
    ERRORMSG ("PartDetect::firstline - virtual dummy");
    doabort();
}

//! Destructor for particle detect
PartDetect::~PartDetect()
{
    fs->close();
}

//! check whether to record (if r_new is InsideDetector and wasn't there before)
inline int PartDetect::run(const TLinkedParticle* part, const gridreal r_new[3])
{
    if (this->InsideDetector(r_new[0],r_new[1],r_new[2]) == false) {
        return 0;
    }
    if (this->InsideDetector2(part->x,part->y,part->z) == true) {
        return 0;
    }
    // now we know that old point was outside while current point is inside: a hit
    save(part, r_new);
    return 1;
}

//! has the geometry of each PartDetect
bool PartDetect::InsideDetector(gridreal x, gridreal y, gridreal z)
{
    ERRORMSG ("PartDetect::InsideDetector - virtual dummy");
    doabort();
    return false;
}

//! only plane detector should have InsideDetector2 different from InsideDetector
bool PartDetect::InsideDetector2(gridreal x, gridreal y, gridreal z)
{
    ERRORMSG ("PartDetect::InsideDetector2 - virtual dummy");
    doabort();
    return false;
}

//! record the particle info to detectionFile - all in SI units, also charge
inline void PartDetect::save(const TLinkedParticle* part, const gridreal r_new[3])
{
    (*fs) << Params::t << " " << part->popid << " " << part->w << " "
          << r_new[0] << " " <<  r_new[1] << " " <<  r_new[2] << " "
          << part->vx << " " << part->vy << " " << part->vz << "\n" << flush;
}

/** \brief Constructor for FieldDetectorSet
 * detectorFUNC input for a Field Detector is
 * point x1 y1 z1
 * point x2 y2 z2
 */
FieldDetectorSet::FieldDetectorSet(DetectorArgs args)
    : Detector(args)
{
    //check that necessary information is given
    if (args.detectorFUNC.given == false && args.coordinateFile.given == false) {
        ERRORMSG ("coordinates for a FieldDetectorSet must be given "
                  "in detectorFUNC or in coordinateFile");
        doabort();
    }
    //warn for unnessary information (args only used by PartDetectorSet)
    if (args.popIdStr.given == true) {
        WARNINGMSG2 ("FieldDetectorSet doesn't use popIdStr",args.popIdStr.value);
    }
    if (args.maxCounts.given == true) {
        WARNINGMSG ("FieldDetectorSet doesn't use maxCounts");
    }
    //warn for unnessary information (args only used by TestParticleSet)
    if (args.m.given == true) {
        WARNINGMSG ("FieldDetectorSet doesn't use mass");
    }
    if (args.q.given == true) {
        WARNINGMSG ("FieldDetectorSet doesn't use charge");
    }
    if (args.testParticleFile.given == true) {
        WARNINGMSG2 ("FieldDetectorSet doesn't use testParticleFile",
                     args.testParticleFile.value);
    }
    detectorType = "field";
    pointCoordinates.clear();
    // Read detectorFUNC args and create field detects
    // all checks are done before creating the detect objects
    if (args.detectorFUNC.given == true) {
        for (unsigned int i=0; i<args.detectorFUNC.name.size(); i++) {
            if (args.detectorFUNC.name[i].compare("point") != 0) {
                ERRORMSG2 ("not a 'point' type field detector - not created",
                           args.detectorFUNC.name[i]);
                //continuing nonetheless by ignoring this index i
                continue;
            }
            if (args.detectorFUNC.funcArgs[i].size()!=3) {
                ERRORMSG ("3 coordinates needed for a point type field detector - not created");
                //continuing nonetheless by ignoring this index i
                continue;
            }
            const gridreal coords[3] = {static_cast<gridreal>(args.detectorFUNC.funcArgs[i][0]),
                                        static_cast<gridreal>(args.detectorFUNC.funcArgs[i][1]),
                                        static_cast<gridreal>(args.detectorFUNC.funcArgs[i][2])
                                       };
            if (Params::insideBox(coords) == false) {
                char errmsg[80];
                sprintf(errmsg, "field detector %u outside the simulation box", i+1);
                // cannot continue! cellintpol will cause abort.
                ERRORMSG (errmsg);
                doabort();
            }
            stringstream nameStr;
            nameStr.clear();
            nameStr <<detectionFile1st << "_" << i+1 << "_"
                    << detectorFunctionNames[i];
            if (detectionFile2nd.compare("")!=0) {
                nameStr << "." << detectionFile2nd;
            }
            const string fileName(nameStr.str());
            ofstream *fs = new ofstream(fileName.c_str(),fstream::out);
            if (!fs->good()) {
                ERRORMSG2 ("unable to open detectionFile",fileName);
                continue; //continue nonetheless!
            }
            detectionFiles.push_back(fileName);
            Tgr3v *point = new Tgr3v(coords);
            pointCoordinates.push_back(*point);
            FieldDetect* fDetectTemp = NULL;
            fDetectTemp = new FieldDetect(fs, *point);
            fieldDetects.push_back(fDetectTemp);
        }
    }
    fieldDetectsFromFunc = fieldDetects.size();
    // Read coordinateFile
    if (args.coordinateFile.given == true) {
        readCoordinateFile(args.coordinateFile.value,pointCoordinates);
        vector<Tgr3v>::iterator i_vect = pointCoordinates.end();
        //create field point detector for all sets of coordinates
        int i_removed=0;
        for (int i=fieldDetectsFromFunc; i<(int) pointCoordinates.size(); i++,i_vect++) {
            stringstream nameStr;
            nameStr.clear();
            nameStr <<detectionFile1st << "_" << i+1 << "_"
                    << detectorFunctionNames[i] << "_fromfile";
            if (detectionFile2nd.compare("")!=0) {
                nameStr << "." << detectionFile2nd;
            }
            const string fileName(nameStr.str());
            ofstream *fs = new ofstream(fileName.c_str(),fstream::out);
            if (!fs->good()) {
                ERRORMSG2 ("unable to open detectionFile",fileName);
                pointCoordinates.erase(i_vect);
                i--;
                i_vect--;
                i_removed++;
                continue; //continue nonetheless!
            }
            detectionFiles.push_back(fileName);
            FieldDetect* fDetectTemp = NULL;
            fDetectTemp = new FieldDetect(fs, pointCoordinates[i]);
            fieldDetects.push_back(fDetectTemp);
        }
        if (pointCoordinates.size()!=fieldDetects.size()) {
            ERRORMSG ("internal error in FieldDetectorSet");
            doabort();
        }
    }
    fieldDetectsFromFile = fieldDetects.size() - fieldDetectsFromFunc;
}

//! Method for reading coordinateFile for point detectors in FieldDetectorSet
void FieldDetectorSet::readCoordinateFile(const string coordinateFile, vector <Tgr3v>& coordinateVector)
{
    FILE *coordFile = fopen(coordinateFile.c_str(), "r" );
    string detType = "point";
    if (coordFile == NULL) {
        ERRORMSG2 ("unable to open coordinateFile for a FieldDetectorSet",coordinateFile);
        return; //proceeding nonetheless!
    }
    gridreal coords[3];
    //coordinateFile has lines with values x y z in meters
    int q = fscanf(coordFile, "%f %f %f \n", &coords[0], &coords[1], &coords[2]);
    while (q==3 && Params::insideBox(coords) == true) { //if all assigned and within simulation box
        Tgr3v *point = new Tgr3v(coords);
        coordinateVector.push_back(*point);
        detectorFunctionNames.push_back(detType);
        q = fscanf(coordFile, "%f %f %f \n", &coords[0], &coords[1], &coords[2]);
    }
    fclose(coordFile);
}

// PARTICLE DETECTORS

//! Particle detector: Into a sphere
PartDetect* newPartDetect_SphereInto(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_SphereInto(fs, PartDetectArgs);
}
string SphereInto = "sphereInto";

//! Particle detector: Inside a sphere
PartDetect* newPartDetect_SphereInside(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_SphereInside(fs, PartDetectArgs);
}
string SphereInside = "sphereInside";

//! Particle detector: Out from a sphere
PartDetect* newPartDetect_SphereOut(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_SphereOut(fs, PartDetectArgs);
}
string SphereOut = "sphereOut";

//! Particle detector: X plane
PartDetect* newPartDetect_XPlane(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_XPlane(fs, PartDetectArgs);
}
string XPlane = "xplane";

//! Particle detector: X plane reverse
PartDetect* newPartDetect_XPlaneReverse(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_XPlaneReverse(fs, PartDetectArgs);
}
string XPlaneReverse = "xplaneReverse";

//! Particle detector: Line
PartDetect* newPartDetect_Line(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_Line(fs, PartDetectArgs);
}
string Line = "line";

//! Particle detector: Inside a line
PartDetect* newPartDetect_LineInside(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_LineInside(fs, PartDetectArgs);
}
string LineInside = "lineInside";

//! Particle detector: Ellipse
PartDetect* newPartDetect_Ellipse(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_Ellipse(fs, PartDetectArgs);
}
string Ellipse = "ellipse";

//! Particle detector: Inside an ellipse
PartDetect* newPartDetect_EllipseInside(ofstream *fs, vector <real> PartDetectArgs)
{
    return new PartDetect_EllipseInside(fs, PartDetectArgs);
}
string EllipseInside = "ellipseInside";

const unsigned int zero = 0;
/** \brief Constructor for ParticleDetectorSet
 *
 * All PartileDetectors have populationIdString as parameter - if "-" collect all!
 *
 * detectorFUNC input for a Particle Detector is
 * type R arg1 arg2 arg3 ...
 *
 * for type "xplane" the specific parameter is xplane[m]
 * for type "sphere" the sphere specific parameters are koordinates x[m] y[m] z[m]
 * for type "line" there are any given number of points (at least 2) with three coordinates each:
 *   x1[m] y1[m] z1[m] ... xN[m] yN[m] zN[m]
 * for type "ellipse" the specific parameters are
 *   a[m] b[m] e_a_x e_a_y e_a_z e_b_x e_b_y e_b_z
 *    here a is the semimajor axis and b is semiminor axis
 *    and e_a[3] is a vector in the direction of the apocenter and
 *    e_a[3] and e_b[3] determine the plane of the ellipse.
 * Note: origin of the grid is always a focal point
 */
ParticleDetectorSet::ParticleDetectorSet(DetectorArgs args)
    : Detector(args), numberOfPartDetectTypes(9) // update this number when adding new types!
{
    static bool setupPartFunc = false;
    static vector <PartDetect* (*) (ofstream*,vector<real>)> newPartDetectFuncs;
    static vector <string> partDetectTypes;
    if (setupPartFunc == false) {
        newPartDetectFuncs.clear();
        partDetectTypes.clear();
        PartDetect* (*npdf) (ofstream*, vector<real>);

        npdf = newPartDetect_SphereInto;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(SphereInto);

        npdf = newPartDetect_SphereInside;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(SphereInside);

        npdf = newPartDetect_SphereOut;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(SphereOut);

        npdf = newPartDetect_XPlane;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(XPlane);

        npdf = newPartDetect_XPlaneReverse;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(XPlaneReverse);

        npdf = newPartDetect_Line;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(Line);

        npdf = newPartDetect_LineInside;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(LineInside);

        npdf = newPartDetect_Ellipse;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(Ellipse);

        npdf = newPartDetect_EllipseInside;
        newPartDetectFuncs.push_back(npdf);
        partDetectTypes.push_back(EllipseInside);
        //add new detector type here
        if (ParticleDetectorSet::numberOfPartDetectTypes != newPartDetectFuncs.size()) {
            ERRORMSG ("inconsistent number of Particle Detect Types");
            doabort();
        }
        setupPartFunc = true;
    }

    //check that necessary information is given
    if (args.popIdStr.given == false) {
        ERRORMSG2 ("particle population identification string must be given "
                   "('-' for all populations)", popIdStr);
        doabort();
    } else if (args.popIdStr.value.compare("-") != 0) { //check that the popIdStr is valid
        bool foundIdStr = false;
        for (int i=0; i<Params::POPULATIONS; i++) {
            if (args.popIdStr.value.compare(Params::pops[i]->getIdStr()) == 0) {
                foundIdStr = true;
                break;
            }
        }
        if (foundIdStr == false) {
            ERRORMSG2 ("invalid particle population identification string for a particle detector set",
                       args.popIdStr.value);
            doabort();
        }
    }
    if (args.maxCounts.given == false) {
        ERRORMSG ("maximum counts to be must be given");
        doabort();
    }
    if (args.detectorFUNC.given == false) {
        ERRORMSG ("no particle detectorFUNC in 'detector particle' group");
        doabort();
    }
    //warn for unnessary information (args only used by TestParticleSet)
    if (args.m.given == true) {
        WARNINGMSG ("FieldDetectorSet doesn't use mass");
    }
    if (args.q.given == true) {
        WARNINGMSG ("FieldDetectorSet doesn't use charge");
    }
    if (args.testParticleFile.given == true) {
        WARNINGMSG2 ("FieldDetectorSet doesn't use testParticleFile",
                     args.testParticleFile.value);
    }
    detectorType = "particle";
    popIdStr = args.popIdStr.value;
    maxCounts = args.maxCounts.value;
    // Read detectorFUNC args and create particle detects
    // all parameter checks are done before creating detect objects
    for (unsigned int i=0; i<args.detectorFUNC.name.size(); i++) {
        if (args.detectorFUNC.name[i].compare("sphere") == 0) {
            args.detectorFUNC.name[i] = SphereInto;
        } //old sphere is sphereInto
        int typeId = -1;
        if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //SphereInto
            if (args.detectorFUNC.funcArgs[i].size()!=4) {
                ERRORMSG2 ("R, x, y, and z needed for a sphereInto type particle detector"
                           "\n - not created", args.detectorFUNC.name[i]);
                //continuing by ignoring this index i
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //SphereInside
            if (args.detectorFUNC.funcArgs[i].size()!=4) {
                ERRORMSG2 ("R, x, y, and z needed for a sphereInside type particle detector"
                           "\n - not created", args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //SphereOut
            if (args.detectorFUNC.funcArgs[i].size()!=4) {
                ERRORMSG2 ("R, x, y, and z needed for a sphereOut type particle detector"
                           "\n - not created", args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //Xplane (or disc)
            if (args.detectorFUNC.funcArgs[i].size()!=2) {
                ERRORMSG2 ("Radius R and coordinate x needed for a X-plane type "
                           "particle detector\n - not created",
                           args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //XplaneReverse
            if (args.detectorFUNC.funcArgs[i].size()!=2) {
                ERRORMSG2 ("Radius R and coordinate x needed for a X-planeReverse type "
                           "particle detector\n - not created",
                           args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //Line
            if ((args.detectorFUNC.funcArgs[i].size()-1)%3!=0 ||
                args.detectorFUNC.funcArgs[i].size()<7) {
                ERRORMSG2 ("R and at least 2 point coordinates needed for "
                           "a line type particle detector\n - not crated",
                           args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //LineInside
            if ((args.detectorFUNC.funcArgs[i].size()-1)%3!=0 ||
                args.detectorFUNC.funcArgs[i].size()<7) {
                ERRORMSG2 ("R and at least 2 point coordinates needed for "
                           "a line type particle detector\n - not crated",
                           args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //Ellipse
            if (args.detectorFUNC.funcArgs[i].size()!=9) {
                ERRORMSG2  ("R a b e_a[3] e_b[3] needed for an ellipse type particle "
                            "detector\n - not created", args.detectorFUNC.name[i]);
                continue;
            }
        } else if (args.detectorFUNC.name[i].compare(partDetectTypes[++typeId]) == 0) {
            //EllipseInside
            if (args.detectorFUNC.funcArgs[i].size()!=9) {
                ERRORMSG2  ("R a b e_a[3] e_b[3] needed for an ellipseInside type particle "
                            "detector\n - not created", args.detectorFUNC.name[i]);
                continue;
            }
            if (args.detectorFUNC.funcArgs[i][1]<args.detectorFUNC.funcArgs[i][2]) {
                ERRORMSG ("semiminor axis larger than semimajor axis\n - not created");
                continue;
            }
            //Add new detector type here by adding another } else if () { ...
        } else {
            ERRORMSG2 ("not a known type for a particle detector - not created",
                       args.detectorFUNC.name[i]);
            //continuing nonetheless by ignoring this index i
            continue;
        }
        if (typeId==-1) {
            ERRORMSG ("internal error in ParticleDetectorSet");
            doabort();
        }
        stringstream nameStr;
        nameStr.clear();
        nameStr <<detectionFile1st << "_" << i+1 << "_"
                << detectorFunctionNames[i];
        if (detectionFile2nd.compare("")!=0) {
            nameStr << "." << detectionFile2nd;
        }
        const string fileName(nameStr.str());
        ofstream *fs = new ofstream(fileName.c_str(),fstream::out);
        if (!fs->good()) {
            ERRORMSG2 ("unable to open detectionFile", fileName);
            continue; //continue nonetheless by ignoring this index i
        } else {
            //Write Header
            (*fs) << "% Particle Detector of type " << partDetectTypes[typeId];
            if (popIdStr.compare("-") == 0) {
                (*fs) << " recording all particle populations.\n"
                      "% Populations:   ";
                for (int popi=0; popi<Params::POPULATIONS; popi++) {
                    string temp = Params::pops[popi]->getIdStr();
                    (*fs) << temp.c_str() << " ";
                }
            } else {
                (*fs) << " recording one population: " << popIdStr;
            }
            (*fs) << "\n% t[s] popid w x[m] y[m] z[m] vx[m/s] vy[m/s] vz[m/s]\n" << flush;
            detectionFiles.push_back(fileName);
            PartDetect* pDetectTemp = NULL;
            pDetectTemp = (*newPartDetectFuncs[typeId])
                          (fs, args.detectorFUNC.funcArgs[i]);
            partDetects.push_back(pDetectTemp);
            pDetectTemp->firstline();
            currentCounts.push_back(zero);
        }
    }
}

//! TestParticleSet constructor
TestParticleSet::TestParticleSet(DetectorArgs args)
    : Detector(args)
{
    stillPropagating = true;
    if (args.popIdStr.given == true) {
        ERRORMSG2 ("testparticle 'detector' doesn't use popIdStr",
                   args.popIdStr.value);
    }
    if (args.maxCounts.given == true) {
        ERRORMSG ("testparticle 'detector' doesn't use maxCounts");
    }
    if (args.detectorFUNC.given == true) {
        ERRORMSG2 ("testparticle 'detector' doesn't use detectorFUNC",
                   args.detectorFUNC.name[0]);
    }
    if (args.testParticleFile.given == false) {
        ERRORMSG ("testParticleFile not given for testparticle 'detector'");
        doabort();
    }
    if (args.q.given == false) {
        ERRORMSG ("charge q not given for testparticle 'detector'");
        doabort();
    }
    if (args.m.given == false) {
        ERRORMSG ("mass m not given for testparticle 'detector'");
        doabort();
    }
    detectorType = "testparticle";
    testParticleFile = args.testParticleFile.value;
    charge = args.q.value;
    mass = args.m.value;
    //read testPartFile and create the TestParticle objects
    FILE *testPartFile = fopen(testParticleFile.c_str(), "r");
    string detType = "testparticle";
    if (testPartFile == NULL) {
        ERRORMSG2 ("unable to open testParticleFile for a test particle 'detector'",
                   testParticleFile);
        doabort();
    }
    gridreal r[3], v[3];
    //check all testparticles with Params::insideBox
    vector <int> outofBounds;
    outofBounds.clear();
    unsigned int tpcount = 0;
    //testParticleFile has lines with values x y z vx vy vz in meters and meters per second
    int q = fscanf(testPartFile, "%f %f %f %f %f %f\n", &r[0],&r[1],&r[2], &v[0],&v[1],&v[2]);
    while (q==6) { //if all assigned and within simulation box
        TestParticle *testpart = new TestParticle(r,v,charge,mass);
        testParts.push_back(testpart);
        if (Params::insideBox(r) == false) {
            outofBounds.push_back(tpcount);
            testParts[testParts.size()-1]->propagate = false;
        }
        tpcount++;
        q = fscanf(testPartFile, "%f %f %f %f %f %f\n", &r[0],&r[1],&r[2], &v[0],&v[1],&v[2]);
    }
    fclose(testPartFile);
    // listing all lines with incorrect (i.e. out of boundaries) testparticle positions
    if (outofBounds.empty() == false) {
        errorlog << "WARNING: Test particles inserted outside the simulation domain\n"
                 << "They are in file " << testParticleFile << " on rows: \n";
        for (unsigned int i=0; i<outofBounds.size(); i++) {
            errorlog << " " << outofBounds[i]+1;
        }
        errorlog << "\n" << flush;
    }
    //open testparticle save file (detectionFile)
    files = new ofstream(detectionFile.c_str(),fstream::out);
    if (!files->good()) {
        ERRORMSG2 ("unable to open detectionFile",detectionFile);
        doabort();
    }
    //Write Header
    // output file format: t x1 y1 z1 vx1 vy1 vz1 x2 y2 z2 vx2 vy2 vz2...
    (*files) << "% Testparticle detectionFile (" << testParts.size()
             << " testparticles)\n% particle mass is "
             << mass/Params::amu << " amu, electric charge is "
             << charge/Params::e << " e\n% t";
    for (unsigned int i=1; i<testParts.size()+1; i++) {
        (*files) << " x"<<i<<", y"<<i<<", z"<<i
                 << ", vx"<<i<<", vy"<<i<<", vz"<<i;
    }
    (*files) << "\n" << flush;
}

//! Store field values for all field detects
void Detector::runFieldDetects()
{
    if (detectorType.compare("field") != 0) {
        return;
    }
    if (Params::t>=detectionTime[0]) {
        if (Params::t<=detectionTime[1]) {
            for (unsigned int i=0; i<fieldDetects.size(); i++) {
                fieldDetects[i]->run();
            }
        } else { // when detection time over
            static bool closed=false;
            if (closed==false) {
                for (unsigned int i=0; i<fieldDetects.size(); i++) {
                    fieldDetects[i]->~FieldDetect();    //destructors
                }
            }
        }
    }
}

//! Run all particle detects
void Detector::runPartDetects(const TLinkedParticle* part, const gridreal r_new[3])
{
    if (detectorType.compare("particle") != 0) {
        return;
    }
    if (Params::t>=detectionTime[0]) {
        if (Params::t<=detectionTime[1]) { //detecting
            if (popIdStr.compare("-")==0 ||
                popIdStr.compare(Params::pops[part->popid]->getIdStr())==0) { //right pop
                for (unsigned int i=0; i<partDetects.size(); i++) {
                    if (currentCounts[i]<maxCounts) { //is the maxCounts already reached
                        currentCounts[i] += partDetects[i]->run(part,r_new);
                    }
                }
            }
        } else { // when detection time over
            static bool closed=false;
            if (closed==false) {
                for (unsigned int i=0; i<partDetects.size(); i++) {
                    partDetects[i]->~PartDetect(); //destructors
                }
            }
        }
    }
}

//! Run test particles
void Detector::runTestParticles(void)
{
    if (detectorType.compare("testparticle") != 0) {
        return;
    }
    if (Params::t>=detectionTime[0]) { //test particles are initialized but will start moving at this point
        if (Params::t<=detectionTime[1] && stillPropagating == true) { // propagating testparticles
            (*files) << Params::t; //saving time
            bool check = false; //check if any of the testparticles are still active
            for (unsigned int i=0; i<testParts.size(); i++) {
                if (testParts[i]->propagate == true) { //is the maxCounts already reached
                    testParts[i]->run();
                    check = true;
                }
                (*files) <<" "<<testParts[i]->x<<" "<<testParts[i]->y<<" "
                         <<testParts[i]->z<<" "<<testParts[i]->vx<<" "
                         <<testParts[i]->vy<<" "<<testParts[i]->vz;
            }
            (*files) << "\n" << flush;
            if (check == false) { //no testParts propagating
                stillPropagating = false;
            }
        }
    }
}

//! Sphere detect constructor: inputs are Radius and point of origin r[3] (Upper class)
PartDetect_Sphere::PartDetect_Sphere(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect(fs1, partDetectArgs)
{
    r[0] = partDetectArgs[1];
    r[1] = partDetectArgs[2];
    r[2] = partDetectArgs[3];
}

//! Dummy firstline function
void PartDetect_Sphere::firstline(void)
{
    ERRORMSG ("PartDetect_Sphere::firstline - virtual dummy");
    doabort();
}

inline bool PartDetect_Sphere::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (sqr(x-r[0]) + sqr(y-r[1]) + sqr(z-r[2]) < Radius2) {
        return true;    //in the sphere
    }
    return false;
}

inline bool PartDetect_Sphere::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    ERRORMSG ("PartDetect_Sphere::InsideDetector2 - virtual dummy");
    doabort();
    return false;
}

//! SphereInto detect constructor -- subclass of Sphere
PartDetect_SphereInto::PartDetect_SphereInto(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_Sphere(fs1, partDetectArgs)
{ }

void PartDetect_SphereInto::firstline(void)
{
    (*fs) << "% SphereInto Detector at x= " << r[0] << ", y= " << r[1] << ", z= " << r[2]
          << " [R_P] with R= " << R/Params::R_P <<" [R_P]\n" << flush;
}

inline bool PartDetect_SphereInto::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    return PartDetect_Sphere::InsideDetector(x,y,z);
}

inline bool PartDetect_SphereInto::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return this->InsideDetector(x,y,z);
}

//! SphereInto detect constructor -- subclass of Sphere
PartDetect_SphereInside::PartDetect_SphereInside(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_Sphere(fs1, partDetectArgs)
{ }

void PartDetect_SphereInside::firstline(void)
{
    (*fs) << "% SphereInside Detector at x= " << r[0] << ", y= " << r[1] << ", z= " << r[2]
          << " [R_P] with R= " << R/Params::R_P <<" [R_P]\n" << flush;
}

inline bool PartDetect_SphereInside::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    return PartDetect_Sphere::InsideDetector(x,y,z);
}

inline bool PartDetect_SphereInside::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return false; // always a hit when r_new inside
}

//! SphereInto detect constructor -- subclass of Sphere
PartDetect_SphereOut::PartDetect_SphereOut(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_Sphere(fs1, partDetectArgs)
{ }

void PartDetect_SphereOut::firstline(void)
{
    (*fs) << "% SphereOut Detector at x= " << r[0] << ", y= " << r[1] << ", z= " << r[2]
          << " [R_P] with R= " << R/Params::R_P <<" [R_P]\n" << flush;
}

inline bool PartDetect_SphereOut::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    return !(PartDetect_Sphere::InsideDetector(x,y,z));
}

inline bool PartDetect_SphereOut::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return this->InsideDetector(x,y,z);
}

//! XPlane detect constructor: inputs are Radius and x coordinate
PartDetect_XPlane::PartDetect_XPlane(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect(fs1, partDetectArgs)
{
    x_plane = partDetectArgs[1];
}

void PartDetect_XPlane::firstline(void)
{
    (*fs) << "% X-Plane Detector at x= " << x_plane/Params::R_P << " [R_P] with R= "
          << R/Params::R_P << " [R_P]\n" << flush;
}

inline bool PartDetect_XPlane::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (x < x_plane && sqr(y)+sqr(z) < Radius2) {
        return true;   //if moving to the sylinder
    }
    return false;
}

inline bool PartDetect_XPlane::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (x < x_plane) {
        return true;    //behind the plane
    }
    return false;
}

//! XPlaneReverse detect constructor -- subclass of XPlane
PartDetect_XPlaneReverse::PartDetect_XPlaneReverse(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_XPlane(fs1, partDetectArgs)
{ }

void PartDetect_XPlaneReverse::firstline(void)
{
    (*fs) << "% X-PlaneReverse Detector at x= " << x_plane/Params::R_P << " [R_P] with R= "
          << R/Params::R_P << " [R_P]\n" << flush;
}

inline bool PartDetect_XPlaneReverse::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (x > x_plane && sqr(y)+sqr(z) < Radius2) {
        return true;   //if moving to the sylinder
    }
    return false;
}

inline bool PartDetect_XPlaneReverse::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (x > x_plane) {
        return true;    //behind the plane
    }
    return false;
}

//! Line(s) detect constructor: inputs are Radius and point coordinates (at least two points)
PartDetect_Line::PartDetect_Line(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect(fs1, partDetectArgs)
{
    Npoints = (unsigned int)((partDetectArgs.size()-1)/3.0);
    Points.clear();
    directions.clear(); //directions between points
    for (unsigned int p_i=0; p_i<Npoints; p_i++) {
        Tgr3v *point = new Tgr3v(partDetectArgs[p_i*3+1], partDetectArgs[p_i*3+2],
                                 partDetectArgs[p_i*3+3]);
        Points.push_back(*point);
    }
    // min and max define the box were detection can happen
    // enables fast check in the insidefunction
    gridreal dcoords[3];
    for (unsigned int i=0; i<3; i++) {
        dcoords[i] = Points[1][i] - Points[0][i];
        mincoord[i] = min2(Points[0][i], Points[1][i]);
        maxcoord[i] = max2(Points[0][i], Points[1][i]);
    }
    normalize(dcoords);
    Tgr3v *dpoint = new Tgr3v(dcoords);
    directions.push_back(*dpoint);
    for (unsigned int p_i=2; p_i<Npoints; p_i++) {
        for (unsigned int i=0; i<3; i++) {
            dcoords[i] = Points[p_i][i] - Points[p_i-1][i];
            mincoord[i] = min2(mincoord[i], Points[p_i][i]);
            maxcoord[i] = max2(maxcoord[i], Points[p_i][i]);
        }
        normalize(dcoords);
        dpoint = new Tgr3v(dcoords);
        directions.push_back(*dpoint);
    }
    for (unsigned int i=0; i<3; i++) {
        mincoord[i] -= R;
        maxcoord[i] += R;
    }
}

void PartDetect_Line::firstline(void)
{
    (*fs) << "% Line Detector (" << Npoints << " lines) with R= "
          << R/Params::R_P << " [R_P]\n" << flush;
}

//! detect is a connected set of cylinders (R=Radius) around the line from point to point
inline bool PartDetect_Line::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    if (x < mincoord[0] || x > maxcoord[0] ||
        y < mincoord[1] || y > maxcoord[1] ||
        z < mincoord[2] || z > maxcoord[2]) {
        return false;
    }
    // check each line segment separately
    // three conditions: 1) on the 'right side' of starting point
    // 2) distance from the line determined by P_a and V
    // 3) on the 'left side' of the ending point of the segment
    bool check;
    gridreal Q[3];
    gridreal PR[3] = {x - Points[0][0], y - Points[0][1], z - Points[0][2]};
    for (unsigned int p_i=0; p_i<Npoints-1; p_i++) {
        check = false;
        // right side of the first line point
        if (dotProduct(directions[p_i].r, PR) > 0) {
            // distance to line
            crossProduct(directions[p_i].r, PR, Q);
            //norm of Q is the distance from the line
            if (vecsqr(Q) < Radius2) {
                check = true;
            }
        }
        PR[0] = x - Points[p_i+1].r[0];
        PR[1] = y - Points[p_i+1].r[1];
        PR[2] = z - Points[p_i+1].r[2];
        if (check) {
            if (dotProduct(directions[p_i].r, PR) < 0) { //left side of the next line point
                return true;
            }
        }
    }
    return false;
}

inline bool PartDetect_Line::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return this->InsideDetector(x,y,z);
}

//! LineInside detect constructor -- Line subclass
PartDetect_LineInside::PartDetect_LineInside(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_Line(fs1, partDetectArgs)
{ }

void PartDetect_LineInside::firstline(void)
{
    (*fs) << "% LineInside Detector (" << Npoints << " lines) with R= "
          << R/Params::R_P << " [R_P]\n" << flush;
}

//! detect is a connected set of cylinders (R=Radius) around the line from point to point
inline bool PartDetect_LineInside::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    return PartDetect_Line::InsideDetector(x,y,z);
}

inline bool PartDetect_LineInside::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return false; // always a hit when r_new inside
}

//! Ellipse detector constructor
PartDetect_Ellipse::PartDetect_Ellipse(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect(fs1, partDetectArgs)
{
    a = partDetectArgs[1];
    b = partDetectArgs[2]; //semimajor and semiminor axes
    e_a[0] = partDetectArgs[3];
    e_a[1] = partDetectArgs[4];
    e_a[2] = partDetectArgs[5];
    // direction of apocenter
    e_b[0] = partDetectArgs[6];
    e_b[1] = partDetectArgs[7];
    e_b[2] = partDetectArgs[8];
    // e_b is another direction determining the elliptical plane together with a_b
    const float checkPerpendicular = dotProduct(e_a,e_b);
    if (checkPerpendicular != 0) {
        char warnmsg[80];
        sprintf(warnmsg, "semimajor and semimajor directions might not be perpendicular"
                " e_a * e_b = %f", checkPerpendicular);
        WARNINGMSG2 (warnmsg, partDetectType);
    }
    crossProduct(e_a,e_b,e_c);  // e_c = e_a x e_b
    normalize(e_a);
    normalize(e_c);
    crossProduct(e_c,e_a,e_b); // now e_x:s are ortonormalized and right-handed;
    ecc = sqrt(1.0 - sqr(b/a));
    ae = a * ecc;
    if (a == b) {
        circle = true; // ecc = 0: a_ecc not used
    } else {
        circle = false;
        a_ecc = 0.5 * a * (1.0 - sqr(ecc)) / ecc;
    }
}

void PartDetect_Ellipse::firstline(void)
{
    (*fs) << "% Ellipse Detector (a= " << a/Params::R_P << " [R_P], eccentr.= "
          << ecc << ") with R= " << R/Params::R_P << " [R_P]\n" << flush;
}

/** /brief Point distance to Ellipse is not so simple - we approximate
 * Point distance to Ellipse is not so simple - we approximate
 *   still VERY time consuming - especially if eccentricity (or R) is large
 *     for large eccentricity
 *  The method uses projected lines from the focal points to estimate the distance
 *  - at worst the error seems to be between fractions 1/1000 and 1/300 of the Radius
 *    this is most obvious at the inner edge near peri- or apoapsis
 */
bool PartDetect_Ellipse::InsideDetector(const gridreal x, const gridreal y,
                                        const gridreal z)
{
    // take coordinate to elliptical system
    const gridreal r[3] = {x, y, z};
    const gridreal z_ell = dotProduct(r, e_c);
    // check distance to the elliptical plane
    if (z_ell > R || z_ell < -R) {
        return false;   // outside the plane
    }
    const gridreal x_ell = dotProduct(r, e_a); // focal point
    const gridreal y_ell = dotProduct(r, e_b);
    //check distance to central point // does this help any? three or four sqr-functions
    // and after this 3x sqrt! + 4 divisions + 11 sqr  - maybe it does.
    const gridreal y2 = sqr(y_ell);
    const gridreal r2 = sqr(x_ell - ae) + y2;
    // Projected distance from the center of ellipse
    if (r2 < sqr(b - R)) { // too close
        return false;
    }
    if (r2 > sqr(a + R)) { // too far
        return false;
    }
    // Final check - distance to the ellipse
    // 1)calculate two points of the ellipse using projected lines from both focal points
    //   towards the particle (x,y)
    // 2)take average (x_q) of the x components of these two intersection poins
    // 3)calculate corresponding y_q = sign(y)* b * sqrt(1 - sqr(x_q/a))
    // 4)calculate distance between (x_q,y_q,0) and particle (x,y,z)

    // ellipse orbit from a focal point r = a * (1 - sqr(ecc)) - ecc * x
    // x_q = (a(1-sqr(ecc)) -r) / ecc, wher r = a(1-sqr(ecc))/(1+ecc*cos(theta))
    //    where cos(theta) = (-/+) x_ell/r_ell
    //   so x_q = a(1-sqr(ecc))/ecc * (-/+)(1 - 1/(1 -/+ ecc*x_ell/r_ell))
    gridreal x_q;
    if (circle == true) { //ecc == 0
        x_q = a * x_ell/sqrt(r2);
    } else {
        const gridreal x_1 = 1.0/(1.0 - ecc * x_ell/sqrt(sqr(x_ell) + y2));
        const gridreal x_e2 = x_ell - 2*ae;
        const gridreal x_2 = 1.0/(1.0 + ecc * x_e2 /sqrt(sqr(x_e2) + y2));
        x_q = a_ecc *(x_1 - x_2); //(-ae +ae) - to central coordinates
        //a_ecc = 0.5 * a * (1-sqr(ecc))/ecc
    }
    const gridreal y_q = sign(y_ell) * b * sqrt(1.0 - sqr(x_q/a)); // in central coord.
    //const gridreal dist=sqrt(sqr(x_ell -ae -x_q) + sqr(y_ell-y_q) + sqr(z_ell));
    if (sqr(x_ell -ae -x_q) + sqr(y_ell-y_q) + sqr(z_ell) > Radius2) return false;
    return true;
}

inline bool PartDetect_Ellipse::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return this->InsideDetector(x,y,z);
}

//! EllipseInside detector constructor  -- subclass of Ellipse
PartDetect_EllipseInside::PartDetect_EllipseInside(ofstream *fs1, vector <real> partDetectArgs)
    : PartDetect_Ellipse(fs1, partDetectArgs)
{ }

void PartDetect_EllipseInside::firstline(void)
{
    (*fs) << "% EllipseInside Detector (a= " << a/Params::R_P << " [R_P], eccentr.= "
          << ecc << ") with R= " << R/Params::R_P << " [R_P]\n" << flush;
}

inline bool PartDetect_EllipseInside::InsideDetector(const gridreal x, const gridreal y,
        const gridreal z)
{
    return PartDetect_Ellipse::InsideDetector(x,y,z);
}

inline bool PartDetect_EllipseInside::InsideDetector2(const gridreal x, const gridreal y,
        const gridreal z)
{
    return false; // always a hit when r_new inside
}


//! Constructor
TestParticle::TestParticle(gridreal r[3], gridreal v[3], real charge, real mass)
{
    x = r[0];
    y = r[1];
    z = r[2];
    vx = v[0];
    vy = v[1];
    vz = v[2];
    propagate = true;
    q = charge;
    m = mass;
}

//! Boundary conditions for test particle
inline bool TestParticle::boundaries(void)
{
    const gridreal coords[3] = {x,y,z};
    return Params::insideBox(coords);
}

#ifndef USE_SPHERICAL_COORDINATE_SYSTEM

//! propagate position and check boundaries
inline bool TestParticle::run(void)
{
    PropagateV();
    x += vx*Params::dt;
    y += vy*Params::dt;
    z += vz*Params::dt;
    if (boundaries() == false) {
        propagate = false;
    }
    return propagate;
}

//! Accelerate particle (Lorentz force) taken from simulation.cpp (modified lines: ***)
void TestParticle::PropagateV(void)
{
    // Particle's centroid coordinates and velocity vectors
    fastreal r[3] = {x, y, z};  //***
    fastreal v[3] = {vx, vy, vz};  //***
    real B[3],Ue[3];
    // Self-consistent B1 field from cell faces + constant B0 field => B(r) = B1(r) + B0(r)
    g.faceintpol(r, Tgrid::FACEDATA_B, B);
    addConstantMagneticField(r, B);
    // Velocity field of the electron fluid from cells
    // Do not pass r => uses saved_cellptr and avoids findcell call
    g.cellintpol(Tgrid::CELLDATA_UE,Ue);
    if(Params::electronPressure==true) {
        real Efield[3],tx,ty,tz,sx,sy,sz,dvx,dvy,dvz,vmx,vmy,vmz,v0x,v0y,v0z,vpx,vpy,vpz,qmideltT2,t2,b2;
        // get the NGP electric field
        // Do not pass r => uses saved_cellptr and avoids findcell call
        g.cellintpol(Tgrid::CELLDATA_TEMP2,Efield);
        Efield[0] += B[1]*Ue[2] - B[2]*Ue[1];
        Efield[1] += B[2]*Ue[0] - B[0]*Ue[2];
        Efield[2] += B[0]*Ue[1] - B[1]*Ue[0];
        qmideltT2= 0.5*q*Params::dt/m; //***
        dvx=qmideltT2*Efield[0];
        dvy=qmideltT2*Efield[1];
        dvz=qmideltT2*Efield[2];
        tx=qmideltT2*B[0];
        ty=qmideltT2*B[1];
        tz=qmideltT2*B[2];
        t2=tx*tx+ty*ty+tz*tz;
        b2=2./(1.+t2);
        sx=b2*tx;
        sy=b2*ty;
        sz=b2*tz;
        vmx=v[0]+dvx;
        vmy=v[1]+dvy;
        vmz=v[2]+dvz;
        v0x=vmx+vmy*tz-vmz*ty;
        v0y=vmy+vmz*tx-vmx*tz;
        v0z=vmz+vmx*ty-vmy*tx;
        vpx=vmx+v0y*sz-v0z*sy;
        vpy=vmy+v0z*sx-v0x*sz;
        vpz=vmz+v0x*sy-v0y*sx;
        v[0]=vpx+dvx;
        v[1]=vpy+dvy;
        v[2]=vpz+dvz;
    } else {
        // Vector: dU = v_i - U_e
        real dU[3] = { v[0]-Ue[0], v[1]-Ue[1], v[2]-Ue[2] };
        // Constant: alpha/2 = q*dt/(2*m)
        const real half_alpha = 0.5*q*Params::dt/m; //***
        // Vector: W = q*dt*B/(2*m)
        real b[3] = {half_alpha*B[0], half_alpha*B[1], half_alpha*B[2]};
        // |W|^2
        const real b2 = sqr(b[0]) + sqr(b[1]) + sqr(b[2]);
        // Constant: beta = 2/(1+|b|^2)
        const real beta = 2.0/(1.0 + b2);
        // Cross products
        real dUxb[3],dUxbxb[3];
        crossProduct(dU,b,dUxb);
        crossProduct(dUxb,b,dUxbxb);
        // Add velocity components
        v[0] += beta*( dUxb[0] + dUxbxb[0] );
        v[1] += beta*( dUxb[1] + dUxbxb[1] );
        v[2] += beta*( dUxb[2] + dUxbxb[2] );
    }
    // Gravity correction
    if(Params::useGravitationalAcceleration == true) {
        real rLength = sqrt( sqr(x) + sqr(y) + sqr(z) ); //***
        real s = -Params::GMdt/cube(rLength);
        v[0] += s*x; //***
        v[1] += s*y; //***
        v[2] += s*z; //***
    }
    const real v2 = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    // Check particle maximum speed (CONSTRAINT)
    if (v2 > Params::vi_max2) {
        const real norm = Params::vi_max/sqrt(v2);
        v[0] *= norm;
        v[1] *= norm;
        v[2] *= norm;
        // Increase particle speed cutting rate counter
        //Params::pops[part.popid]->counter.cutRateV += 1.0; //***
    }
    vx = v[0]; //***
    vy = v[1]; //***
    vz = v[2]; //***
}

#else // spherical version

//! (SPHERICAL)
inline bool TestParticle::run(void)
{
    PropagateV();
    x += vx*Params::dt;
    y += vy*Params::dt;
    z += vz*Params::dt;
    // Transformation positions from hybrid to shperical coordinates
    gridreal r[3] = {x,y,z};
    sph_transf_H2S_R(r);
    // Cyclic condition for phi
    if (r[2] < 0.0)    r[2] = r[2] + 2.0*pi;
    if (r[2] > 2.0*pi) r[2] = r[2] - 2.0*pi;
    // Back to hybrid coordinates and velocities
    sph_transf_S2H_R(r);
    x = r[0];
    y = r[1];
    z = r[2];
    if (boundaries() == false) {
        propagate = false;
    }
    return propagate;
}

//! (SPHERICAL) Accelerate particle (Lorentz force) taken from simulation.cpp (modified lines: ***)
void TestParticle::PropagateV(void)
{
    // Particle's centroid coordinates and velocity vectors
    fastreal r[3] = {x, y, z};  //***
    fastreal v[3] = {vx, vy, vz};  //***
    real B[3],Ue[3];
    // Self-consistent B1 field from cell faces + constant B0 field => B(r) = B1(r) + B0(r)
    //!g.faceintpol(r, Tgrid::FACEDATA_B, B);
    g.sph_faceintpol(r, Tgrid::FACEDATA_B, B);
    addConstantMagneticField(r, B);
    // Velocity field of the electron fluid from cells
    // Do not pass r => uses saved_cellptr and avoids findcell call
    g.cellintpol(Tgrid::CELLDATA_UE,Ue);
    if(Params::electronPressure==true) {
        real Efield[3],tx,ty,tz,sx,sy,sz,dvx,dvy,dvz,vmx,vmy,vmz,v0x,v0y,v0z,vpx,vpy,vpz,qmideltT2,t2,b2;
        // get the NGP electric field
        // Do not pass r => uses saved_cellptr and avoids findcell call
        g.cellintpol(Tgrid::CELLDATA_TEMP2,Efield);
        Efield[0] += B[1]*Ue[2] - B[2]*Ue[1];
        Efield[1] += B[2]*Ue[0] - B[0]*Ue[2];
        Efield[2] += B[0]*Ue[1] - B[1]*Ue[0];
        qmideltT2= 0.5*q*Params::dt/m; //***
        dvx=qmideltT2*Efield[0];
        dvy=qmideltT2*Efield[1];
        dvz=qmideltT2*Efield[2];
        tx=qmideltT2*B[0];
        ty=qmideltT2*B[1];
        tz=qmideltT2*B[2];
        t2=tx*tx+ty*ty+tz*tz;
        b2=2./(1.+t2);
        sx=b2*tx;
        sy=b2*ty;
        sz=b2*tz;
        vmx=v[0]+dvx;
        vmy=v[1]+dvy;
        vmz=v[2]+dvz;
        v0x=vmx+vmy*tz-vmz*ty;
        v0y=vmy+vmz*tx-vmx*tz;
        v0z=vmz+vmx*ty-vmy*tx;
        vpx=vmx+v0y*sz-v0z*sy;
        vpy=vmy+v0z*sx-v0x*sz;
        vpz=vmz+v0x*sy-v0y*sx;
        v[0]=vpx+dvx;
        v[1]=vpy+dvy;
        v[2]=vpz+dvz;
    } else {
        // Vector: dU = v_i - U_e
        real dU[3] = { v[0]-Ue[0], v[1]-Ue[1], v[2]-Ue[2] };
        // Constant: alpha/2 = q*dt/(2*m)
        const real half_alpha = 0.5*q*Params::dt/m; //***
        // Vector: W = q*dt*B/(2*m)
        real b[3] = {half_alpha*B[0], half_alpha*B[1], half_alpha*B[2]};
        // |W|^2
        const real b2 = sqr(b[0]) + sqr(b[1]) + sqr(b[2]);
        // Constant: beta = 2/(1+|b|^2)
        const real beta = 2.0/(1.0 + b2);
        // Cross products
        real dUxb[3],dUxbxb[3];
        crossProduct(dU,b,dUxb);
        crossProduct(dUxb,b,dUxbxb);
        // Add velocity components
        v[0] += beta*( dUxb[0] + dUxbxb[0] );
        v[1] += beta*( dUxb[1] + dUxbxb[1] );
        v[2] += beta*( dUxb[2] + dUxbxb[2] );
    }
    // Gravity correction
    if(Params::useGravitationalAcceleration == true) {
        real rLength = abs(x);                 //***
        real s = -Params::GMdt/cube(rLength);
        v[0] += s*x;                           //***
    }
    // To calculate v2 in real space we need to transform velocities and positions from hybrid to shperical coordinates
    sph_transf_H2S_R(r);
    sph_transf_H2S_V(v);
    //const real v2 = sqr(v[0]) + sqr(v[1]) + sqr(v[2]);
    const real v2 = sqr(v[0]) + sqr(r[0]*v[1]) + sqr(r[0]*sin(r[1])*v[2]);
    // Check particle maximum speed (CONSTRAINT)
    if (v2 > Params::vi_max2) {
        const real norm = Params::vi_max/sqrt(v2);
        v[0] *= norm;
        v[1] *= norm;
        v[2] *= norm;
        // Increase particle speed cutting rate counter
        //Params::pops[part.popid]->counter.cutRateV += 1.0; //***
    }
    // Back to hybrid coordinates and velocities
    sph_transf_S2H_R(r);
    sph_transf_S2H_V(v);
    x = r[0];
    y = r[1];
    z = r[2];
    vx = v[0]; //***
    vy = v[1]; //***
    vz = v[2]; //***
}

#endif

