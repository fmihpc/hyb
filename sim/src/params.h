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

#ifndef PARAMS_H
#define PARAMS_H

#include <string>
#include <vector>
#include "definitions.h"
#include "population.h"
#include "detector.h"
#include "refinement.h"
#include "resistivity.h"
#include "forbidsplitjoin.h"
#include "backgroundcharge.h"
#include "magneticfield.h"
#include "splitjoin.h"
#include "diagnostics.h"
#include "particle.h"

struct PopulationArgs;
class Population;
class PopulationFactory;
struct DetectorArgs;
class Detector;
class DetectorFactory;
struct ProcessArgs;

//Some macros to ease the use of addVariable() method
#define ADD_INT(varName, comment) addVariable((void *)& varName, TYPE_INT, #varName, comment)
#define ADD_INT_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_INT, #varName, comment, tableSize)
#define ADD_BOOL(varName, comment) addVariable((void *)& varName, TYPE_BOOL, #varName, comment)
#define ADD_BOOL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_BOOL, #varName, comment, tableSize)
#define ADD_GRIDREAL(varName, comment) addVariable((void *)& varName, TYPE_GRIDREAL, #varName, comment)
#define ADD_GRIDREAL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_GRIDREAL, #varName, comment, tableSize)
#define ADD_DATAREAL(varName, comment) addVariable((void *)& varName, TYPE_DATAREAL, #varName, comment)
#define ADD_DATAREAL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_DATAREAL, #varName, comment, tableSize)
#define ADD_TPDF_IDL(varName, comment) addVariable((void *)& varName, TYPE_TPDF_IDL, #varName, comment)
#define ADD_TPDF_IDL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_TPDF_IDL, #varName, comment, tableSize)
#define ADD_SHORTREAL(varName, comment) addVariable((void *)& varName, TYPE_SHORTREAL, #varName, comment)
#define ADD_SHORTREAL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_SHORTREAL, #varName, comment, tableSize)
#define ADD_FASTREAL(varName, comment) addVariable((void *)& varName, TYPE_FASTREAL, #varName, comment )
#define ADD_FASTREAL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_FASTREAL, #varName, comment, tableSize)
#define ADD_REAL(varName, comment) addVariable((void *)& varName, TYPE_REAL, #varName, comment)
#define ADD_REAL_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_REAL, #varName, comment, tableSize)
#define ADD_STRING(varName, comment) addVariable((void *)& varName, TYPE_STRING, #varName, comment)
#define ADD_STRING_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_STRING, #varName, comment, tableSize)
#define ADD_CHAR(varName, comment) addVariable((void *)& varName, TYPE_CHAR, #varName, comment)
#define ADD_CHAR_TBL(varName, comment, tableSize) addVariable((void *)& varName, TYPE_CHAR, #varName, comment, tableSize)
#define ADD_FUNCTION(varName, comment) addVariable((void *)& varName, TYPE_FUNCTION, #varName, comment)
#define ADD_POPULATION(varName, comment) addVariable((void *)& varName, TYPE_POPULATION, #varName, comment)
#define ADD_DETECTOR(varName, comment) addVariable((void *)& varName, TYPE_DETECTOR, #varName, comment)
#define ADD_PROCESS(varName, comment) addVariable((void *)& varName, TYPE_PROCESS, #varName, comment)

//! List of allowed types
enum DATATYPE {TYPE_INT, TYPE_BOOL, TYPE_GRIDREAL, TYPE_DATAREAL, TYPE_TPDF_ID, TYPE_SHORTREAL, TYPE_FASTREAL, TYPE_REAL, TYPE_STRING, TYPE_CHAR, TYPE_FUNCTION, TYPE_POPULATION, TYPE_DETECTOR, TYPE_PROCESS};

/** \brief Pointers to variables with types.
 *
 * All assignements should be done through pointer of proper type.
 * This ensures valid behaviour of the variables.
 */
union varPtr {
    int	       *ptr_int;
    bool	       *ptr_bool;
    gridreal    *ptr_gridreal;
    datareal    *ptr_datareal;
    TPDF_ID     *ptr_TPDF_ID;
    shortreal   *ptr_shortreal;
    fastreal    *ptr_fastreal;
    real        *ptr_real;
    std::string *ptr_string;
    char        *ptr_char;
    void        *ptr_void;
};

//! Structure to contain pointer to the variable, and some other necessary information.
struct dynamicVar {
    union varPtr varPtr;
    DATATYPE type;
    int tableSize;
    char *name;
    char *comment;
    bool constant;
    bool initconstant;
    bool dumpping;
    bool updatedFromFile;
    void (*action)(void *);
    char **expressions; 		//ptr to tbl of expressions
    struct dynamicVar *next;
};

/** \brief Input parameter handling and simulation parameters
 *
 * Class Params includes simulation parameters, constants
 * and initial value constants and their cross dependencies (e.g. T
 * <-> v_th). The class also handles simulation input parameters and
 * their dynamics (dynamic configuration).
 */
class Params
{
public:
    Params();
    ~Params();
    void init(const bool dumpVariables = false);
    static std::string getSimuTimeStr();
    static bool onlyOneObject;
    static const std::string codeVersion;
    static bool stoppingPhase;
    const static char* configFileName;
    static bool wsFileGiven;
    static std::string wsFileName;
    // physical constants
    static const real c;
    static const real e;
    static const real k_B;
    static const real mu_0;
    static const real eps_0;
    static const real G;
    static const real AU;
    static const real gamma;
    static const real amu;
    static const real m_e;
    static const real m_p;
    static const real m_H;
    static const real m_H2;
    static const real m_2H;
    static const real m_3H;
    static const real m_3He;
    static const real m_He;
    static const real m_N;
    static const real m_N2;
    static const real m_O;
    static const real m_O2;
    static const real m_CH3;
    static const real m_CH4;
    static const real m_Na;
    static const real m_C2H5;
    static const real m_C2H6;
    static const real m_C3H8;
    static const real M_Me;
    static const real R_Me;
    static const real M_V;
    static const real R_V;
    static const real M_Mo;
    static const real R_Mo;
    static const real M_Ma;
    static const real R_Ma;
    static const real M_T;
    static const real R_T;
    static const real M_Pl;
    static const real R_Pl;
    // config file template variables
    static real tempRealA;
    static real tempRealB;
    static real tempRealC;
    static int tempIntA;
    static int tempIntB;
    static int tempIntC;
    // profiles
    static GridRefinementProfile gridRefinementFunction;
    static ResistivityProfile resistivityFunction;
    static ForbidSplitAndJoinProfile forbidSplitAndJoinFunction;
    static BackgroundChargeDensityProfile bgChargeDensityFunction;
    static std::vector<MagneticFieldProfile> initialMagneticFieldProfile;
    static std::vector<MagneticFieldProfile> constantMagneticFieldProfile;
    // particle populations
    static const int MAX_POPULATIONS;
    static int POPULATIONS;
    static std::vector<Population*> pops;
    static PopulationFactory popFactory;
    // TEMPLATE POPULATION VARIABLES
    static std::string population;
    static std::string idStr;
    static std::string hcFilePrefix;
    static std::string boundaryFUNC;
    static real m;
    static real q;
    static real T;
    static real vth;
    static real macroParticlesPerDt;
    static bool propagateV;
    static bool accumulate;
    static bool split;
    static bool join;
    static bool logParams;
    static real n;
    static real V;
    static real backWallWeight;
    static real R;
    static real totalRate;
    static std::string distFunc;
    // population arguments
    void clearPopulationVars();
    PopulationArgs getPopulationArgsStruct();
    // diagnostics and detectors
    static Diagnostics diag;
    static const int MAX_DETECTORS;
    static int DETECTORS;
    static std::vector<Detector*> detectors;
    static DetectorFactory detectorFactory;
    // TEMPLATE DETECTOR VARIABLES
    static std::string detector;
    static std::string popIdStr;
    static std::string detectionFile;
    static real detectionTime[2];
    static real maxCounts;
    static std::string coordinateFile;
    static std::string testParticleFile;
    static std::string detectorFUNC;
    // detector arguments
    void clearDetectorVars();
    DetectorArgs getDetectorArgsStruct();
    // TEMPLATE PROCESS VARIABLES
    static std::string process;
    static std::string procIdStr;
    static std::string incidentIonIdStr;
    static std::string exoNeutralCoronaIdStr;
    static std::string ENAIdStr;
    static std::string slowIonIdStr;
    static fastreal crossSection;
    static fastreal weightFactor;
    static fastreal slowIonVth;
    static real N_limitHeavyReactions;
    static real probLimitHeavyReactions;
    // /process arguments
    void clearProcessVars();
    ProcessArgs getProcessArgsStruct();
    // simulation parameters
    static int objectIdHWA;
    static real R_P;
    static real M_P;
    static real R_zeroFields;
    static real R_zeroFields2;
    static real R_zeroPolarizationField;
    static real dtField;
    static bool propagateField;
    static int nx;
    static int ny;
    static int nz;
    static real dx;
    static fastreal box_xmin;
    static fastreal box_xmax;
    static fastreal box_ymin;
    static fastreal box_ymax;
    static fastreal box_zmin;
    static fastreal box_zmax;
    static fastreal box_eps;
    static bool insideBox(const gridreal coords[3]);
    static bool insideBoxTight(const TLinkedParticle*);
    static bool insideBoxTightFrontWall(const TLinkedParticle*);
    static bool insideBoxTightBackWall(const TLinkedParticle*);
    static bool insideBoxTightSideWall(const TLinkedParticle*);
    static bool insideBoxTight(const gridreal coords[3]);
    static fastreal box_X;
    static fastreal box_Y;
    static fastreal box_Z;
    static fastreal box_V;
    static fastreal box_xmin_tight;
    static fastreal box_xmax_tight;
    static fastreal box_ymin_tight;
    static fastreal box_ymax_tight;
    static fastreal box_zmin_tight;
    static fastreal box_zmax_tight;
    static fastreal box_X_tight;
    static fastreal box_Y_tight;
    static fastreal box_Z_tight;
    static real SW_Bx;
    static real SW_By;
    static real SW_Bz;
    static real B_limit;
    static real Ecut;
    static bool Bboundaries[];
    static real boundary_Bx;
    static real boundary_By;
    static real boundary_Bz;
    static std::string initialMagneticFieldFUNC;
    static std::string constantMagneticFieldFUNC;
    static real dt;
#ifdef USE_PARTICLE_SUBCYCLING
    static int subcycleMaxLevel;
    static int subcycleType;
    static real dt_psub[256];
    static real accum_psubfactor[256];
#endif
    static real t;
    static real t_max;
    static int cnt_dt;
    static real saveInterval;
    static int saveHC;
    static int saveVTK;
    static bool averaging;
    static bool plasma_hcfile;
    static bool dbug_hcfile;
    static bool bg_in_avehcfile;
    static bool saveExtraHcFiles;
    static real wsDumpInterval[2];
    static real inputInterval;
    static real logInterval;
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    static bool saveParticlesAlongOrbit;
    static std::string saveParticlesAlongOrbitFile;
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    static real spectraEmin_eV;
    static real spectraEmax_eV;
    static int spectraNbins;
    static bool spectraLogBins;
    static bool spectraEminAll;
    static bool spectraEmaxAll;
    static int spectraMethod;
    static std::string spectraUnit;
    static std::vector<real> spectraEnergyBins_eV;
    static std::vector<real> spectra_dE_eV;
    static std::vector< std::vector<real> > spectraV2BinsPerPop;
#endif
    static bool fieldPredCor;
    static bool electronPressure;
    static real Te;
    static bool useGravitationalAcceleration;
    static real GMdt;
    static int macroParticlesPerCell;
    static bool useMacroParticleSplitting;
    static bool useMacroParticleJoining;
    static Split splittingFunction;
    static Join joiningFunction;
    static std::string splitFUNC;
    static std::string joinFUNC;
    static std::string gridRefinementFUNC;
    static int maxGridRefinementLevel;
    static int currentGridRefinementLevel;
    static std::string forbidSplitAndJoinFUNC;
    static std::string bgChargeDensityFUNC;
    static real splitJoinDeviation[2];
    static real splitjoin_a;
    static real vi_max;
    static real vi_max2;
    static real Ue_max;
    static real Ue_max2;
    static real rho_q_min;
    static real maxVw;
    static int densitySmoothingNumber;
    static int electricFieldSmoothingNumber;
    static std::string resistivityFUNC;
    // jstag
    static bool useJstag;
    static bool useNodeUe;
    // Titan specific
    static real SaturnLocalTime;
    static real SubSolarLatitude;
    // functions
    void addVariable(void *ptr, DATATYPE type, const char *name, const char *comment="", int tableSize=1);
    void addAction(const char varName[], void (*action)(void *value)=(void(*)(void*))0);
    void makeConstant(const char varName[]);
    void makeInitConstant(const char varName[]);
    void setVarDumppingOff(const char varName[]);
    void readAndUpdateVariables(const char *dumpFile);
    void readAndInitVariables(const char *dumpFile);
    void dumpVars(const char *dumpFile);
    void outParams();
    bool getFunctionNamesAndArgs(const char varName[], std::vector<std::string> &funcNames, std::vector< std::vector<real> > &args);
    std::string getFunctionName(const char varName[]);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    // simulation variables
    static real sph_dy;
    static real sph_dz;
    static real sph_dtheta;
    static real sph_dphi;
    static fastreal sph_alpha_theta;
    static fastreal sph_alpha_phi;
    static fastreal sph_theta_min;
    static fastreal sph_theta_max;
    static fastreal sph_phi_min;
    static fastreal sph_phi_max;
    static fastreal sph_box_theta;
    static fastreal sph_box_phi;
    static fastreal sph_box_V;
    static int sph_BC_type;
    static int sph_BC_use_ghost_cell;
    static int sph_propagation_type;
    static int sph_propagation_dir;
    static int sph_coordinate_grid_visual;
#endif
private:
    std::vector<real> parseFunctionNameAndArguments(std::string str, std::string &funcName);
    struct dynamicVar *varList;
    static bool initPhase;
    static bool readingPopulation;
    static bool readingDetector;
    static bool readingProcess;
    struct dynamicVar * getLast();
    struct dynamicVar * lookupVar(const char *name);
    static void setInitialValues();
    static void updateDependantParameters();
    static std::vector<std::string> idStrTbl;
    static bool checkIdStrAlreadyFound(std::string str);
    static std::vector<std::string> procIdStrTbl;
    static bool checkProcIdStrAlreadyFound(std::string str);
    void initVariables();
    int readAndUpdateVariable(std::ifstream &varDump, struct dynamicVar *var, int index);
    std::string readFunctionTypeVar(std::ifstream &fileStream, struct dynamicVar *var);
    void resetUpdatedFromFileFlags();
    int getSize(DATATYPE type);
    void setVar(struct dynamicVar *var,  void * value, int index=0);
    real getRealValue(const char*varName, int index=0);
    std::string getStringValue(const char*varName, int index=0);
    double readAndEvaluateRPNExpression(const char *expression);
    std::string getType(DATATYPE);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    static void sph_setInitialValues();
#endif
};

#endif

