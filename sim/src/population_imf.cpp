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

#include <sstream>
#include "params.h"
#include "population.h"
#include "population_imf.h"
#include "simulation.h"
#include "random.h"
#include "definitions.h"

using namespace std;

extern Tgrid g;
extern Params simuConfig;

//! Constructor
PopulationIMF::PopulationIMF(PopulationArgs args) : Population(args)
{
    // Check the necessary information is given
    if(args.n.given == false) {
        ERRORMSG2("number density must be given",idStr);
        doabort();
    }
    if(args.V.given == false) {
        ERRORMSG2("bulk speed must be given",idStr);
        doabort();
    }
    // Initialize variable values
    n = 0;
    V = 0;
    rotVec = Tr3v(0, 0, 0);
    rotAngle = 0;
    backWallWeight = 1;
    negativeV = false;
    A0 = 0;
    usw = 0;
    uswCrossBdirMagn = 0;
    k = 0;
    v1 = 0;
    pitch_cutoff = 0;
    pitch_width = 0;
    pitch_distr_type = 0;
    xmin = 0;
    xmax = 0;
    ymin = 0;
    ymax = 0;
    zmin = 0;
    zmax = 0;
    box_x = 0;
    box_y = 0;
    box_z = 0;
    V_avg = 0;
    conserveE = false;
    t0 = 0;
    updateDistributions();
}

//! Destructor
PopulationIMF::~PopulationIMF() { }

//! Initialize imf population
void PopulationIMF::initialize()
{
    if(logParams == true && Params::t <= 0) {
        writeLog();
    }
}

//! Update distributions
bool PopulationIMF::updateDistributions()
{
    bool uswgiven, egiven, boxgiven;
    uswgiven = false;
    egiven = false, boxgiven = false;
    if(args.distFunc.name.size() <= 0) {
        ERRORMSG2("distribution functions not found",idStr);
        doabort();
    }
    for(unsigned int i = 0; i < args.distFunc.name.size(); ++i) {
        if(args.distFunc.name[i].compare("") == 0) {
            ERRORMSG2("trying to set empty distribution function name",idStr);
            doabort();
        } else if(args.distFunc.name[i].compare("usw") == 0) {
            if(args.distFunc.funcArgs[i].size() != 1) {
                ERRORMSG("usw line requires one argument (Usw..)");
                doabort();
            } else {
                usw = args.distFunc.funcArgs[i][0];
                uswgiven = true;
            }
        } else if(args.distFunc.name[i].compare("energyPowerLaw") == 0) {
            if(args.distFunc.funcArgs[i].size() != 2) {
                ERRORMSG("Energy_PowerLaw requires two arguments (upper v limit, spectral index)");
                doabort();
            } else {
                v1 = args.distFunc.funcArgs[i][0];
                k = args.distFunc.funcArgs[i][1];
                egiven = true;
                E_distr_type = 1;
            }
        } else if(args.distFunc.name[i].compare("Box") == 0) {
            if(args.distFunc.funcArgs[i].size() != 6) {
                ERRORMSG("Box requires six arguments (extents)");
                doabort();
            } else {
                xmin = args.distFunc.funcArgs[i][0];
                xmax = args.distFunc.funcArgs[i][1];
                ymin = args.distFunc.funcArgs[i][2];
                ymax = args.distFunc.funcArgs[i][3];
                zmin = args.distFunc.funcArgs[i][4];
                zmax = args.distFunc.funcArgs[i][5];
                box_x = xmax - xmin;
                box_y = ymax - ymin;
                box_z = zmax - zmin;
                boxgiven = true;
            }
        } else if(args.distFunc.name[i].compare("pitchAngleScatter") == 0) {
            if(args.distFunc.funcArgs[i].size() != 3) {
                ERRORMSG("pitchAngleScatter requires three arguments (type: 1-Gaussian 2-Constant 3-Isotropic; width/deg; cutoff/deg)");
                doabort();
            } else {
                pitch_distr_type = args.distFunc.funcArgs[i][0];
                if(pitch_distr_type > 3) {
                    ERRORMSG("no such pitch distribution");
                    doabort();
                }
                pitch_width = args.distFunc.funcArgs[i][1]/180*pi; //does nothing for const or isotropic
                pitch_cutoff = args.distFunc.funcArgs[i][2]/180*pi;
            }
        } else if(args.distFunc.name[i].compare("conserveE") == 0) {
            conserveE = true;
        } else if(args.distFunc.name[i].compare("injectStart") == 0) {
            t0 = args.distFunc.funcArgs[i][0];
        } else if(args.distFunc.name[i].compare("noBackWall") == 0) {
            backWallWeight = 0;
        }
    }
    if(!uswgiven) {
        usw = 0;
    }
    if(!egiven) {
        E_distr_type = 0;
    }
    return boxgiven;
}

//! Create solar wind population particles
void PopulationIMF::createParticles()
{
    if(Params::t > t0) {
        for (int i = 0; i < probround(macroParticlesPerDt); ++i) {
            newParticle();
        }
    }
}

//! Power law E spectrum, v*f(v)
real PopulationIMF::vf(real x, real pitch)
{
    return pow(x,2*k + 1)/sqr(x)/sin(pitch)/(2*pi) * pitchf(pitch);
}

//! Pitch angle
real PopulationIMF::pitchf(real a)
{
    if(pitch_distr_type == 1) { //Gaussian
        if(a > pitch_cutoff) {
            return 0;
        } else {
            return 1/sqrt(pitch_width*2*pi)*exp(-0.5*sqr(a)/(pitch_width));
        }
    } else if(pitch_distr_type == 2) {
        if(a > pitch_cutoff) { //constant
            return 0;
        } else {
            return 1;
        }
    } else if(pitch_distr_type ==3) { //isotropic
        if(a > pitch_cutoff) {
            return 0;
        } else {
            return sin(a);
        }
    } else {
        ERRORMSG("Dummy pitch angle distribution called, aborting");
        doabort();
    }
    return 0;
}


//! Midpoint integrator for EP distribution
real PopulationIMF::EnergyPitchIntegrator(Tr3v uswVec, real parkerAngle, real clockAngle) //implement: functions/tables to be passed
{
    double dv = 100e3, dpitch = pitch_cutoff/20, dphase = 2*pi/40;
    double vMin = V, vMax = v1, pitchMin = 0+0.5*dpitch, pitchMax = pitch_cutoff-0.5*dpitch, phaseMin = 0.5*dphase, phaseMax = 2*pi-0.5*dphase;
    Tr3v vVec, tmpVec;
    double result = 0, norm = 0;
    //norm over pitch angle
    if(pitch_distr_type == 1) { // gauss
        norm = -erf(-pitch_cutoff / sqrt(pitch_width * 2));
    } else if(pitch_distr_type == 2) { // const
        norm = pitch_cutoff;
    } else if(pitch_distr_type == 3) { // sin (isotropic)
        norm = cos(pitch_cutoff);
    } else {
        norm = 1;
    }
    for(double v = vMin; v < vMax; v += dv) {
        for(double pitch = pitchMin; pitch <= pitchMax; pitch += dpitch) {
            for(double phase = phaseMin; phase <= phaseMax; phase += dphase) {
                //init vVec
                vVec[0] = v;
                vVec[1] = 0;
                vVec[2] = 0;
                //rot to pitch
                tmpVec[0] = cos(pitch)*vVec(0);
                tmpVec[1] = sin(pitch)*vVec(0);
                tmpVec[2] = vVec(2);
                //rot around pitch axis (axis 0)
                vVec[0] = tmpVec(0);
                vVec[1] = cos(phase)*tmpVec(1);
                vVec[2] = sin(phase)*tmpVec(1);
                //rot to match B
                //tmpVec[0] = cos(parkerAngle)*vVec(0)-sin(parkerAngle)*vVec(1);
                //tmpVec[1] = sin(parkerAngle)*vVec(0)+cos(parkerAngle)*vVec(1);
                //tmpVec[2] = vVec(2);
                //with RotateVector
                tmpVec = RotateVector(tmpVec, rotAngle, rotVec);
                //translate to lab
                tmpVec = tmpVec + uswVec;
                //add to result
                result += tmpVec.magn()* vf(v,pitch) * sqr(v)*dv*dpitch*dphase*sin(pitch);
            }
        }
    }
    return 2*result/norm; //2x from dE->dv
}

//! Update imf population arguments
void PopulationIMF::updateArgs()
{
    bool boxgiven;
    double vdummy1 = 0, vdummy2 = 0;
    Population::updateBaseClassArgs(args);
    boxgiven = updateDistributions();
    //smaller injection box, if required for performance
    if(!boxgiven) {
        xmin = Params::box_xmin_tight;
        xmax = Params::box_xmax_tight;
        ymin = Params::box_ymin_tight;
        ymax = Params::box_ymax_tight;
        zmin = Params::box_zmin_tight;
        zmax = Params::box_zmax_tight;
        box_x = Params::box_X_tight;
        box_y = Params::box_Y_tight;
        box_z = Params::box_Z_tight;
    }
    // Set number density
    if(args.n.given == true) {
        if(args.n.value >= 0) {
            this->n = args.n.value;
        } else {
            ERRORMSG2("trying to set n < 0",idStr);
            doabort();
        }
    }
    Tr3v vVec;	//vVec for macrosperdt
    // Set bulk speed strictly along IMF
    if(args.V.given == true) {
        if(args.V.value >= 0) {
            this->V = args.V.value;
            negativeV = false;
        } else {
            this->V = -args.V.value;
            negativeV = true;
        }
        vVec = Tr3v(V, 0, 0);
        Tr3v B = Tr3v(Params::SW_Bx, Params::SW_By, Params::SW_Bz);
        rotVec = Cross(Tr3v(1, 0, 0), B/B.magn()); //rot in the plane perp to B and usw
        if(B(0) > 0) {
            B = -1*B;	// Flip B here to point antisunward for simplicity
        }
        rotAngle = acos(dot(Tr3v(1, 0, 0),B)/B.magn());
        vVec = RotateVector(vVec, rotAngle, rotVec);
        vdummy2 = v1;
        vdummy1 = V;
        //ERRORMSG2("RotAngle ", rotAngle/pi*180);	ERRORMSG2("RotVec(0) ", rotVec(0));
        //ERRORMSG2("RotVec(1) ", rotVec(1));	ERRORMSG2("RotVec(2) ", rotVec(2));
        if(vVec(0) > 0) {
            negativeV = true;
        }
    }
    // Set backwall weight
    if(args.backWallWeight.given == true) {
        if(args.backWallWeight.value >= 0) {
            this->backWallWeight = args.backWallWeight.value;
        } else {
            ERRORMSG2("trying to set backWallWeight < 0",idStr);
            doabort();
        }
    }
    //calculate box cross-section for principal direction
    real XYa, YZa, XZa, Atot;
    XYa = abs(vVec(2)/V) * box_x*box_y;
    YZa = abs(vVec(0)/V) * box_y*box_z;
    XZa = abs(vVec(1)/V) * box_x*box_z;
    Atot = XYa + YZa + XZa;
    this->A0 = Atot;
    // Find average velocity
    if( E_distr_type == 1) {
        double densityIntegral, vMoment;
        if(conserveE) {
            if(k != -1) {
                densityIntegral = (pow(vdummy2,(2*k+2))-pow(vdummy1,(2*k+2)))/(2*k+2);
            } else {
                densityIntegral = log(vdummy2)-log(vdummy1);
            }
            if(k != -1.5) {
                vMoment = (pow(vdummy2,(2*k+3))-pow(vdummy1,(2*k+3)))/(2*k+3);
            } else {
                vMoment = log(vdummy2)-log(vdummy1);
            }
        } else {
            if(k != -1) {
                densityIntegral = 1/(2*k+2)*(pow(vdummy2,(2*k+2))-pow(vdummy1,(2*k+2)));
            } else {
                densityIntegral = log(vdummy2)-log(vdummy1);
            }
            vMoment = EnergyPitchIntegrator(Tr3v(usw,0,0),rotAngle,0);
        }
        V_avg = vMoment/densityIntegral;
    } else {
        V_avg = V;
    }
    // Set macroParticlesPerDt
    int swPopsN = Params::popFactory.getNumberOfPopulations("imf") + Params::popFactory.getNumberOfPopulations("solarwind");
    if(args.macroParticlesPerDt.given == true) {
        this->macroParticlesPerDt = args.macroParticlesPerDt.value;
    } else {
        macroParticlesPerDt =  Params::nx*Params::ny*Params::nz*Params::macroParticlesPerCell/swPopsN * (box_x*box_y*box_z)/Params::box_V;
        macroParticlesPerDt *= Params::dt / (box_x*box_y*box_z) * V_avg * Atot;
    }
    // Set macroParticleStatisticalWeight so that principal direction gives correct results as is - corretion for deviations in individual particle weights
    if(macroParticlesPerDt <= 0) {
        macroParticleStatisticalWeight = 1;
    } else {
        real macroParticlesInFreeFlow = Params::nx*Params::ny*Params::nz*Params::macroParticlesPerCell * (box_x*box_y*box_z)/Params::box_V;
        real temp_sw_macros = macroParticlesInFreeFlow/swPopsN * Params::dt / (box_x*box_y*box_z) * V_avg * Atot;
        macroParticleStatisticalWeight = n*(box_x*box_y*box_z)*swPopsN/(macroParticlesInFreeFlow);
        // if macroParticlePerDt given in config
        if(macroParticlesPerDt != temp_sw_macros) {
            macroParticleStatisticalWeight *= temp_sw_macros/macroParticlesPerDt;
        }
    }
    // Write parameter log
    if(logParams == true && Params::t > 0) {
        writeLog();
    }
}

//! New imf population particle
void PopulationIMF::newParticle()
{
    real x,y,z;
    real XYa, YZa, XZa, Atot, rndPitch = 0, rndClock = 0;	// face areas and random numbers
    real BxExtra;
    Tr3v vVec_tmp, vVec;
    double a, b, V_tmp, v0, randomn;
    v0 = 0;
    real weight;
    weight = 1;
    // Launch ion along B_sw from front and side walls
    // create appropriate vel distr and rotate
    if (E_distr_type == 1) {	// Population distributed in energy as n = c*E^k
        if(k == -1.5) {
            b = pow(v1,2*k+3);
            randomn = uniformrnd();
            V_tmp = V*exp(randomn*log(v1/V));
        } else {
            a = pow(V,2*k+3);
            b = pow(v1,2*k+3);
            randomn = uniformrnd();
            V_tmp = pow((b-a)*randomn +a,1/(2*k+3));
        }
    } else { 	// Default: all particles with the same energy
        V_tmp = V;
    }
    if(usw != 0 && conserveE) {
        BxExtra = asin( uswCrossBdirMagn/V_tmp );
    } else {
        BxExtra = 0;
    }
    vVec_tmp = Tr3v(V_tmp, 0, 0); //Gaussian distr in popSW: vth*derivgaussrnd(V/vth), vth*gaussrnd(), vth*gaussrnd();
    // Rotate to initial distribution
    if(pitch_distr_type == 1) {
        do {
            rndPitch = abs(gaussrnd())*pitch_width;
        } while(abs(rndPitch) > pitch_cutoff);
    }
    if(pitch_distr_type == 2) {
        do {
            rndPitch = uniformrnd()*pitch_cutoff;
        } while(abs(rndPitch) > pitch_cutoff);
    }
    if(pitch_distr_type == 3) {
        do {
            rndPitch = (pi-acos(uniformrnd()*2-1));
        } while(abs(rndPitch) > pitch_cutoff);
    }
    vVec_tmp = Tr3v(cos(rndPitch)*vVec_tmp(0),sin(rndPitch)*vVec_tmp(0),vVec_tmp(2)); //initially points to +X
    rndClock = uniformrnd()*2*pi;
    vVec_tmp = Tr3v(vVec_tmp(0),cos(rndClock)*vVec_tmp(1),sin(rndClock)*vVec_tmp(1)); // no initial Z component
    if(conserveE) {
        // Rotation to corrected IMF direction
        //vVec_tmp = Tr3v(cos(rotAngle +BxExtra)*vVec_tmp(0)-sin(rotAngle +BxExtra)*vVec_tmp(1),sin(rotAngle +BxExtra)*vVec_tmp(0)+cos(rotAngle +BxExtra)*vVec_tmp(1),vVec_tmp(2));
        vVec = RotateVector(vVec_tmp, rotAngle+BxExtra, rotVec);
    } else {
        // Rotation to IMF direction and translation to co-moving coords
        //vVec_tmp = Tr3v(cos(rotAngle)*vVec_tmp(0)-sin(rotAngle)*vVec_tmp(1),sin(rotAngle)*vVec_tmp(0)+cos(rotAngle)*vVec_tmp(1),vVec_tmp(2));
        vVec = RotateVector(vVec_tmp,rotAngle,rotVec);
        vVec_tmp = vVec;
        vVec[0] += usw;
    }
    // Set up flow from walls chosen by V direction
    // face areas as seen from direction of V
    //vx_tmp = vVec_tmp(0); vy_tmp = vVec_tmp(1); vz_tmp = vVec_tmp(2);
    XYa = abs(vVec(2)/vVec.magn()) * box_x*box_y;
    YZa = abs(vVec(0)/vVec.magn()) * box_y*box_z;
    XZa = abs(vVec(1)/vVec.magn()) * box_x*box_z;
    Atot = XYa + YZa + XZa;
    //normalize
    XYa /= Atot;
    YZa /= Atot;
    XZa /= Atot;
    // choose face
    weight *= macroParticleStatisticalWeight;
    if(conserveE) {
        randomn = uniformrnd();
        if(randomn < XYa) {
            if(vVec(2) < 0) {
                z = zmax;
            } else {
                z = zmin;
            }
            y = ymin + uniformrnd()*box_y;
            x = xmin + uniformrnd()*box_x;
        } else if(randomn < XYa + YZa) {
            if(vVec(0) < 0) {
                x = xmax;
            } else {
                x =xmin;
            }
            y = ymin + uniformrnd()*box_y;
            z = zmin + uniformrnd()*box_z;
        } else {
            if(vVec(1) < 0) {
                y = ymax;
            } else {
                y = ymin;
            }
            z = zmin + uniformrnd()*box_z;
            x = xmin + uniformrnd()*box_x;
        }
        weight *= Atot/A0;
    } else {
        randomn = uniformrnd();
        if(randomn < XYa) {
            if(vVec(2) < 0) {
                z = zmax;
            } else {
                z = zmin;
            }
            y = ymin + uniformrnd()*box_y;
            x = xmin + uniformrnd()*box_x;
        } else if(randomn < XYa + YZa) {
            if(vVec(0) < 0) {
                x = xmax;
            } else {
                x = xmin;
                if(backWallWeight == 0) {
                    return;
                }
            }
            y = ymin + uniformrnd()*box_y;
            z = zmin + uniformrnd()*box_z;
        } else {
            if(vVec(1) < 0) {
                y = ymax;
            } else {
                y = ymin;
            }
            z = zmin + uniformrnd()*box_z;
            x = xmin + uniformrnd()*box_x;
        }
        v0 = vVec.magn()/vVec_tmp.magn();
        weight *= Atot/A0 * v0;
    }
    g.addparticle(x,y,z,vVec(0),vVec(1),vVec(2),weight,popid);
}

//! Nothing to write
void PopulationIMF::writeExtraHcFile() { }

//! Write population log
void PopulationIMF::writeLog()
{
    // Write particle log header lines
    if(logHeaderWritten == false) {
        // Create logger + logfile
        string logfile = "population_";
        logfile.append(idStr);
        logfile.append(".log");
        populationlog = Logger(logfile.c_str(),100000,false);
        // Initialize log (logfile is created at this point)
        populationlog.init();
        populationlog << "% particle population logfile (IMF)\n";
        populationlog << "% idStr = " << idStr << "\n";
        populationlog << "% popid = " <<  static_cast<int>(popid) << "\n";
        populationlog << "% m = " << m << "\n";
        populationlog << "% q = " << q << "\n";
        // File structure
        int nn = 1;
        string separator(1,'\t');
        populationlog << "% ";
        populationlog << int2string(nn++,2) << ". t [s]        " << separator;
        populationlog << int2string(nn++,2) << ". T [K]          " << separator;
        populationlog << int2string(nn++,2) << ". vth [m/s]      " << separator;
        populationlog << int2string(nn++,2) << ". mples/dt [#]   " << separator;
        populationlog << int2string(nn++,2) << ". statWeight [#] " << separator;
        populationlog << int2string(nn++,2) << ". split [-]      " << separator;
        populationlog << int2string(nn++,2) << ". join [-]       " << separator;
        populationlog << int2string(nn++,2) << ". n [#/m^3]      " << separator;
        populationlog << int2string(nn++,2) << ". V [m/s]        " << separator;
        populationlog << int2string(nn++,2) << ". backWall [-]   " << separator;
        populationlog << "\n";
        populationlog << scientific << showpos;
        populationlog.precision(10);
        logHeaderWritten = true;
    }
    string separator(2,' ');
    separator.append(1,'\t');
    populationlog << Params::t << separator;
    populationlog << T << separator;
    populationlog << vth << separator;
    populationlog << macroParticlesPerDt << separator;
    populationlog << macroParticleStatisticalWeight << separator;
    populationlog << static_cast<real>(split) << separator;
    populationlog << static_cast<real>(join) << separator;
    populationlog << n << separator;
    populationlog << V << separator;
    populationlog << backWallWeight << separator;
    populationlog << "\n";
}

//! Dump example population config (not implemented)
string PopulationIMF::configDump()
{
    return "";
}

//! Returns string representation of the population
string PopulationIMF::toString()
{
    stringstream ss;
    if(negativeV == false) {
        ss << "type = IMF\n";
    } else {
        ss << "type = IMF (from the back wall)\n";
    }
    ss << Population::toStringGeneral();
    ss << "n = " << n/1e6 << " #/cm^3\n";
    // ss << "V = " << vVec/1e3 << "km/s\n";
    /* if(negativeV == false)
       {
    ss << "V = " << V/1e3 << " km/s\n";
       }
     else
       {
    ss << "V = " << -V/1e3 << " km/s\n";
       }  */
    ss << "backWallWeight = " << backWallWeight << "\n";
    return ss.str();
}

