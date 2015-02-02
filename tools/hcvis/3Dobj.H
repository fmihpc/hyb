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

#ifndef THREEDOBJ_H

#ifdef __GNUC__
#  pragma interface
#endif

class T3DObject {
public:
	T3DObject *next;
	virtual const char *type() const = 0;
	virtual ~T3DObject() {}
};

class T3DSliceObject : public T3DObject {
private:
	int d;
	double xval;
public:
	bool DontShowIn3D;
	T3DSliceObject(int d1, double xval1) : DontShowIn3D(false) {d=d1; xval=xval1;}
	virtual const char *type() const {return "Slice";}
	int SlicedDim() const {return d;}
	double Slice_Xval() const {return xval;}
	virtual ~T3DSliceObject() {}
};

class T3DSphereObject : public T3DObject {
public:
	double x0,y0,z0;
	double r;
	double deltatheta;
	double deltatheta_grid;
public:
	T3DSphereObject(double r1, double x1, double y1, double z1, double deltatheta1, double deltatheta_grid1) {
		r=r1; x0=x1; y0=y1; z0=z1; deltatheta=deltatheta1; deltatheta_grid=deltatheta_grid1;
	}
	virtual const char *type() const {return "Sphere";}
	virtual ~T3DSphereObject() {}
};

enum TVectorField {VectorFieldRhoV, VectorFieldB, VectorFieldB0, VectorFieldB1, VectorFieldJ,
				   VectorFieldPoynting, VectorFieldEnergyFlux};
enum TTraceDirection {PositiveTraceDirection, NegativeTraceDirection, BothTraceDirections};
enum TFieldLineDistribution {UniformDistribution, InvSqrtDistribution, InvDistribution};
enum TLoopThresholdType {LocalSpacingLoopThreshold,BasegridSpacingLoopThreshold,AbsoluteLoopThreshold};

class T3DFieldLineBunchObject : public T3DObject {
public:
	double r1[3],r2[3];
	int n;
	TVectorField VectorField;
	TTraceDirection TraceDirection;
	TFieldLineDistribution Distribution;
	TLoopThresholdType LoopThresholdType;
	double LoopThreshold;
	T3DFieldLineBunchObject(int,const double[3],const double[3],
							TVectorField,TTraceDirection,TFieldLineDistribution,TLoopThresholdType,double);
	virtual const char *type() const {return "FieldLineBunch";}
	virtual ~T3DFieldLineBunchObject() {}
};

#define THREEDOBJ_H

#endif
