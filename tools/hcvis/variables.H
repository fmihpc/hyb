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

#ifndef VARIABLES_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <vector>
#include "metagrid.H"

class Tvariable {
private:
	int var;
	bool logflag;
	bool pseudobackgroundflag;
	double gamma,invmu0,mass;
	bool isSpectraFile;
	double *mintab, *maxtab;
	bool *minmaxexiststab;
public:
	Tvariable();
	int Nvars() const;
	bool select(const char *, double gamma, double invmu0, double mass);
	void select(int i);
	bool check(const char *);
	void logarithmic(bool flag) {logflag=flag;}
	void pseudobackground(bool flag) {pseudobackgroundflag=flag;}
	const char *selected() const;
	const char *description() const;
	double get(Tmetagrid& g, const Tdimvec& X) const;
	vector<double> getSpectra(Tmetagrid& g, const Tdimvec& X) const;
	void clearAllMinMax();
	void setMinMax(double m, double M) {mintab[var]=m; maxtab[var]=M; minmaxexiststab[var]=true;}
	void getMinMax(double& m, double& M) const;
	bool hasMinMax() const {return minmaxexiststab[var];}
	~Tvariable();
};

extern void ComputeCurl(Tmetagrid& g, const Tdimvec& X, smallnat a0, real& jx, real& jy, real& jz);

#define VARIABLES_H

#endif
