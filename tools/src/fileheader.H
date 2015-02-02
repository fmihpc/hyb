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

#ifndef FILEHEADER_H

#ifdef __GNUC__
#  pragma interface
#endif

#include <iostream>
using namespace std;

// Names in Theader are not case sensitive

class Theader {
public:
	enum TItemType {INT,REAL,STRING};
private:
	struct THeaderDataItem {
		THeaderDataItem *next;
		char *name;
		TItemType type;
		long i;
		double x;
		char *s;
	};
public:
	class print {
	private:
		TItemType t;
	public:
		print(TItemType t1) {t=t1;}
		friend ostream& operator<<(ostream& o, const print& p);
	};
private:
	bool dirty;
	THeaderDataItem *lst;
	THeaderDataItem *find(const char *name) const;
	bool check(const char *name, TItemType t, bool probeflag=false) const;
public:
	Theader() {lst=0; dirty=true;}
	Theader(const char *fn);
	bool good() const {return !dirty;}
	bool exists(const char *name) const {return check(name,INT,true);}
	TItemType gettype(const char *name) const;
	long getint(const char *name) const;
	char *getstr(const char *name) const;
	double getreal(const char *name) const;
	friend ostream& operator<<(ostream& o, const Theader& h);
	friend istream& operator>>(istream& o, Theader& h);
	~Theader();
};

#define FILEHEADER_H

#endif
