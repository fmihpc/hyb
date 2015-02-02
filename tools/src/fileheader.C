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

#ifdef __GNUC__
#  pragma implementation "fileheader.H"
#endif

#include "fileheader.H"
#include <cstdlib>
#include <limits.h>		/* to get LONG_MIN */
#include <cstring>
#ifdef AIX
#  include <strings.h>
#endif
#include <cctype>
#include <fstream>
using namespace std;

istream& operator>>(istream& o, Theader& h)
{
	int i,j;
	h.dirty = false;
	const int maxlinelen = 1025;
	char buff[maxlinelen],name[maxlinelen],val[maxlinelen];
	while (1 /*o.getline(buff,maxlinelen-1)*/ ) {
		o.getline(buff,maxlinelen-1);
		if (!o.good()) break;
		if (buff[0] == '#') continue;
		for (i=0; buff[i]; i++) buff[i] = tolower(buff[i]);
		if (!strcmp(buff,"eoh")) break;
		for (i=0; buff[i] && buff[i] != '='; i++) name[i] = buff[i];
		name[i] = '\0';
		int len = strlen(name);
		len--;
		while (isspace(name[len])) {name[len] = '\0'; len--;}
		if (buff[i] != '=') {
			cerr << "*** Theader::Theader: syntax error: no '=' char in input line:\n";
			cerr << "    " << buff << "\n";
			h.dirty = true;
			break; /*continue;*/
		}
		i++;
		while (isspace(buff[i])) i++;
		int will_be_string = 0;
		if (buff[i] == '"') {
			i++;
			for (j=0; buff[i] && buff[i] != '"'; i++,j++) {
				char ch = buff[i];
				if (ch == '\\') {
					ch = buff[++i];
					if (ch == 'n')
						ch = '\n';
					else if (ch == 't')
						ch = '\t';
				}
				val[j] = ch;
			}
			will_be_string = 1;
		} else
			for (j=0; buff[i]; i++,j++) val[j] = buff[i];
		val[j] = '\0';
		Theader::THeaderDataItem *p = new Theader::THeaderDataItem;
		p->next = h.lst;
		p->name = strdup(name);
		h.lst = p;
		char *illegal1 = 0;
		long L = strtol(val,&illegal1,10);
		if (!will_be_string && *illegal1 == '\0' && L != LONG_MIN && L != LONG_MAX) {
			p->type = Theader::INT;
			p->i = L;
		} else {
			double d = strtod(val,&illegal1);
			if (!will_be_string && *illegal1 == '\0') {
				p->type = Theader::REAL;
				p->x = d;
			} else {
				p->type = Theader::STRING;
				p->s = strdup(val);
			}
		}
	}
	return o;
}

ostream& operator<<(ostream& o, const Theader::print& p)
{
	switch (p.t) {
	case Theader::INT:
		o << "int";
		break;
	case Theader::REAL:
		o << "real";
		break;
	case Theader::STRING:
		o << "string";
		break;
	}
	return o;
}

Theader::THeaderDataItem *Theader::find(const char *name) const
{
	THeaderDataItem *p;
	for (p=lst; p; p=p->next)
		if (!strcasecmp(p->name,name)) return p;
	return 0;
}

Theader::TItemType Theader::gettype(const char *name) const
{
	THeaderDataItem *p = find(name);
	if (!p) {
		cerr << "*** Theader::gettype: name '" << name << "' does not exist in header\n";
		return INT;
	}
	return p->type;
}

bool Theader::check(const char *name, TItemType t, bool probeflag) const
// If probeflag==true, t is ignored
{
	THeaderDataItem *p = find(name);
	if (probeflag) return (p != 0);
	if (!p) {
		cerr << "*** Theader::check: required name '" << name << "' does not appear in header\n";
		return false;
	}
	if (p->type != t && !(p->type == INT && t == REAL)) {
		cerr << "*** Theader::check: " << name << " should be " << print(t) << ", not " << print(p->type) << "\n";
		return false;
	}
	return true;
}

long Theader::getint(const char *name) const
{
	if (!check(name,INT)) return 0;
	THeaderDataItem *p = find(name);
	return p->i;
}

double Theader::getreal(const char *name) const
{
	if (!check(name,REAL)) return 0;
	THeaderDataItem *p = find(name);
	return (p->type == INT ? p->i : p->x);
}

char *Theader::getstr(const char *name) const
{
	if (!check(name,STRING)) return 0;
	THeaderDataItem *p = find(name);
	return p->s;
}

ostream& operator<<(ostream& o, const Theader& h)
{
	Theader::THeaderDataItem *p;
	for (p=h.lst; p; p=p->next) {
		o << "'" << p->name << "' = ";
		switch (p->type) {
		case Theader::INT:
			o << p->i << " (int)";
			break;
		case Theader::REAL:
			o << p->x << " (real)";
			break;
		case Theader::STRING:
			o << '"' << p->s << '"';
			break;
		}
		o << '\n';
	}
	return o;
}

Theader::Theader(const char *fn)
{
	lst = 0;
	ifstream f(fn);
	f >> *this;
}

Theader::~Theader()
{
	while (lst) {
		THeaderDataItem *p = lst;
		lst = lst->next;
		free(p->name);
		if (p->type == STRING) free(p->s);
		delete p;
	}
}
