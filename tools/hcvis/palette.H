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
#  pragma interface
#endif

class Tpalette {
private:
	int N;
	unsigned char *pal;
	void default_init();
public:
	Tpalette() {default_init();}
	void load(const char *fn);
	void grayscale() {delete [] pal; default_init();}
	void invert();
	void band(int n=10);
	void dim();		// multiply all by 0.85
	Tpalette(const char *fn) {default_init(); load(fn);}
	void RGB(float x, float& R, float& G, float& B) const;
	~Tpalette() {delete [] pal;}
};

inline void Tpalette::RGB(float x, float& R, float& G, float& B) const
{
	int i = (int)(x*N);
	if (i < 3) i = 3;
	if (i > N-1) i = N-1;
	R = pal[3*i  ]/256.0;
	G = pal[3*i+1]/256.0;
	B = pal[3*i+2]/256.0;
}
