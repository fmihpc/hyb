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

#include <string>
#include <vector>
#include <iostream>
#include <string.h>
#include <stdlib.h>
using namespace std;



inline void eatWhite( string &str)
{
    while(isspace(str[0]))
        str.erase(0,1);
}


inline void getToken(string &str, vector<string> &tokens)
{
    string part;
    while(!str.empty() and !isspace(str[0]))
    {
        part.push_back(str[0]);
        str.erase(0,1);
    }
    if(part.length()>0)
        tokens.push_back(part);
}


void tokenize(string line, vector<string> &tokens) 
{
    while(line.length())
    {
        eatWhite(line);
        getToken(line,tokens);
    }
}


void getDouble( string &str, double &val)
{
    char *end;
    val = strtod( str.c_str(), &end);
    if(end != NULL and end[0] != 0)
        cerr << "Error: non double argument. (Improve this message...)"  << endl, exit(-1);
}


void getLong( string &str, long int &val)
{
    char *end;
    val = strtol( str.c_str(), &end, 0);
    if(end != NULL and end[0] != 0)
        cerr << "Error: non integer argument. (Improve this message...)"  << endl, exit(-1);
}


bool isDouble(string  &str)
{
    char *end;
    strtod(str.c_str(), &end);
    if(end != NULL and end[0] != 0)
        return false;
    return true;
}


void itos(int n, string &str)
{
    int m = abs(n);
    do{
        div_t q = div(m,10);
        m=q.quot;
        str.insert(0,1,(char)q.rem+'0');
    }while(m!=0);

    if(n<0)
        str.insert(0,1,'-');
}

