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

#include <iostream>
#include <cstdlib>
#include "definitions.h"
#include "simulation.h"

using namespace std;

#ifdef LITTLE_ENDIAN
//! Byte conversion for little endian systems
void ByteConversion(int sz, unsigned char *x, int n)
{
    int i,ind,j;
    unsigned int b[16];
    if (sz > 16) {
        errorlog << "ERROR [ByteConversion]: internal error\n";
        doabort();
    }
    for (i=ind=0; i<n; i++,ind+=sz) {
        for (j=0; j<sz; j++) {
            b[j] = x[ind+j];
        }
        for (j=0; j<sz; j++) {
            x[ind+sz-1-j] = b[j];
        }
    }
}
#endif

//! Stop the program execution
void doabort()
{
    errorlog << "ABORT: doabort() called\n";
    cerr     << "ABORT: doabort() called\n";
    abort();
}

//! Handle the terminate signal
void TermHandler(int)
{
    static int cnt = 0;
    cnt++;
    errorlog << "STOPPING: kill or ctrl-c caught\n";
    cerr     << "STOPPING: kill or ctrl-c caught\n";
    Params::stoppingPhase = true;
    if(cnt > 1) {
        errorlog << "Signal caught for the 2nd time, forcing termination\n";
        cerr     << "Signal caught for the 2nd time, forcing termination\n";
        doabort();
    }
}

//! Integer to string using zeros in front (for example: int2string(75,4) = "0075")
string int2string(int nn, int zeros)
{
    string result;
    // Integer lengths
    int nLength = 0;
    if(abs(nn) < 10) {
        nLength = 1;
    } else if(abs(nn) < 100) {
        nLength = 2;
    } else if(abs(nn) < 1000) {
        nLength = 3;
    } else if(abs(nn) < 10000) {
        nLength = 4;
    } else if(abs(nn) < 100000) {
        nLength = 5;
    } else if(abs(nn) < 1000000) {
        nLength = 6;
    }
    stringstream ss;
    ss << nn;
    if(zeros <= nLength || nLength <= 0) {
        result = ss.str();
    } else {
        string addZeros(zeros - nLength,'0');
        result = addZeros + ss.str();
    }
    return result;
}

//! String to double
real string2double(string str)
{
    // Check the string includes only numeric characters
    const unsigned int N = 15;
    string xx[N] = {"0","1","2","3","4","5","6","7","8","9","e","E",".","+","-"};
    for (unsigned int i=0; i < str.length(); ++i) {
        string sub = str.substr(i,1);
        bool valueOk = false;
        for (unsigned int j=0; j < N; ++j) {
            if(sub.compare(xx[j]) == 0) {
                valueOk = true;
            }
        }
        if (valueOk == false) {
            errorlog << "ERROR [string2double(" << str << ")]: not a valid number\n";
            doabort();
        }
    }
    return atof( str.c_str() );
}

//! Returns the time (in seconds) elapsed since to the first call of this function.
real getExecutionSecs()
{
    //! Timestamp of the first call of this function.
    static time_t start = time(&start);
    time_t end;
    time (&end);
    return difftime(end,start);
}

//! Returns the time (in seconds) difference between the last two calls of this function.
real getLastIntervalSecs()
{
    static time_t start;
    static time_t end;
    static bool started = false;
    real t;
    if(started == false) {
        time(&start);
        t = 0.0;
        started = true;
    } else {
        //! Time now
        time(&end);
        //! Time difference since the last call of this function.
        t = difftime(end,start);
        // Start counting from this call.
        time(&start);
    }
    return t;
}

//! Returns a string representation (XXX hrs YY mins ZZ secs) of the given seconds
string secsToTimeStr(real t)
{
    double intpart;
    modf(t/3600.0,&intpart);
    int hrs  = static_cast<int>(intpart);
    modf((t - hrs*3600.0) / 60.0,&intpart);
    int mins = static_cast<int>(intpart);
    modf(t - hrs*3600.0 - mins*60.0,&intpart);
    int secs = static_cast<int>(intpart);
    stringstream ss;
    // hour string
    if (hrs < 10) {
        ss << "000";
    } else if (hrs < 100) {
        ss << "00";
    } else if (hrs < 1000) {
        ss << "0";
    } else {
        ss << "";
    }
    ss << hrs << " hrs ";
    // minute string
    if (mins < 10) {
        ss << "0";
    }
    ss << mins << " mins ";
    // second string
    if (secs < 10) {
        ss << "0";
    }
    ss << secs << " secs";
    return ss.str();
}

//! Remove spaces from the beginning and end of a string
string cropPrecedingAndTrailingSpaces(string str)
{
    // Crop preceding/trailing " ", "\n" or "\t"
    const unsigned int nStart = str.find_first_not_of(" \n\t");
    const unsigned int nEnd = str.find_last_not_of(" \n\t");
    // Almost empty string
    if(nStart >= nEnd-1) {
        return str;
    }
    str = str.substr(nStart,nEnd-nStart+1);
    // Replace "\t" -> " "
    string::size_type nn = str.find_first_of("\t");
    while (nn != string::npos) {
        string temp = str.substr(nn,1);
        if(temp.compare("\t") == 0) {
            str.replace(nn,1," ");
        }
        nn = str.find_first_of("\t",nn);
    }
    // Crop extra spaces: "  "  -> " "
    nn = str.find_first_of(" ");
    while (nn != string::npos) {
        string temp = str.substr(nn+1,1);
        while(temp.compare(" ") == 0) {
            str.erase(nn,1);
            temp = str.substr(nn+1,1);
        }
        nn = str.find_first_of(" ",nn+1);
    }
    // " \n" -> "\n"
    nn = str.find_first_of(" ");
    while (nn != string::npos) {
        string temp = str.substr(nn+1,1);
        while(temp.compare("\n") == 0) {
            str.erase(nn,1);
            temp = str.substr(nn+1,1);
        }
        nn = str.find_first_of(" ",nn+1);
    }
    // Remove "\n\n", "\n "
    nn = str.find_first_of("\n");
    while (nn != string::npos) {
        string temp = str.substr(nn+1,1);
        while(temp.compare(" ") == 0 || temp.compare("\n") == 0) {
            str.erase(nn+1,1);
            temp = str.substr(nn+1,1);
        }
        nn = str.find_first_of("\n",nn+1);
    }
    return str;
}

//! Skip braced content (handles also inner braces)
void skipBracedContent(ifstream &fileStream)
{
    const unsigned int maxCh = 1024;
    char tempA[maxCh];
    // Get characters before the next opening brace and place pointer
    fileStream.getline(tempA,maxCh,'{');
    fileStream.unget();
    // Check there are no illegal characters
    string str(tempA);
    if (str.find_first_not_of(" \n\t") != string::npos) {
        errorlog << "ERROR [skipBracedContent]: bad character before the opening brace (" << str << ")\n";
        doabort();
    } else if (fileStream.eof() == true) {
        errorlog << "ERROR [skipBracedContent]: no opening brace found\n";
        doabort();
    }
    // Go thru the characters
    int cnt = 0;
    bool loop = true;
    while(loop == true) {
        char ch = fileStream.get();
        string str(1,ch);
        if(str.compare("{") == 0) {
            ++cnt;
        } else if(str.compare("}") == 0) {
            --cnt;
        }
        if(cnt == 0) {
            loop = false;
        } else if(fileStream.eof() == true) {
            errorlog << "ERROR [skipBracedContent]: braces should appear in pairs\n";
            doabort();
        }
    }
}

