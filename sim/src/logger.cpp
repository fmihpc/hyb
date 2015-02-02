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
#include <fstream>
#include <sstream>
#include <ctime>
#include <cmath>
#include "logger.h"
#include "params.h"

using namespace std;

unsigned long int Logger::totalLogLines = 1;
int Logger::lineHeaderChars = 14;

//! Constructor
Logger::Logger(const char* filename, const unsigned long int maxLines, const bool headerAndLineNumberingg) : logName(filename)
{
    maxLogLines = maxLines;
    logfile = new std::fstream(NULL,std::fstream::out);
    if (headerAndLineNumberingg == true) {
        headerAndLineNumbering = true;
    } else {
        headerAndLineNumbering = false;
    }
}

//! Destructor
Logger::~Logger()
{
    if (headerAndLineNumbering == true) {
        writeExit();
    }
    logfile->close();
    delete logfile;
}

//! Initialize logger
void Logger::init()
{
    logfile = new std::fstream(logName,std::fstream::out);
    if (headerAndLineNumbering == true) {
        writeInit(logName);
    }
    currentCounter = 0;
    logLines = 1;
    firstLine = true;
    maxLimitReached = false;
    doFileLimitAbort = false;
    lineHeader = ": ";
    loggerTime = -1;
    // Initialize counters
    for(int i = 0; i < MAX_COUNTERS; ++i) {
        counter[i] = 0;
        counterInterval[i] = 1;
        counterModeRepetitive[i] = true;
        if (i < 10) {
            counterModeLogarithmic[i] = false;
        }
        // Counters >= 10 are logarithmic by default
        else {
            counterModeLogarithmic[i] = true;
        }
    }
}

//! Write initial log lines
void Logger::writeInit(const char* filename)
{
    time_t rawtime;
    time(&rawtime);
    (*logfile) << "|-----------------------------------------------|\n";
    (*logfile) << "| File name: " << filename << "\n";
    (*logfile) << "| Line limit: " << maxLogLines << "\n";
    (*logfile) << "| Opened at " << ctime(&rawtime);
    (*logfile) << "| Logging: " << Params::codeVersion << "\n";
    (*logfile) << "|-----------------------------------------------|\n\n";
}

//! Write exit lines
void Logger::writeExit()
{
    time_t rawtime;
    time(&rawtime);
    (*logfile) << "\n\n";
    (*logfile) << "|-----------------------------------------------|\n";
    (*logfile) << "| File lines: ~" << std::flush << logLines << "\n";
    (*logfile) << "| File size: ~" << std::flush << ceil(getFileSize()/1024.0) << " kB\n";
    (*logfile) << "| Closed at " << ctime(&rawtime);
    (*logfile) << "|-----------------------------------------------|\n\n";
}

char Logger::fill() const
{
    return (*logfile).fill();
}

char Logger::fill(char fillch)
{
    return (*logfile).fill(fillch);
}

ios_base::fmtflags Logger::flags() const
{
    return (*logfile).flags();
}

ios_base::fmtflags Logger::flags(ios_base::fmtflags fmtfl)
{
    return (*logfile).flags(fmtfl);
}

streamsize Logger::precision() const
{
    return (*logfile).precision();
}

streamsize Logger::precision(streamsize prec)
{
    return (*logfile).precision(prec);
}

ios_base::fmtflags Logger::setf(ios_base::fmtflags fmtfl)
{
    return (*logfile).setf(fmtfl);
}

ios_base::fmtflags Logger::setf(ios_base::fmtflags fmtfl,ios_base::fmtflags mask)
{
    return (*logfile).setf(fmtfl,mask);
}

void Logger::unsetf(ios_base::fmtflags mask)
{
    (*logfile).unsetf(mask);
}

streamsize Logger::width() const
{
    return (*logfile).width();
}

streamsize Logger::width(streamsize wide)
{
    return (*logfile).width(wide);
}

//! Input operator for manipulators
Logger& Logger::operator<<(std::ostream& (*pf)(std::ostream& ))
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    // Check whether the stream manipulator is std::endl
    if (pf == static_cast<std::ostream& (*)(std::ostream&)>(std::endl)) {
        (*logfile) << insertLineHeaders("\n");
    } else {
        (*logfile) << pf;
    }
    (*logfile) << std::flush;
    return *this;
}

Logger& Logger::operator<<(std::ios_base& (*pf)(std::ios_base& ))
{
    (*logfile) << pf;
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for single character
Logger& Logger::operator<<(const char msg)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    // Construct a string from char
    string msgStr(1, msg);
    (*logfile) << insertLineHeaders(msgStr);
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for character array
Logger& Logger::operator<<(const char* msg)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    string msgStr(msg);
    (*logfile) << insertLineHeaders(msgStr);
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for string
Logger& Logger::operator<<(const string msg)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    (*logfile) << insertLineHeaders(msg);
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for float number
Logger& Logger::operator<<(const float val)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    (*logfile) << val;
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for double number
Logger& Logger::operator<<(const double val)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    (*logfile) << val;
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for unsigned long interger number
Logger& Logger::operator<<(const unsigned long int val)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    (*logfile) << val;
    (*logfile) << std::flush;
    return *this;
}

//! Input operator for integer number
Logger& Logger::operator<<(const int val)
{
    doChecks();
    if(checkCounter() == false) {
        return *this;
    }
    (*logfile) << val;
    (*logfile) << std::flush;
    return *this;
}

void Logger::flush()
{
    logfile->flush();
}

void Logger::doChecks()
{
    checkFirstLine();
    checkFileLines();
    checkTime();
}

void Logger::checkFirstLine()
{
    if(firstLine == true && headerAndLineNumbering == true) {
        (*logfile) << getLineHeaderStr();
        ++Logger::totalLogLines;
        ++logLines;
        firstLine = false;
    }
}

//! Check whether file line limit is reached and set file pointer to NULL or abort program execution
void Logger::checkFileLines()
{
    if(maxLimitReached == false && logLines > maxLogLines) {
        (*logfile) << "\n\nWARNING [Logger::checkFileLines]: logfile line limit (" <<  maxLogLines << ") reached.. ";
        if(doFileLimitAbort == true) {
            (*logfile) << "aborting program execution\n";
            writeExit();
            (*logfile) << std::flush;
            logfile->close();
            delete logfile;
            doabort();
        } else {
            (*logfile) << "setting file pointer to NULL\n";
            writeExit();
            (*logfile) << std::flush;
            logfile->close();
            delete logfile;
            logfile = new std::fstream(NULL,std::fstream::out);
            maxLimitReached = true;
        }
    }
}

//! Check current counter condition
bool Logger::checkCounter()
{
    if (counterModeLogarithmic[currentCounter] == false) {
        // Repetitive counter mode
        if (counterModeRepetitive[currentCounter] == true) {
            if( fmod(counter[currentCounter], counterInterval[currentCounter]) == 0.0 ) {
                return true;
            } else {
                return false;
            }
        }
        // Basic non-repetitive counter mode
        else {
            if(counter[currentCounter] < counterInterval[currentCounter]) {
                return true;
            } else {
                return false;
            }
        }
    }
    // Logarithmic counter mode
    else {
        if( fmod(counter[currentCounter], counterInterval[currentCounter]) == 0.0 ) {
            return true;
        } else {
            return false;
        }
    }
}

//! Check if the simulation time has changed since the last log entry and if it has changed write time&timestep to logfile
void Logger::checkTime()
{
    if(loggerTime != Params::t && checkCounter() == true && headerAndLineNumbering == true) {
        stringstream ss;
        ss << "|----------------------------------|\n";
        ss << "| Simulation time = " << Params::t << " s\n";
        ss << "| Timesteps taken = " << Params::cnt_dt << "\n";
        ss << "|----------------------------------|\n";
        (*logfile) << insertLineHeaders(ss.str());
        loggerTime = Params::t;
    }
}

//! Get current file size
long int Logger::getFileSize()
{
    ifstream filee(logName, ios::in | ios::binary);
    filee.seekg(0, ios::end);
    long int fileSize = filee.tellg();
    filee.close();
    return fileSize;
}

//! Return current line header as a string in the format: "0001: "
string Logger::getLineHeaderStr()
{
    int n = Logger::totalLogLines;
    stringstream ss;
    string lineNumberStr;
    ss << n;
    ss >> lineNumberStr;
    // lineHeaderChars when n < 10000:
    // "1234lineHeaderStr" => m = 4 + lineHeaderStr
    int m = 4 + lineHeader.size();
    Logger::lineHeaderChars = m;
    if(n < 10) {
        // "x" -> "000x"
        lineNumberStr.insert(0,"000");
    } else if(n < 100) {
        // "xy" -> "00xy"
        lineNumberStr.insert(0,"00");
    } else if(n < 1000) {
        // "xyz" -> "0xyz"
        lineNumberStr.insert(0,"0");
    } else if(n < 10000) {
        // xyza -> "xyxa"
    } else if(n < 100000) {
        Logger::lineHeaderChars = m+1;
    } else if(n < 1000000) {
        Logger::lineHeaderChars = m+2;
    } else if(n < 10000000) {
        Logger::lineHeaderChars = m+3;
    } else if(n < 100000000) {
        Logger::lineHeaderChars = m+4;
    }
    return lineNumberStr + lineHeader;
}

//! Inserts line header before each linebreak ('\n')
string Logger::insertLineHeaders(string str)
{
    if (headerAndLineNumbering == false) {
        return str;
    }
    string newStr = str;
    // Find the location of the first endl char
    string::size_type endlLoc = str.find('\n', 0);
    int n = 1;
    // Go through all '\n' chars and insert line numbers
    // in front of each new lines
    while(endlLoc != string::npos) {
        newStr.insert(endlLoc+n, getLineHeaderStr());
        n = n + Logger::lineHeaderChars;
        ++Logger::totalLogLines;
        ++logLines;
        endlLoc =  str.find('\n', endlLoc+1);
    }
    return newStr;
}

//! Sets counter interval for current counter
void Logger::setCounterInterval(const int N)
{
    if (N > 0) {
        counterInterval[currentCounter] = N;
    }
}

//! Set current counter to zero
void Logger::zeroCounter()
{
    counter[currentCounter] = 0;
    if (counterModeLogarithmic[currentCounter] == true) {
        counterInterval[currentCounter] = 1;
    }
}

//! Set all counters to zero
void Logger::zeroAllCounters()
{
    for(int i = 0; i < MAX_COUNTERS; ++i) {
        counter[i] = 0;
        if (counterModeLogarithmic[i] == true) {
            counterInterval[i] = 1;
        }
    }
}

//! Set current counter mode to repetitive (= write only every Nth line, where N = counter interval)
void Logger::setCounterModeRepetitive(const bool val)
{
    counterModeRepetitive[currentCounter] = val;
}

//! Set current counter mode to logarithmic (= write only every 10^Nth line, where N = 0,1,2..)
void Logger::setCounterModeLogarithmic(const bool val)
{
    counterModeLogarithmic[currentCounter] = val;
}

//! Set counter number n in use
void Logger::useCounter(const int n)
{
    if(n >= 0 && n < MAX_COUNTERS) {
        currentCounter = n;
    }
}

//! Get current counter in use
int Logger::getCounter()
{
    return currentCounter;
}

//! Get current counter state
double Logger::getCounterState()
{
    return counter[currentCounter];
}

//! Write counter start to logfile, e.g. (1e3)
void Logger::writeCounterState()
{
    doChecks();
    if(checkCounter() == true) {
        streamsize prec= precision();
        ios_base::fmtflags fmtfl = flags();
        precision(0);
        (*logfile) << "(" << getCounterState() << ")" << std::flush;
        precision(prec);
        flags(fmtfl);
    }
}

//! Increase the value of current counter
void Logger::increaseCounter(const int n)
{
    if(n > 0) {
        // Increase logarithmic interval
        if (counterModeLogarithmic[currentCounter] == true && checkCounter() == true) {
            counterInterval[currentCounter] *= 10;
        }
        counter[currentCounter] += n;
    }
}

//! Abort program execution when maximum file limit reached ( uses abort() )
void Logger::setDoFileLimitAbort(const bool val)
{
    doFileLimitAbort = val;
}

