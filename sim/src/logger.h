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

#ifndef LOGGER_H
#define LOGGER_H

#define MAX_COUNTERS 20

//! Log file writer
class Logger
{
private:
    static unsigned long int totalLogLines; //!< Number of total counted log lines
    static int lineHeaderChars; //!< Number of header chars in each line
    const char* logName; //!< Log file name
    std::fstream* logfile; //!< Log file stream
    unsigned long int logLines;
    unsigned long int maxLogLines;
    bool firstLine;
    bool maxLimitReached;
    bool doFileLimitAbort;
    bool headerAndLineNumbering;
    int currentCounter;
    std::string lineHeader;
    double loggerTime;
    double counter[MAX_COUNTERS];
    double counterInterval[MAX_COUNTERS];
    bool counterModeRepetitive[MAX_COUNTERS];
    bool counterModeLogarithmic[MAX_COUNTERS];
    void writeInit(const char*);
    void writeExit();
    void doChecks();
    void checkFirstLine();
    void checkFileLines();
    bool checkCounter();
    void checkTime();
    long int getFileSize();
    std::string getLineHeaderStr();
    std::string insertLineHeaders(std::string);
public:
    Logger(const char* filename = "logfile.log", const unsigned long int maxLines = 100000, const bool headerAndLineNumberingg = true);
    ~Logger();
    void init();
    char fill () const;
    char fill (char fillch);
    std::ios_base::fmtflags flags() const;
    std::ios_base::fmtflags flags(std::ios_base::fmtflags fmtfl);
    std::streamsize precision() const;
    std::streamsize precision(std::streamsize prec);
    std::ios_base::fmtflags setf(std::ios_base::fmtflags fmtfl);
    std::ios_base::fmtflags setf(std::ios_base::fmtflags fmtfl,std::ios_base::fmtflags mask);
    void unsetf(std::ios_base::fmtflags mask);
    std::streamsize width() const;
    std::streamsize width(std::streamsize wide);
    void flush();
    void setCounterInterval(const int);
    void zeroCounter();
    void zeroAllCounters();
    void useCounter(const int);
    int getCounter();
    double getCounterState();
    void writeCounterState();
    void increaseCounter(const int);
    void setDoFileLimitAbort(const bool);
    void setCounterModeRepetitive(const bool);
    void setCounterModeLogarithmic(const bool);
    Logger& operator<<(std::ostream& (*pf)(std::ostream&));
    Logger& operator<<(std::ios_base& (*pf)(std::ios_base&));
    Logger& operator<<(std::ostream& os);
    Logger& operator<<(const char);
    Logger& operator<<(const char*);
    Logger& operator<<(const std::string);
    Logger& operator<<(const float);
    Logger& operator<<(const double);
    Logger& operator<<(const int);
    Logger& operator<<(const unsigned long int);
};

#endif

