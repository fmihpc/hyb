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

#ifndef VIS_DB_VTK_H
#define VIS_DB_VTK_H

#include <string>
#include "cpp_utils.h"
#include "vis_db.h"
#include "vis_data_source.h"

//! VTK visualization database
class VTKVisDB : public VisDB
{
public:
    enum FileFormat {LEGACY_ASCII, LEGACY_BINARY};
    /** \brief How to remove duplicate coordinates
     *
     * DONT_REMOVE: don't remove duplicate coordinates.
     * CALCULATE_ONCE: remove duplicate coordinates, and calculate which
     * coordinates are duplicates when writeVisValues is called first time,
     * and use this information in subsequent calls.
     * CALCULATE_ALWAYS: remove duplicate coordinates, and calculate which
     * coordinates are duplicates every time when writeVisValues is called
     */
    enum RemoveDuplicatesMode {DONT_REMOVE, CALCULATE_ONCE, CALCULATE_ALWAYS};
    VTKVisDB(const VisDataSource& dataSource, FileFormat format = LEGACY_ASCII,
             RemoveDuplicatesMode duplMode = DONT_REMOVE);
    virtual void writeVisValues(const std::string& filename,
                                FloatingPrecision binFloatingPrec = FLOAT,
                                std::vector<Particles> writeParticles
                                = std::vector<Particles>(0));
private:
    FileFormat m_format;
    RemoveDuplicatesMode m_duplMode;
    const VisDataSource& m_dataSource;
    SharedPtr<std::map<Coordinates<3>, std::size_t> > m_coordIndexes;
};

/*class VTKDumpDB : public DumpDB {
public:
	VTKDumpDB(CartesianDataSource& dataSource);
	void dumpState();
	void readState();
private:
	CartesianDataSource& mDataSource;
};*/

#endif

