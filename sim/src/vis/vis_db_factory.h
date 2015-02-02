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

#ifndef VIS_DB_FACTORY_H
#define VIS_DB_FACTORY_H

#include <stdexcept>
#include "../params.h"
#include "vis_data_source.h"
#include "vis_db.h"
#include "vis_db_vtk.h"

enum VisFormat { VTK_LEGACY_ASCII = 0, VTK_LEGACY_BINARY = 1 };

//! Visualization database factory
class VisDBFactory
{
public:
    /** \brief Creates and returns (dynamically allocated) instance
     *         of subclass of VisDB. Class is of type that writes format "format".
     */
    static VisDB& getVisDB(VisFormat format, const VisDataSource& dataSource) {
        switch(format) {
        case VTK_LEGACY_ASCII:
            return *new VTKVisDB(dataSource, VTKVisDB::LEGACY_ASCII,
                                 VTKVisDB::CALCULATE_ONCE);
        case VTK_LEGACY_BINARY:
            return *new VTKVisDB(dataSource, VTKVisDB::LEGACY_BINARY,
                                 VTKVisDB::CALCULATE_ONCE);
        }
        throw std::invalid_argument("Unrecognized file format");
    }
};

#endif

