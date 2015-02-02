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

#include <stdexcept>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include "cpp_utils.h"
#include "container.h"
#include "vis_db_vtk.h"
#include "../definitions.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "../transformations.h"
#endif

using namespace std;

namespace
{

//! Writes tokens to stream.
class DataWriter
{
public:
    /** \brief How to write tokens.
     *
     * FORMATTED:        formatted write with operator <<
     * BINARY:           unformatted write
     * BYTECONV_BINARY:  unformatted write, convert endianness
     */
    enum WriteMode {FORMATTED, BINARY, BYTECONV_BINARY};
    DataWriter(ostream& os, WriteMode mode = FORMATTED,
               bool useTokenSeparator = false,
               std::string tokenSeparator = " ", unsigned int rowLength = 0,
               DB::FloatingPrecision binFloatingPrec = DB::FLOAT)
        : m_mode(mode), m_useTokenSeparator(useTokenSeparator),
          m_tokenSeparator(tokenSeparator), m_rowLength(rowLength),
          m_os(os), m_tokensOnCurrentRow(0),
          m_binFloatingPrec(binFloatingPrec) { }
    template <class T>
    DataWriter& write(const T& data);
    void setMode(WriteMode mode) {
        m_mode = mode;
    }
    void setUseTokenSeparator(bool useTokenSeparator) {
        m_useTokenSeparator = useTokenSeparator;
    }
    void setRowLength(unsigned int rowLength) {
        m_rowLength = rowLength;
        m_tokensOnCurrentRow = 0;
    }
private:
    WriteMode m_mode;
    bool m_useTokenSeparator;
    string m_tokenSeparator;
    unsigned int m_rowLength;
    ostream& m_os;
    unsigned int m_tokensOnCurrentRow;
    DB::FloatingPrecision m_binFloatingPrec;
    template <class T>
    void writeFloating(const T& data) {
        switch (m_mode) {
        case FORMATTED:
            writeWithMode(data);
            break;
        case BINARY:
        case BYTECONV_BINARY:
            switch (m_binFloatingPrec) {
            case DB::SAME:
                writeWithMode(data);
                break;
            case DB::DOUBLE:
                writeWithMode(static_cast<double>(data));
                break;
            case DB::FLOAT:
                writeWithMode(static_cast<float>(data));
                break;
            }
        }
    }
    template <class T>
    void writeWithMode(const T& data) {
        switch(m_mode) {
        case FORMATTED:
            writeFormatted(data);
            break;
        case BINARY:
            writeBinary(data);
            break;
        case BYTECONV_BINARY:
            writeByteConvBinary(data);
            break;
        }
    }
    template <class T>
    void writeFormatted(const T& data) {
        m_os << data;
        if (m_rowLength != 0 && ++m_tokensOnCurrentRow >= m_rowLength) {
            m_os << "\n";
            m_tokensOnCurrentRow = 0;
        } else if (m_useTokenSeparator) {
            m_os << m_tokenSeparator;
        }
    }
    template <class T>
    void writeBinary(const T& data) {
        m_os.write(reinterpret_cast<const char*>(&data), sizeof(T));
    }
    template <class T>
    void writeByteConvBinary(const T& data) {
        T wData = data;
        ByteConversion(sizeof(T), reinterpret_cast<unsigned char*>(&wData), 1);
        m_os.write(reinterpret_cast<const char*>(&wData), sizeof(T));
    }
};

template <class T>
DataWriter& DataWriter::write(const T& data)
{
    writeWithMode(data);
    return *this;
}

template <>
DataWriter& DataWriter::write(const double& data)
{
    writeFloating(data);
    return *this;
}

template <>
DataWriter& DataWriter::write(const float& data)
{
    writeFloating(data);
    return *this;
}

string getVTKType(DB::FloatingPrecision binFloatingPrec)
{
    switch (binFloatingPrec) {
    case DB::SAME:
        return "double";
    case DB::DOUBLE:
        return "double";
    case DB::FLOAT:
        return "float";
    }
    return "";
}


vector<std::size_t> getPatchDimensions(const VisData::Patch& p)
{
    vector<std::size_t> dims(p.startCoordinates.size());
    for (unsigned int dim = 0; dim < dims.size(); ++dim)
        dims[dim] = p.endCoordinates[dim] - p.startCoordinates[dim] + 1;
    return dims;
}

/** \brief Writes cells as structured points in VTK format.
 *
 * Cannot be used with AMR mesh.
 */
void writeCellsAsStructuredPoints(VisData& data, ostream& os,
                                  VTKVisDB::FileFormat format,
                                  DB::FloatingPrecision binFloatingPrec)
{
    DataWriter::WriteMode dataWriteMode;
    std::string dataFormatString;
    std::string dataTypeString;
    switch(format) {
    case VTKVisDB::LEGACY_ASCII:
        dataWriteMode = DataWriter::FORMATTED;
        dataFormatString = "ASCII";
        dataTypeString = "double";
        break;
    case VTKVisDB::LEGACY_BINARY:
        dataWriteMode = DataWriter::BYTECONV_BINARY;
        dataFormatString = "BINARY";
        dataTypeString = getVTKType(binFloatingPrec);
        break;
    default:
        throw invalid_argument("Unrecognized file format");
    }
    DataWriter writer(os, dataWriteMode, true, " ", 10, binFloatingPrec);
    vector<std::size_t> dims = getPatchDimensions(data.patches[0]);
    const CommonDataSourceData::Patch& patch = data.patches[0];
    real cellSize[] = { data.refRatio1CellSize[0] / patch.refinementRatio,
                        data.refRatio1CellSize[1] / patch.refinementRatio,
                        data.refRatio1CellSize[2] / patch.refinementRatio
                      };
    std::size_t nCells = dims[0] * dims[1] * dims[2];
    os << "# vtk DataFile Version 2.0\n";
    os << "Hybridcode visualization data\n";
    os << dataFormatString << "\n";
    os << "DATASET STRUCTURED_POINTS\n";
    os << "DIMENSIONS " << (dims[0] + 1) << " " << (dims[1] + 1)
       << " " << (dims[2] + 1) << "\n";
    os << "ORIGIN "
       << data.startCoordinates[0] + patch.startCoordinates[0] * cellSize[0]
       << " "
       << data.startCoordinates[1] + patch.startCoordinates[1] * cellSize[1]
       << " "
       << data.startCoordinates[2] + patch.startCoordinates[2] * cellSize[2]
       << "\n";
    os << "SPACING " << cellSize[0] << " " << cellSize[1] << " " << cellSize[2]
       << "\n";
    os << "CELL_DATA " << nCells;
    for (ConstSequenceHandle<VisData::VectorVariable>::const_iterator var
         = data.cellVectorVariables.begin();
         var != data.cellVectorVariables.end(); ++var) {
        os << "\nSCALARS " << var->name << " " << dataTypeString << " "
           << var->components.size() << "\n";
        os << "LOOKUP_TABLE default\n";
        for (std::size_t valIdx = 0; valIdx < nCells; ++valIdx) {
            for (SequenceHandle<VisData::VectorVariable::Component>::const_iterator
                 comp = var->components.begin(); comp != var->components.end();
                 ++comp) {
                writer.write(comp->values[0][valIdx]);
            }
        }
    }
}

/** \brief Writes cells as unstructured grid in VTK format.
 *
 * Can be used with AMR mesh.
 */
void writeCellsAsUnstructuredGrid
(VisData& data, ostream& os, VTKVisDB::FileFormat format,
 SharedPtr<std::map<Coordinates<3>, std::size_t> > coordIndexes,
 DB::FloatingPrecision binFloatingPrec)
{
    DataWriter::WriteMode dataWriteMode;
    std::string dataFormatString;
    std::string dataTypeString;
    switch(format) {
    case VTKVisDB::LEGACY_ASCII:
        dataWriteMode = DataWriter::FORMATTED;
        dataFormatString = "ASCII";
        dataTypeString = "double";
        break;
    case VTKVisDB::LEGACY_BINARY:
        dataWriteMode = DataWriter::BYTECONV_BINARY;
        dataFormatString = "BINARY";
        dataTypeString = getVTKType(binFloatingPrec);
        break;
    default:
        throw invalid_argument("Unrecognized file format");
    }
    DataWriter writer(os, dataWriteMode, true, " ", 10, binFloatingPrec);
    unsigned int maxRefRatio
        = data.patches[data.patches.size() - 1].refinementRatio;
    os << "# vtk DataFile Version 2.0\n";
    os << "Hybridcode visualization data\n";
    os << dataFormatString << "\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
    // Write cell node points
    unsigned int nPoints = 0, nCells = 0;
    for (ConstSequenceHandle<VisData::Patch>::const_iterator patch
         = data.patches.begin(); patch != data.patches.end(); ++patch) {
        vector<std::size_t> dims = getPatchDimensions(*patch);
        nCells += dims[0] * dims[1] * dims[2];
    }
    real spacing[] = { data.refRatio1CellSize[0] / maxRefRatio,
                       data.refRatio1CellSize[1] / maxRefRatio,
                       data.refRatio1CellSize[2] / maxRefRatio
                     };
    if (coordIndexes                 // don't remove duplicate coordinates
        == SharedPtr<std::map<Coordinates<3>, std::size_t> >(0)) {
        for (ConstSequenceHandle<VisData::Patch>::const_iterator patch
             = data.patches.begin(); patch != data.patches.end(); ++patch) {
            vector<std::size_t> dims = getPatchDimensions(*patch);
            nPoints += (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1);
        }
        os << "POINTS " << nPoints << " " << dataTypeString << "\n";
        for (ConstSequenceHandle<VisData::Patch>::const_iterator patch
             = data.patches.begin(); patch != data.patches.end(); ++patch) {
            vector<std::size_t> dims = getPatchDimensions(*patch);
            const unsigned int coordRatio = maxRefRatio / patch->refinementRatio;
            for (std::size_t z = patch->startCoordinates[2] * coordRatio;
                 z <= (patch->endCoordinates[2] + 1) * coordRatio; z += coordRatio)
                for (std::size_t y = patch->startCoordinates[1] * coordRatio;
                     y <= (patch->endCoordinates[1] + 1) * coordRatio;
                     y += coordRatio)
                    for (std::size_t x = patch->startCoordinates[0] * coordRatio;
                         x <= (patch->endCoordinates[0] + 1) * coordRatio;
                         x += coordRatio) {
                        // Spherical coordinates test for VisIt
                        /*
                        double xx = x * spacing[0] + data.startCoordinates[0];
                        double yy = y * spacing[1] + data.startCoordinates[1];
                        double zz = z * spacing[2] + data.startCoordinates[2];
                        double r = -1*xx + Params::box_xmax;
                        double theta = M_PI*(yy - Params::box_ymin)/Params::box_Y;
                        double phi = 2*M_PI*(zz - Params::box_zmin)/Params::box_Z;
                        writer.write(r*sin(theta)*cos(phi));
                        writer.write(r*sin(theta)*sin(phi));
                        writer.write(r*cos(theta));
                        */
                        gridreal crd[3];
                        crd[0] = x*spacing[0] + data.startCoordinates[0];
                        crd[1] = y*spacing[1] + data.startCoordinates[1];
                        crd[2] = z*spacing[2] + data.startCoordinates[2];
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
                        sph_transf_H2S_VTK(crd);
#endif
                        writer.write(crd[0]);
                        writer.write(crd[1]);
                        writer.write(crd[2]);
                    }
        }
    } else {           // remove duplicate coordinates
        nPoints = coordIndexes->size();
        os << "POINTS " << nPoints << " " << dataTypeString << "\n";
        std::size_t idx = 0;
        for (std::map<Coordinates<3>, std::size_t>::iterator coords
             = coordIndexes->begin(); coords != coordIndexes->end(); ++coords) {
            coords->second = idx;
            // Spherical coordinates test for VisIt
            /*
            double xx = coords->first[0] * spacing[0] + data.startCoordinates[0];
            double yy = coords->first[1] * spacing[1] + data.startCoordinates[1];
            double zz = coords->first[2] * spacing[2] + data.startCoordinates[2];
            double r = -1*xx + Params::box_xmax;
            double theta = M_PI*(yy - Params::box_ymin)/Params::box_Y;
            double phi = 2*M_PI*(zz - Params::box_zmin)/Params::box_Z;
            writer.write(r*sin(theta)*cos(phi));
            writer.write(r*sin(theta)*sin(phi));
            writer.write(r*cos(theta));
            */
            gridreal crd[3];
            for(int i=0; i<3; ++i) {
                crd[i] = coords->first[i] * spacing[i] + data.startCoordinates[i];
            }
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
            sph_transf_H2S_VTK(crd);
#endif
            writer.write(crd[0]);
            writer.write(crd[1]);
            writer.write(crd[2]);
            ++idx;
        }
    }
    os << "\n";
    // Write "grid structure"
    os << "CELLS " << nCells << " " << nCells * 9 << "\n";
    writer.setRowLength(9);
    if (coordIndexes            // duplicate coordinates aren't removed
        == SharedPtr<std::map<Coordinates<3>, std::size_t> >(0)) {
        unsigned int firstPatchPointIdx = 0;
        for (ConstSequenceHandle<VisData::Patch>::const_iterator patch
             = data.patches.begin(); patch != data.patches.end(); ++patch) {
            vector<std::size_t> dims = getPatchDimensions(*patch);
            int xJump = 1, yJump = dims[0] + 1,
                zJump = (dims[0] + 1) * (dims[1] + 1);
            for (std::size_t z = 0; z < dims[2]; ++z)
                for (std::size_t y = 0; y < dims[1]; ++y)
                    for (std::size_t x = 0; x < dims[0]; ++x) {
                        std::size_t startIdx = x + y * (dims[0] + 1) +
                                               z * (dims[0] + 1) * (dims[1] + 1);
                        int cornerIdx = startIdx + firstPatchPointIdx;
// Swapped the index walk for hexahedron configuration: 3<->4, 6<->7
// the original in comments
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
                        writer.write<int>(int(8))
                        .write<int>(cornerIdx)
                        .write<int>(cornerIdx + xJump)
                        .write<int>(cornerIdx + xJump + yJump)
                        .write<int>(cornerIdx + yJump)
                        .write<int>(cornerIdx + zJump)
                        .write<int>(cornerIdx + zJump + xJump)
                        .write<int>(cornerIdx + xJump + yJump + zJump)
                        .write<int>(cornerIdx + zJump + yJump);
#else
                        writer.write<int>(int(8))
                        .write<int>(cornerIdx)
                        .write<int>(cornerIdx + xJump)
                        .write<int>(cornerIdx + yJump)
                        .write<int>(cornerIdx + xJump + yJump)
                        .write<int>(cornerIdx + zJump)
                        .write<int>(cornerIdx + zJump + xJump)
                        .write<int>(cornerIdx + zJump + yJump)
                        .write<int>(cornerIdx + xJump + yJump + zJump);
#endif
                    }
            firstPatchPointIdx += (dims[0] + 1) * (dims[1] + 1) * (dims[2] + 1);
        }
    } else {                     // duplicate coordinates are removed
        for (ConstSequenceHandle<VisData::Patch>::const_iterator patch
             = data.patches.begin(); patch != data.patches.end(); ++patch) {
            const unsigned int coordRatio = maxRefRatio / patch->refinementRatio;
            Interval<3> patchInterval
                = getNodeIntervalFromPatch<3>(*patch, 1, false);
            for (Interval<3>::const_iterator coords
                 = patchInterval.begin(); coords != patchInterval.end();
                 ++coords) {
                Coordinates<3> co = (*coords) * coordRatio;
// Swapped the index walk for hexahedron configuration: 3<->4, 6<->7
// the original in comments
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
                writer.write<int>(int(8))
                .write<int>((*coordIndexes)
                            [co])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, 0, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, coordRatio, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, coordRatio, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, 0, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, 0, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, coordRatio, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, coordRatio, coordRatio)]);
#else
                writer.write<int>(int(8))
                .write<int>((*coordIndexes)
                            [co])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, 0, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, coordRatio, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, coordRatio, 0)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, 0, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, 0, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(0, coordRatio, coordRatio)])
                .write<int>((*coordIndexes)
                            [co + Coordinates<3>(coordRatio, coordRatio,
                                                 coordRatio)]);
#endif
            }
        }
    }
    os << "\n";
    os << "CELL_TYPES " << nCells << "\n";
    writer.setRowLength(1);
    // here convert to cell type VTK_WEDGE (13) when at pole lines, to VTK_HEXAHEDRON (12) otherwise
    for (std::size_t cellIdx = 0; cellIdx < nCells; ++cellIdx)
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        writer.write<int>(int(12));
#else
        writer.write<int>(int(11)); //old VTK_VOXEL (11) commented
#endif
    os << "\n";
    // Write variables
    os << "CELL_DATA " << nCells;
    writer.setRowLength(10);
    for (ConstSequenceHandle<VisData::VectorVariable>::const_iterator var
         = data.cellVectorVariables.begin();
         var != data.cellVectorVariables.end(); ++var) {
        os << "\nSCALARS " << var->name << " " << dataTypeString << " "
           << var->components.size() << "\n";
        os << "LOOKUP_TABLE default\n";
        for (unsigned int patchIdx = 0;
             patchIdx < var->components[0].values.size(); ++patchIdx)
            for (std::size_t valIdx = 0;
                 valIdx < data.patches[patchIdx].cells.size(); ++valIdx) {
                for (SequenceHandle<VisData::VectorVariable::Component>::
                     const_iterator comp = var->components.begin();
                     comp != var->components.end(); ++comp)
                    writer.write(comp->values[patchIdx][valIdx]);
            }
    }
    os << "\n";
}

//! Writes particles as unstructured grid in VTK format
void writeParticlesAsUnstructuredGrid(VisData& data, ostream& os,
                                      const VTKVisDB::Particles& particles)
{
    os << "# vtk DataFile Version 2.0\n";
    os << "Hybridcode visualization data\n";
    os << "ASCII\n";
    os << "DATASET UNSTRUCTURED_GRID\n";
    std::size_t nParticles = data.particles.size();
    // Write particle coordinates
    os << "POINTS " << nParticles << " float" << "\n";
    const unsigned int rowLength = 10;
    unsigned int rowPosition = 0;
    for (unsigned int i = 0; i < nParticles; ++i) {
        os << data.particles[i].coordinates[0] << " "
           << data.particles[i].coordinates[1] << " "
           << data.particles[i].coordinates[2];
        ++rowPosition;
        if (rowPosition < rowLength) {
            os << "  ";
        } else {
            os << "\n";
            rowPosition = 0;
        }
    }
    // Write "grid structure"
    os << "\n";
    os << "CELLS " << nParticles << " " << nParticles * 2 << "\n";
    rowPosition = 0;
    for (unsigned int i = 0; i < nParticles; ++i) {
        os << "1 " << i;
        ++rowPosition;
        if (rowPosition < rowLength) {
            os << "  ";
        } else {
            os << "\n";
            rowPosition = 0;
        }
    }
    os << "CELL_TYPES " << nParticles << "\n";
    rowPosition = 0;
    for (unsigned int i = 0; i < nParticles; ++i) {
        os << "2";
        ++rowPosition;
        if (rowPosition < rowLength) {
            os << " ";
        } else {
            os << "\n";
            rowPosition = 0;
        }
    }
    // Write variables
    os << "CELL_DATA " << nParticles << "\n";
    for (ConstSequenceHandle<VisData::VectorVariable>::const_iterator var
         = data.particleVectorVariables.begin();
         var != data.particleVectorVariables.end(); ++var) {
        os << "SCALARS " << var->name << " float "
           << var->components.size() << "\n";
        os << "LOOKUP_TABLE default\n";
        for (std::size_t valIdx = 0; valIdx < nParticles; ++valIdx) {
            for (SequenceHandle<VisData::VectorVariable::Component>::const_iterator
                 comp = var->components.begin(); comp != var->components.end();
                 ++comp) {
                os << comp->values[0][valIdx] << " ";
            }
            os << "\n";
        }
    }
}

}

VTKVisDB::VTKVisDB(const VisDataSource& dataSource, FileFormat format,
                   RemoveDuplicatesMode duplMode)
    : m_format(format), m_duplMode(duplMode), m_dataSource(dataSource),
      m_coordIndexes(0) { }

void VTKVisDB::writeVisValues(const std::string& filename,
                              FloatingPrecision binFloatingPrec,
                              vector<Particles> writeParticles)
{
    VisData data = m_dataSource.getVisData();
    // Remove duplicate coordinates, if requested
    switch(m_duplMode) {
    case DONT_REMOVE:
        m_coordIndexes = SharedPtr<std::map<Coordinates<3>, std::size_t> >(0);
        break;
    case CALCULATE_ONCE:
        if (m_coordIndexes == SharedPtr<std::map<Coordinates<3>, std::size_t> >(0))
            m_coordIndexes = getCoordinatesWithoutDuplicates<3>(data);
        break;
    case CALCULATE_ALWAYS:
        m_coordIndexes = getCoordinatesWithoutDuplicates<3>(data);
        break;
    }
    string fName = filename + ".vtk";
    ofstream of(fName.c_str());
    // If data has only one AMR patch, write as structured points
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    if (data.patches.size() == 1) {
        writeCellsAsStructuredPoints(data, of, m_format, binFloatingPrec);
    } else {
        writeCellsAsUnstructuredGrid(data, of, m_format, m_coordIndexes, binFloatingPrec);
    }
#else
    writeCellsAsUnstructuredGrid(data, of, m_format, m_coordIndexes, binFloatingPrec);
#endif
    of.close();
    // write particles
    for (vector<Particles>::const_iterator particles = writeParticles.begin();
         particles != writeParticles.end(); ++particles) {
        ostringstream pFilename;
        pFilename << filename << "_particles_popIDs_";
        if (particles->allPopulations)
            pFilename << "all_";
        else
            pFilename << particles->startPopID << "-" << particles->endPopID
                      << "_";
        pFilename << "pIdxs_";
        if (particles->allParticles)
            pFilename << "all";
        else
            pFilename << particles->startParticleIdx << "-"
                      << particles->endParticleIdx;
        pFilename << ".vtk";
        ofstream pof(pFilename.str().c_str());
        writeParticlesAsUnstructuredGrid(data, pof, *particles);
        pof.close();
        // save memory by removing unneeded reference
        if (m_duplMode == CALCULATE_ALWAYS)
            m_coordIndexes = SharedPtr<std::map<Coordinates<3>, std::size_t> >(0);
    }
}

