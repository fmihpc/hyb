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

#ifndef VIS_DATA_SOURCE_H
#define VIS_DATA_SOURCE_H

#include <string>
#include <vector>
#include <map>
#include "cpp_utils.h"
#include "container.h"
#include "../definitions.h"
#include "../random.h"

/** \brief Interface from which visualization and dump data writers
 *         can read data that they should write.
 *
 * CommonDataSourdeData contains common parts of interface, and
 * its subclasses contain more specialized parts.
 */
struct CommonDataSourceData {
    typedef real* RealArray;
    struct Patch {
        //! Some data of one cell.
        struct Cell {
            /** \brief How many values cell vector variables have in this cell.
             *
             * Number of values can be bigger than 1, because sides
             * of cell can be refined (usually they are refined in dump data
             * and unrefined in visualization data). Index of outer sequence is index
             * of variable (same index as in sequence VisData::cellVectorVariables
             * or DumpData::cellVectorVariables)
             * and inner sequence is component of variable.
             */
            SequenceHandle<SequenceHandle<unsigned int> > nCellVectorVariableValues;
            /** \brief Like nCellVectorVariableValues, but for particle vector
             *              variables.
             */
            SequenceHandle<SequenceHandle<unsigned int> >
            nParticleVectorVariableValues;
        };
        /** \brief Refinement ratio of patch. N times bigger value means
         *         N times shorter edge of cell.
         */
        unsigned int refinementRatio;
        //! Start coordinates of patch (inclusive)
        SequenceHandle<size_t> startCoordinates;
        //! End coordinates of patch (inclusive)
        SequenceHandle<size_t> endCoordinates;
        /** \brief Some data of cells in patch.
         *
         * Order of patches in this sequence is defined by member
         * function flatIndex.
         */
        ConstSequenceHandle<Cell> cells;
        /** \brief Converts multidimensinal coordinates of cell to
         *         one dimensional index in sequence #cells.
         *
         * This function is currently untested.
         */
        std::size_t flatIndex
        (const std::vector<std::size_t>& cellCoordinates) const;
    };
    /** \brief Data of one vector variable.
     */
    struct VectorVariable {
        //! Data of one component of vector variable.
        struct Component {
            //! Name of component
            std::string name;
            /** \brief Values of component at each cell or particle.
             *
             * Index of outer sequence is index of patch (same index as
             * in sequence CommonDataSourceData::patches) if variable is
             * cell variable and 0 if variable is particle variable.
             * Index of inner sequence is index of cell in patch
             * (same index as in sequence CommonDataSourceData::Patch::cells)
             * if variable is cell variable or index of particle
             * (same index as in sequence VisData::particles or DumpData::particles)
             * if variable is particle variable.
             */
            SequenceHandle<RealArray> values;
        };
        //! name of variable
        std::string name;
        //! Whether there is also an expression defining this variable
        bool hasExpression;
        //! Components of variable
        SequenceHandle<Component> components;
    };
    /** \brief Data of one expression
     */
    struct Expression {
        std::string name;     //! Name of expression
        std::string formula;  //! Formula of expression
    };
    /** \brief Data of one particle
     */
    struct Particle {
        SequenceHandle<real> coordinates;   //! Coordinates of particle
    };
    /** \brief Coordinates of data in space.
     *
     * This is location of cell with coordinates (for example,
     * CommonDataSourceData::Cell::startCoordinates) of which all
     * components are 0.
     */
    ConstSequenceHandle<real> startCoordinates;
    /** \brief Length of edge of cell on refinement ratio 1.
     */
    ConstSequenceHandle<real> refRatio1CellSize;
    /** \brief AMR patches of data.
     *
     * These patches don't have "holes": they are areas of same refinement
     * ratio.
     */
    ConstSequenceHandle<Patch> patches;
    //! Current time in simulation
    real time;
};

//! Parts of data that are specialized to visualization.
struct VisData : public CommonDataSourceData {
    ConstSequenceHandle<Particle> particles;
    ConstSequenceHandle<VectorVariable> cellVectorVariables;
    ConstSequenceHandle<VectorVariable> particleVectorVariables;
    ConstSequenceHandle<Expression> expressions;
};

//! Parts of data that are specialized to breakpoints.
struct DumpData : public CommonDataSourceData {
    SequenceHandle<Particle> particles;
    SequenceHandle<VectorVariable> cellVectorVariables;
    SequenceHandle<VectorVariable> particleVectorVariables;
    int cnt_dt;
    Tportrand portrand;
};

//! Interface for components that produce visualization data.
class VisDataSource
{
public:
    virtual ~VisDataSource() { }
    //! Retuns visualization data
    virtual VisData getVisData() const = 0;
};

//! Interface for components that produce breakpoint data.
class DumpDataSource
{
public:
    virtual ~DumpDataSource() { }
    //! Returns breakpoint data
    virtual DumpData getDumpData() = 0;
};

inline std::size_t CommonDataSourceData::Patch::flatIndex
(const std::vector<std::size_t>& cellCoordinates) const
{
    std::size_t multiplier = 1, idx = 0;
    unsigned int sizeIdx = 0;
    for (std::vector<std::size_t>::const_iterator coord =
             cellCoordinates.begin();
         coord != cellCoordinates.end(); ++coord) {
        idx += *coord * multiplier;
        multiplier *= endCoordinates[sizeIdx] - startCoordinates[sizeIdx] + 1;
        ++sizeIdx;
    }
    return idx;
}

template<Dimension nDims>
Interval<nDims> getNodeIntervalFromPatch
(const CommonDataSourceData::Patch& patch, unsigned int coordRatio = 1,
 bool includeLast = true)
{
    Interval<nDims> i;
    for (Dimension dim = 0; dim < nDims; ++dim) {
        i.getBegin()[dim] = patch.startCoordinates[dim] * coordRatio;
        i.getEnd()[dim] = (patch.endCoordinates[dim] + (includeLast ? 1 : 0) )
                          * coordRatio + 1;
    }
    return i;
}

/** \brief Returns node coordinates from given data source without duplicates.
 *
 * If coordinates of data source are traversed, there are duplicate node
 * coordinates at adjacent sides of different patches.
 * This function returns all node coordinates of data source without
 * these duplicates.
 *
 * Coordinates are returned in std::map where keys are coordinates
 * and values are set to zero (caller can set this values for its
 * own purposes.
 *
 * Template argument nDims means dimensions of data source data.
 */
template <Dimension nDims>
SharedPtr<std::map<Coordinates<nDims>, std::size_t> >
getCoordinatesWithoutDuplicates(const CommonDataSourceData& data)
{
    SharedPtr<std::map<Coordinates<nDims>, std::size_t> > coordIndexes
    (new std::map<Coordinates<nDims>, std::size_t>);
    for (ConstSequenceHandle<CommonDataSourceData::Patch>::const_iterator patch
         = data.patches.begin(); patch != data.patches.end(); ++patch) {
        const unsigned int maxRefRatio
            = data.patches[data.patches.size() - 1].refinementRatio;
        const unsigned int coordRatio = maxRefRatio / patch->refinementRatio;
        Interval<nDims> patchInterval
            = getNodeIntervalFromPatch<nDims>(*patch);
        for (typename Interval<nDims>::const_iterator coords
             = patchInterval.begin(); coords != patchInterval.end(); ++coords) {
            (*coordIndexes)[(*coords) * coordRatio] = 0;
        }
    }
    return coordIndexes;
}

#endif

