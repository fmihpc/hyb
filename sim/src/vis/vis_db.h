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

#ifndef VIS_DB_H
#define VIS_DB_H

#include <vector>
#include <string>

//! Base class of simulation database
class DB
{
public:
    /** \Precision of written binary numeric values.
     *
     * SAME: same precision than in program's internal presentation
     * FLOAT: precision of float type
     * DOUBLE: precision of double type
     */
    enum FloatingPrecision {SAME, FLOAT, DOUBLE};
    /** \brief Struct that is used to define which particles from all particles
     *         are selected.
     *
     * Allows selection of particle populations and particles in these
     * populations. Particle indexes are in arbitrary order.
     */
    struct Particles {
        bool allPopulations;           //! Select all populations
        unsigned int startPopID;       //! Smallest population ID that is selected
        unsigned int endPopID;         //! Biggest population id that is selected
        bool allParticles;             //! Select all particles
        std::size_t startParticleIdx;  //! Smallest particle index that is selected
        std::size_t endParticleIdx;    //! Biggest particle index that is selected
    };
};

//! Interface for visualization file writers
class VisDB : public DB
{
public:
    virtual ~VisDB() { }
    /** \brief Write visualization data to file with given name.
     *
     * \param filename Name of file to which visualization data is written
     * \param writeParticles Which particles (if any) are written to file
     */
    virtual void writeVisValues(const std::string& filename,
                                FloatingPrecision binFloatingPrec = FLOAT,
                                std::vector<Particles> writeParticles
                                = std::vector<Particles>(0)) = 0;
};

//! Interface for breakpoint file writers
class DumpDB : public DB
{
public:
    virtual ~DumpDB() { }
    virtual void dumpState(const std::string& filename) = 0;
    virtual void readState(const std::string& filename) = 0;
};

#endif

