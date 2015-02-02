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

#include <cassert>
#include <iostream>
#include <sstream>
#include <algorithm>
#include "cpp_utils.h"
#include "container.h"
#include "vis_data_source_simulation.h"
#include "../grid.h"
#include "../magneticfield.h"
#include "../templates.h"
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
#include "../transformations.h"
#endif

using std::string;
using std::vector;

class GridPrivatesHelper
{
protected:
    typedef Tgrid::Tcell Tcell;
};

namespace
{

//! Class that allows other classes to use some private members of Tgrid
class GridPriv : private GridPrivatesHelper
{
public:
    typedef GridPrivatesHelper::Tcell Tcell;
};

//! Replaces occurrences of string "repl" in string "s" by string "by"
string replace(const string& s, const string& repl, const string& by)
{
    string c = s;
    std::size_t idx = 0;
    while (true) {
        std::size_t occIdx = c.find(repl, idx);
        if (occIdx == string::npos)
            break;
        c.replace(occIdx, repl.size(), by);
        idx = occIdx + by.size();
    }
    return c;
}

/** \brief Creates visualization or dump data filename (without suffix)
 *         from current time.
 */
string createTimeFilename()
{
    std::ostringstream oss;
    oss << static_cast<long int>(1000*Params::t + 0.5);
    return oss.str();
}

//! Return 2^power
inline unsigned int twoPower(unsigned int power)
{
    return 1 << power;
}

enum CoordinateID {X = 0, Y = 1, Z = 2};

/** \brief Returns an Interval describing surface of cube with given
 *         start corner ("begin"), distance from that corner and direction
 *         of surface normal.
 */
Interval<3> cubeSurface(const Coordinates<3>& begin,
                        std::size_t distance,
                        CoordinateID direction)
{
    return Interval<3>(Coordinates<3>
                       (begin[0] + (direction == X ? distance : 0),
                        begin[1] + (direction == Y ? distance : 0),
                        begin[2] + (direction == Z ? distance : 0)),
                       Coordinates<3>
                       (begin[0] + distance + 1,
                        begin[1] + distance + 1,
                        begin[2] + distance + 1));
}


//! Struct that contains formulas for cell variables
struct CellFormulas {
    struct B1Formula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(0.5*(cell.faceave(comp, 0, Tgrid::FACEDATA_B) +
                                       cell.faceave(comp, 1, Tgrid::FACEDATA_B)));
            return results;
        }
    };
    struct B1AveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(0.5*(cell.faceave(comp, 0, Tgrid::FACEDATA_AVEB) +
                                       cell.faceave(comp, 1, Tgrid::FACEDATA_AVEB)));
            return results;
        }
    };
    struct B0Formula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            real B0[3] = {0.0, 0.0, 0.0};
            addConstantMagneticField(cell.centroid, B0);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(B0[comp]);
            return results;
        }
    };
    struct vFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results(3);
            cell.plist.calc_avev(results[0], results[1],
                                 results[2], popId);
            return results;
        }
    };
    struct nFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.calc_weight(popId)/(cell.size*cell.size*cell.size));
        }
    };
    struct nTotAveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.ave_nc);
        }
    };
#ifdef SAVE_POPULATION_AVERAGES
    struct vAveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results(3);
            real n=0,vx=0,vy=0,vz=0;
            for(unsigned int i=0; i< popId.size(); i++) {
                int j = popId[i];
                n += cell.pop_ave_n[j];
                vx += cell.pop_ave_n[j]*cell.pop_ave_vx[j];
                vy += cell.pop_ave_n[j]*cell.pop_ave_vy[j];
                vz += cell.pop_ave_n[j]*cell.pop_ave_vz[j];
            }
            real inv_n=0;
            if(n>0) {
                inv_n = 1.0/n;
            }
            vx *= inv_n;
            vy *= inv_n;
            vz *= inv_n;
            results[0] = vx;
            results[1] = vy;
            results[2] = vz;
            return results;
        }
    };
    struct nAveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            real n=0;
            for(unsigned int i=0; i< popId.size(); i++) {
                int j = popId[i];
                n += cell.pop_ave_n[j];
            }
            return vector<real> (1,n);
        }
    };
#endif
    struct rhoFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.calc_mass(popId) /
                                 (cell.size*cell.size*cell.size));
        }
    };
    struct TFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            real vx0=0,vy0=0,vz0=0;
            cell.plist.calc_avev(vx0,vy0,vz0,popId);

            return vector<real> (1, 2*cell.plist.calc_avemv2(vx0,vy0,vz0,popId)/(3*Params::k_B));
        }
    };
    struct EConvectiveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            const real Uex = cell.celldata[Tgrid::CELLDATA_UE][0];
            const real Uey = cell.celldata[Tgrid::CELLDATA_UE][1];
            const real Uez = cell.celldata[Tgrid::CELLDATA_UE][2];
            real Bx = cell.celldata[Tgrid::CELLDATA_B][0];
            real By = cell.celldata[Tgrid::CELLDATA_B][1];
            real Bz = cell.celldata[Tgrid::CELLDATA_B][2];
            if(Params::constantMagneticFieldProfile.size() > 0) {
                real B0[3] = {0.0, 0.0, 0.0};
                addConstantMagneticField(cell.centroid, B0);
                Bx += B0[0];
                By += B0[1];
                Bz += B0[2];
            }
            results.push_back(-(Uey*Bz-Uez*By));
            results.push_back(-(Uez*Bx-Uex*Bz));
            results.push_back(-(Uex*By-Uey*Bx));
            return results;
        }
    };
    struct EResistiveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            const real Jx = cell.celldata[Tgrid::CELLDATA_J][0];
            const real Jy = cell.celldata[Tgrid::CELLDATA_J][1];
            const real Jz = cell.celldata[Tgrid::CELLDATA_J][2];
            const real eta = Params::resistivityFunction.getValue(cell.centroid);
            results.push_back(eta*Jx);
            results.push_back(eta*Jy);
            results.push_back(eta*Jz);
            return results;
        }
    };
    struct EPolarizationFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            const real Epx = cell.celldata[Tgrid::CELLDATA_TEMP2][0];
            const real Epy = cell.celldata[Tgrid::CELLDATA_TEMP2][1];
            const real Epz = cell.celldata[Tgrid::CELLDATA_TEMP2][2];
            results.push_back(Epx);
            results.push_back(Epy);
            results.push_back(Epz);
            return results;
        }
    };
    struct UeFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(cell.celldata[Tgrid::CELLDATA_UE][comp]);
            return results;
        }
    };
    struct MacroParticlesFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.Nparticles());
        }
    };
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    struct SaveParticleFlagFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.save_particles);
        }
    };
#endif
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    // Spherical version of "B1Formula"
    struct sph_B1Formula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal B[3] = {cell.celldata[Tgrid::CELLDATA_B][0], cell.celldata[Tgrid::CELLDATA_B][1], cell.celldata[Tgrid::CELLDATA_B][2]};
            ////////////////////////////////////////////////////////////////////////////////////
            //gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            //gridreal B[3] = {0, 0, 0};
            //for (unsigned int comp = 0; comp < 3; ++comp)
            //    {
            //     B[comp] = 0.5*(cell.faceave(comp, 0, Tgrid::FACEDATA_B) + cell.faceave(comp, 1, Tgrid::FACEDATA_B));
            //    }
            ////////////////////////////////////////////////////////////////////////////////////
            sph_transf_S2C_A1(r, B);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(B[comp]);
            return results;
        }
    };
    // Spherical version of "B0Formula"
    struct sph_B0Formula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            real B0[3] = {0.0, 0.0, 0.0};
            addConstantMagneticField(cell.centroid, B0);

            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            sph_transf_S2C_A(r, B0);

            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(B0[comp]);
            return results;
        }
    };
    // Total magnetic field B0 + B1
    struct sph_B_totFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            real B0[3] = {0.0, 0.0, 0.0};
            addConstantMagneticField(cell.centroid, B0);
            real B[3] = {cell.celldata[Tgrid::CELLDATA_B][0], cell.celldata[Tgrid::CELLDATA_B][1], cell.celldata[Tgrid::CELLDATA_B][2]};
            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            int d;
            for (d = 0; d < 3; ++d) B[d] += B0[d];
            sph_transf_S2C_A(r, B);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(B[comp]);
            return results;
        }
    };
    // Average value of magnetic field
    struct sph_B1AveFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            ////////////////////////////////////////////////////////////////////////////////////
            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal B[3] = {0, 0, 0};
            for (unsigned int comp = 0; comp < 3; ++comp) {
                //B[comp] = 0.5*(cell.faceave(comp, 0, Tgrid::FACEDATA_AVEB) + cell.faceave(comp, 1, Tgrid::FACEDATA_AVEB));
                B[comp] = cell.faceave(comp, 1, Tgrid::FACEDATA_AVEB); // It gives the better result
            }
            ////////////////////////////////////////////////////////////////////////////////////
            sph_transf_S2C_A1(r, B);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(B[comp]);
            return results;
        }
    };
    // dB/dt = - rot(E) value
    struct sph_dBdt_Formula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            gridreal r[3]    = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal dBdt[3] = {cell.celldata[Tgrid::CELLDATA_B_TEMP1][0], cell.celldata[Tgrid::CELLDATA_B_TEMP1][1], cell.celldata[Tgrid::CELLDATA_B_TEMP1][2]};
            /*gridreal dBdt[3] = {0, 0, 0};
            for (unsigned int comp = 0; comp < 3; ++comp)
                {
                 dBdt[comp] = 0.5*(cell.faceave(comp, 0, Tgrid::FACEDATA_MINUSDB) + cell.faceave(comp, 1, Tgrid::FACEDATA_MINUSDB));
                }*/
            sph_transf_S2C_A1(r, dBdt);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(dBdt[comp]);
            return results;
        }
    };
    // E celldata
    struct sph_EFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal E[3] = {cell.celldata[Tgrid::CELLDATA_E][0], cell.celldata[Tgrid::CELLDATA_E][1], cell.celldata[Tgrid::CELLDATA_E][2]};
            sph_transf_S2C_A1(r, E);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(E[comp]);
            return results;
        }
    };
    // Spherical version of "vFormula"
    struct sph_vFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            vector<real> results(3);
            cell.plist.calc_avev(results[0], results[1], results[2], popId);
            return results;
        }
    };
    // Spherical version of "nFormula"
    struct sph_nFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.calc_weight(popId)/(cell.sph_dV));
            //return vector<real> (1, cell.nc);
        }
    };
    // Spherical version of "rhoFormula"
    struct sph_rhoFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.calc_mass(popId)/(cell.sph_dV));
        }
    };
    // rho_q calculation.
    struct sph_rhoqFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.plist.calc_charge(popId)/(cell.sph_dV));
            //return vector<real> (1, cell.rho_q);
        }
    };
    // rho_q_tot calculation.
    struct sph_rhoq_totFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            return vector<real> (1, cell.rho_q);
        }
    };
    // Spherical version of "EPolarizationFormula"
    struct sph_EPolarizationFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            gridreal r[3]  = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal Ep[3] = {cell.celldata[Tgrid::CELLDATA_TEMP2][0], cell.celldata[Tgrid::CELLDATA_TEMP2][1], cell.celldata[Tgrid::CELLDATA_TEMP2][2]};
            sph_transf_S2C_A1(r, Ep);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(Ep[comp]);
            return results;
        }
    };
    // Spherical version of "UeFormula"
    struct sph_UeFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            // CELLDATA_UE. Transformation from spherical to cartesian coordinates. To be represented in visit as a spherical vector
            gridreal r[3]  = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal Ue[3] = {cell.celldata[Tgrid::CELLDATA_UE][0], cell.celldata[Tgrid::CELLDATA_UE][1], cell.celldata[Tgrid::CELLDATA_UE][2]};
            sph_transf_S2C_A1(r, Ue);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(Ue[comp]);
            return results;
        }
    };
    // Spherical version of "JFormula"
    struct sph_JFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            // CELLDATA_J. Transformation from spherical to cartesian coordinates. To be represented in visit as a spherical vector
            gridreal r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal J[3] = {cell.celldata[Tgrid::CELLDATA_J][0], cell.celldata[Tgrid::CELLDATA_J][1], cell.celldata[Tgrid::CELLDATA_J][2]};
            sph_transf_S2C_A1(r, J);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(J[comp]);
            return results;
        }
    };
    // Spherical version of "JiFormula"
    struct sph_JiFormula {
        vector<real> operator()(const GridPriv::Tcell& cell,
                                const vector<int>& popId) {
            // If we calculate Ji in Cartesian coordinates we dont have to transform Ji.
            gridreal  r[3] = {cell.sph_centroid[0], cell.sph_centroid[1], cell.sph_centroid[2]};
            gridreal Ji[3] = {cell.celldata[Tgrid::CELLDATA_Ji][0], cell.celldata[Tgrid::CELLDATA_Ji][1], cell.celldata[Tgrid::CELLDATA_Ji][2]};
            //sph_transf_H2S_V(Ji);
            //sph_transf_S2C_V(r, Ji);
            vector<real> results;
            for (unsigned int comp = 0; comp < 3; ++comp)
                results.push_back(Ji[comp]);
            return results;
        }
    };
#endif
};

/** \brief Struct that contains formulas for particle variables.
 */
struct ParticleFormulas {
    /** \brief Retuns population ID of given particle.
     */
    struct popIDFormula {
        vector<real> operator()(const TLinkedParticle& particle,
                                const vector<int>& popId) {
            return vector<real> (1, particle.popid);
        }
    };
};

/** \brief Struct that defines what kind of variable are written to file.
 */
struct VariableDefinitions {
    /** \brief Struct that defines common attributes of vector variables.
     */
    struct VectorVariableCommon {
        //! Name of the variable
        string name;
        //! Number of components in variable
        unsigned int nComponents;
        //! Whether there is also an expression defining this variable
        bool hasExpression;
        /** \brief Creates new VectorVariable with given attributes.
         */
        VectorVariableCommon(const string& name, unsigned int nComponents,
                             bool hasExpression)
            : name(name), nComponents(nComponents), hasExpression(hasExpression) { }
    };
    /** \brief Struct that defines attributes of cell vector variables.
     */
    struct CellVectorVariable : public VectorVariableCommon {
        //! Whether variable is population variable
        bool isPopulationVariable;
        //! Formula that returns value of the variable in given cell
        Function<vector<real> (const GridPriv::Tcell& cell,
                               const vector<int>& popId)> formula;
        /** \brief Creates new CellVectorVariable with given attributes.
         */
        CellVectorVariable
        (const string& name, unsigned int nComponents, bool hasExpression,
         bool isPopulationVariable,
         const Function<vector<real> (const GridPriv::Tcell& cell,
                                      const vector<int>& popId)>& formula)
            : VectorVariableCommon(name, nComponents, hasExpression),
              isPopulationVariable(isPopulationVariable), formula(formula) { }
    };
    /** \brief Struct that defines common attributes of particle
     *         vector variables.
     */
    struct ParticleVectorVariable : public VectorVariableCommon {
        //! Formula that returns value of the variable in given cell
        Function<vector<real> (const TLinkedParticle& particle,
                               const vector<int>& popId)> formula;
        /** \brief Creates new ParticleVectorVariable with given attributes.
         */
        ParticleVectorVariable
        (const string& name, unsigned int nComponents, bool hasExpression,
         const Function<vector<real> (const TLinkedParticle& particle,
                                      const vector<int>& popId)>& formula)
            : VectorVariableCommon(name, nComponents, hasExpression),
              formula(formula) { }
    };
    //! Cell vector variables that are written to file
    vector<CellVectorVariable> cellVectorVariables;
    //! Particle vector variables that are written to file
    vector<ParticleVectorVariable> particleVectorVariables;
};

/** \brief Struct that defines what kinds of expressions are written to file.
 */
struct ExpressionDefinitions {
    struct Expression {
        //! Name of the expression
        string name;
        //! Whether expression is population expression
        bool isPopulationExpression;
        //! Formula of expression
        string formula;
        /** \brief Creates new Expression with given attributes.
         */
        Expression(const string& name, bool isPopulationExpression,
                   const string& formula)
            : name(name), isPopulationExpression(isPopulationExpression),
              formula(formula) { }
    };
    vector<Expression> expressions;
};

/** \brief Class that gives a "view" to array of grid's Tcells. Contains
 *         iterators for traversing cells and slicing operations.
 */
class GridCells
{
public:
    /** \brief Creates new GridCells from given cells.
     * \param cells Array of Tcells to which GridCells gives a view.
     * \param dimensions dimensions (size) of given array of Tcells
     */
    GridCells(GridPriv::Tcell** cells, const Coordinates<3>& dimensions)
        : m_cells(cells), m_dimensions(dimensions),
          m_slice(Coordinates<3>(0, 0, 0), Coordinates<3>(dimensions)) {
        for (int i = 0; i <= Params::maxGridRefinementLevel; ++i)
            m_refLevels.push_back(i);
    }
    /** \brief Iterator that iterates through all cells, even those
     *         that are empty in current refinement level.
     * Values that iterator points to are pointers to Tcells, not
     * references. Points to null, if current cell is empty on
     * current refinement level.
     */
    class allIterator
        : public std::iterator<std::forward_iterator_tag, GridPriv::Tcell*>
    {
    public:
        allIterator(GridCells& cells, const Coordinates<3>& startCell,
                    unsigned int startRefLevelIdx)
            : m_cells(&cells), m_currentCell(startCell),
              m_currentRefLevelIdx(startRefLevelIdx),
              m_currentInterval(cells.intervalAtRefLevelIdx(startRefLevelIdx)),
              m_endCoordinates(cells.intervalAtRefLevelIdx
                               (cells.m_refLevels.size() - 1).getEnd() +
                               Coordinates<3>(1, 0, 0)) { }
        allIterator(GridCells& cells, IteratorPosition startPosition)
            : m_cells(&cells),
              m_endCoordinates(cells.intervalAtRefLevelIdx
                               (cells.m_refLevels.size() - 1).getEnd() +
                               Coordinates<3>(1, 0, 0)) {
            switch (startPosition) {
            case BEGIN: {
                m_currentCell = cells.m_slice.getBegin();
                m_currentRefLevelIdx = 0;
                m_currentInterval = cells.m_slice;
                break;
            }
            case END: {
                m_currentCell = m_endCoordinates;
                m_currentRefLevelIdx = cells.m_refLevels.size() - 1;
                m_currentInterval = cells.intervalAtRefLevelIdx(m_currentRefLevelIdx);
                break;
            }
            }
        }
        allIterator& operator++() {
            assert (m_currentCell != m_endCoordinates);
            if (!incrementCoordinatesInInterval(m_currentCell, m_currentInterval)) {
                if (m_currentRefLevelIdx == m_cells->m_refLevels.size() - 1)
                    m_currentCell = m_endCoordinates;
                else {
                    ++m_currentRefLevelIdx;
                    unsigned int refRatio =
                        twoPower(m_cells->m_refLevels[m_currentRefLevelIdx]);
                    m_currentCell = m_currentInterval.getBegin() * refRatio;
                    m_currentInterval = m_cells->m_slice * refRatio;
                }
            }
            return *this;
        }
        allIterator operator++(int) {
            allIterator prev = *this;
            ++(*this);
            return prev;
        }
        bool operator==(const allIterator& other) const {
            return m_currentCell == other.m_currentCell;
        }
        bool operator!=(const allIterator& other) const {
            return !(*this == other);
        }
        GridPriv::Tcell* operator*() const {
            return m_cells->getCell(m_currentCell, m_currentRefLevelIdx);
        }
        GridPriv::Tcell* operator->() const {
            return m_cells->getCell(m_currentCell, m_currentRefLevelIdx);
        }
        /** \brief Returns coordinates of current cell in
         *         current refinement level.
         */
        Coordinates<3> getCoordinates() const {
            return m_currentCell - m_currentInterval.getBegin();
        }
        /** \brief Returns refinement level of current cell.
         */
        unsigned int getRefLevel() const {
            return m_cells->m_refLevels[m_currentRefLevelIdx];
        }
    private:
        GridCells* m_cells;
        Coordinates<3> m_currentCell;
        unsigned int m_currentRefLevelIdx;
        Interval<3> m_currentInterval;
        const Coordinates<3> m_endCoordinates;
    };
    /** \brief const allIterator
     */
    class const_allIterator
        : public std::iterator<std::input_iterator_tag, const GridPriv::Tcell*>
    {
    public:
        const_allIterator(const GridCells& cells, const Coordinates<3>& startCell,
                          unsigned int startRefLevelIdx)
            : m_cells(&cells), m_currentCell(startCell),
              m_currentRefLevelIdx(startRefLevelIdx),
              m_currentInterval(cells.intervalAtRefLevelIdx(startRefLevelIdx)),
              m_endCoordinates(cells.intervalAtRefLevelIdx
                               (cells.m_refLevels.size() - 1).getEnd() +
                               Coordinates<3>(1, 0, 0)) { }
        const_allIterator(const GridCells& cells, IteratorPosition startPosition)
            : m_cells(&cells),
              m_endCoordinates(cells.intervalAtRefLevelIdx
                               (cells.m_refLevels.size() - 1).getEnd() +
                               Coordinates<3>(1, 0, 0)) {
            switch (startPosition) {
            case BEGIN: {
                m_currentCell = cells.m_slice.getBegin();
                m_currentRefLevelIdx = 0;
                m_currentInterval = cells.m_slice;
                break;
            }
            case END: {
                m_currentCell = m_endCoordinates;
                m_currentRefLevelIdx = cells.m_refLevels.size() - 1;
                m_currentInterval = cells.intervalAtRefLevelIdx(m_currentRefLevelIdx);
                break;
            }
            }
        }
        const_allIterator& operator++() {
            assert (m_currentCell != m_endCoordinates);
            if (!incrementCoordinatesInInterval(m_currentCell, m_currentInterval)) {
                if (m_currentRefLevelIdx == m_cells->m_refLevels.size() - 1)
                    m_currentCell = m_endCoordinates;
                else {
                    ++m_currentRefLevelIdx;
                    unsigned int refRatio =
                        twoPower(m_cells->m_refLevels[m_currentRefLevelIdx]);
                    m_currentCell = m_currentInterval.getBegin() * refRatio;
                    m_currentInterval = m_cells->m_slice * refRatio;
                }
            }
            return *this;
        }
        const_allIterator operator++(int) {
            const_allIterator prev = *this;
            ++(*this);
            return prev;
        }
        bool operator==(const const_allIterator& other) const {
            return m_currentCell == other.m_currentCell;
        }
        bool operator!=(const const_allIterator& other) const {
            return !(*this == other);
        }
        const GridPriv::Tcell* operator*() const {
            return m_cells->getCell(m_currentCell, m_currentRefLevelIdx);
        }
        const GridPriv::Tcell* operator->() const {
            return m_cells->getCell(m_currentCell, m_currentRefLevelIdx);
        }
        /** \brief Returns coordinates of current cell in
         *         current refinement level.
         */
        Coordinates<3> getCoordinates() const {
            return m_currentCell - m_currentInterval.getBegin();
        }
        /** \brief Returns refinement level of current cell.
         */
        unsigned int getRefLevel() const {
            return m_cells->m_refLevels[m_currentRefLevelIdx];
        }
    private:
        const GridCells* m_cells;
        Coordinates<3> m_currentCell;
        unsigned int m_currentRefLevelIdx;
        Interval<3> m_currentInterval;
        Coordinates<3> m_endCoordinates;
    };
    /** \brief Iterator that iterates through only those cells that
     *	       are not empty on current refinement level.
     */
    class iterator
        : public std::iterator<std::forward_iterator_tag, GridPriv::Tcell>
    {
    public:
        iterator(GridCells& cells, const Coordinates<3>& startCell,
                 unsigned int startRefLevelIdx)
            : m_currentPosition(GridCells::allIterator
                                (cells, startCell, startRefLevelIdx)),
            m_endPosition(GridCells::allIterator(cells, END)) { }
        iterator(GridCells& cells, IteratorPosition startPosition)
            : m_currentPosition
            (GridCells::allIterator(cells, startPosition)),
            m_endPosition(GridCells::allIterator(cells, END)) { }
        iterator& operator++() {
            assert(m_currentPosition != m_endPosition);
            do {
                ++m_currentPosition;
            } while (m_currentPosition != m_endPosition && *m_currentPosition == 0);
            return *this;
        }
        iterator operator++(int) {
            iterator prev = *this;
            ++(*this);
            return prev;
        }
        bool operator==(const iterator& other) const {
            return m_currentPosition == other.m_currentPosition;
        }
        bool operator!=(const iterator& other) const {
            return !(*this == other);
        }
        GridPriv::Tcell& operator*() const {
            return *(*m_currentPosition);
        }
        GridPriv::Tcell* operator->() const {
            return *m_currentPosition;
        }
        /** \brief Returns coordinates of current cell in
         *         current refinement level.
         */
        Coordinates<3> getCoordinates() const {
            return m_currentPosition.getCoordinates();
        }
        /** \brief Returns refinement level of current cell.
         */
        unsigned int getRefLevel() const {
            return m_currentPosition.getRefLevel();
        }
    private:
        GridCells::allIterator m_currentPosition;
        GridCells::allIterator m_endPosition;
    };
    /** \brief Const iterator
     */
    class const_iterator
        : public std::iterator<std::input_iterator_tag, const GridPriv::Tcell>
    {
    public:
        const_iterator(const GridCells& cells, const Coordinates<3>& startCell,
                       unsigned int startRefLevelIdx)
            : m_currentPosition(GridCells::const_allIterator
                                (cells, startCell, startRefLevelIdx)),
            m_endPosition(GridCells::const_allIterator(cells, END)) { }
        const_iterator(const GridCells& cells, IteratorPosition startPosition)
            : m_currentPosition
            (GridCells::const_allIterator(cells, startPosition)),
            m_endPosition(GridCells::const_allIterator(cells, END)) { }
        const_iterator& operator++() {
            assert (m_currentPosition != m_endPosition);
            do {
                ++m_currentPosition;
            } while (m_currentPosition != m_endPosition && *m_currentPosition == 0);
            return *this;
        }
        const_iterator operator++(int) {
            const_iterator prev = *this;
            ++(*this);
            return prev;
        }
        bool operator==(const const_iterator& other) const {
            return m_currentPosition == other.m_currentPosition;
        }
        bool operator!=(const const_iterator& other) const {
            return !(*this == other);
        }
        const GridPriv::Tcell& operator*() const {
            return *(*m_currentPosition);
        }
        const GridPriv::Tcell* operator->() const {
            return *m_currentPosition;
        }
        /** \brief Returns coordinates of current cell in
         *         current refinement level.
         */
        Coordinates<3> getCoordinates() const {
            return m_currentPosition.getCoordinates();
        }
        /** \brief Returns refinement level of current cell.
         */
        unsigned int getRefLevel() const {
            return m_currentPosition.getRefLevel();
        }
    private:
        GridCells::const_allIterator m_currentPosition;
        GridCells::const_allIterator m_endPosition;
    };
    /** \brief Returns GridCells that is a subset of this GridCells.
     *
     * All refinement levels are in new subset.
     *
     * \param interval Which cells are in new subset. Coordinates of interval
     *                 are coordinates of highest refinement level of this
     *                 GridCells.
     */
    GridCells slice(const Interval<3>& interval) const {
        return slice(interval, m_refLevels);
    }
    /** \brief Returns GridCells that is a subset of this GridCells.
     *
     * \param interval Which cells are in new subset. Coordinates of interval
     *                 are coordinates of highest refinement level of this
     *                 GridCells.
     * \param refLevel Which refinement level is in new subset.
     */
    GridCells slice(const Interval<3>& interval,
                    unsigned int refLevel) const {
        return slice(interval, vector<unsigned int>(1, refLevel));
    }
    /** \brief Returns GridCells that is a subset of this GridCells.
     *
     * \param interval Which cells are in new subset. Coordinates of interval
     *                 are coordinates of highest refinement level of this
     *                 GridCells.
     * \param refLevels Which refinement levels are in new subset.
     */
    GridCells slice(const Interval<3>& interval,
                    const vector<unsigned int>& refLevels) const {
        GridCells newCells(m_cells, m_dimensions);
        unsigned int factor = twoPower(refLevels[0] - m_refLevels[0]);
        Coordinates<3> currentPos = m_slice.getBegin() * factor;
        Interval<3> newSlice = Interval<3>(interval.getBegin() + currentPos,
                                           interval.getEnd() + currentPos);
        assert(isSubSet(newSlice, m_slice * factor));
        newCells.m_slice =
            Interval<3>(interval.getBegin() + currentPos,
                        interval.getEnd() + currentPos);
        newCells.m_refLevels = refLevels;
        return newCells;
    }
    /** \brief Returns Interval that shows which cells are in current
     *         subset of all cells.
     *
     * Coordinates of Interval are coordinates of highest refinement
     * level of current subset.
     */
    const Interval<3>& getSliceInterval() const {
        return m_slice;
    }
    /** \brief Returns dimensions (size) of current subset.
     *
     * Coordinates of size are coordinates in highest refinement
     * level of current subset.
     */
    Coordinates<3> getDimensions() const {
        return m_slice.getEnd() - m_slice.getBegin();
    }
    /** \brief Returns refinement levels that are in current subset
     *         in ascending order.
     */
    vector<unsigned int> getRefinementLevels() const {
        return m_refLevels;
    }
    /** \brief Returns iterator to begin of cells.
     */
    iterator begin() {
        return iterator(*this, BEGIN);
    }
    /** \brief Returns const_iterator to begin of cells.
     */
    const_iterator begin() const {
        return const_iterator(*this, BEGIN);
    }
    /** \brief Returns iterator to end of cells.
     */
    iterator end() {
        return iterator(*this, END);
    }
    /** \brief Returns const_iterator to end of cells.
     */
    const_iterator end() const {
        return const_iterator(*this, END);
    }
    /** \brief Returns allIterator to begin of cells.
     */
    allIterator allBegin() {
        return allIterator(*this, BEGIN);
    }
    /** \brief Returns const_allIterator to begin of cells.
     */
    const_allIterator allBegin() const {
        return const_allIterator(*this, BEGIN);
    }
    /** \brief Returns allIterator to end of cells.
     */
    allIterator allEnd() {
        return allIterator(*this, END);
    }
    /** \brief Returns const_allIterator to end of cells.
     */
    const_allIterator allEnd() const {
        return const_allIterator(*this, END);
    }
private:
    GridPriv::Tcell** m_cells;
    const Coordinates<3> m_dimensions;
    Interval<3> m_slice;                    //! Current subset of cells
    std::vector<unsigned int> m_refLevels;  //! Ref. levels in current subset
    std::size_t flatIndex(Coordinates<3> coords) const {
        return (coords[0] * m_dimensions[1] + coords[1]) * m_dimensions[2] +
               coords[2];
    }
    Interval<3> intervalAtRefLevelIdx(unsigned int idx) const {
        return m_slice * twoPower(m_refLevels[idx] - m_refLevels[0]);
    }
    /** \brief Returns pointer to cell at given coordinates and refinement
     *         level, or null, if cell is empty at given refinement level.
     */
    GridPriv::Tcell* getCell(const Coordinates<3>& coords,
                             unsigned int refLevelIdx) const;
};

GridPriv::Tcell* GridCells::getCell(const Coordinates<3>& coords,
                                    unsigned int refLevelIdx) const
{
    const unsigned int refLevel = m_refLevels[refLevelIdx];
    // Optimization for case refLevel == 0
    if (refLevel == 0) {
        GridPriv::Tcell* cell = m_cells[flatIndex(coords)];
        if (!cell->haschildren)
            return cell;
        else
            return 0;
    }
    unsigned int refRatio = twoPower(refLevel);
    Coordinates<3> ancestorCellCoords = coords / refRatio;
    Coordinates<3> restCoords = coords - ancestorCellCoords * refRatio;
    GridPriv::Tcell* cell = m_cells[flatIndex(ancestorCellCoords)];
    while (refRatio > 1) {
        if (!cell->haschildren)
            return 0;
        refRatio /= 2;
        ancestorCellCoords = restCoords / refRatio;
        cell = cell->child[ancestorCellCoords[0]][ancestorCellCoords[1]]
               [ancestorCellCoords[2]];
        restCoords -= ancestorCellCoords * refRatio;
    }
    if (!cell->haschildren)
        return cell;
    else
        return 0;
}

/** \brief Returns true if given cell is empty, otherwise returns false.
 */
struct EmptyCell {
    bool operator()(const GridPriv::Tcell* cell) {
        return cell == 0;
    }
};

Coordinates<3> findPatch(const GridCells& cells,
                         const Coordinates<3>& begin,
                         unsigned int refLevel,
                         const vector<vector<Interval<3> > >& readyPatches);
/** \brief Finds and returns AMR patches (regions with same refinement level)
 *         of given cells.
 */
vector<vector<Interval<3> > > findAMRPatches(const GridCells& cells)
{
    vector<vector<Interval<3> > > patches
        = vector<vector<Interval<3> > >
          (cells.getRefinementLevels()[cells.getRefinementLevels().size() - 1] + 1);
    // If cells have only one refinement level, every cell can be in same patch
    if (cells.getRefinementLevels().size() == 1) {
        patches[cells.getRefinementLevels()[0]].push_back
        (Interval<3>(Coordinates<3>(0, 0, 0), cells.getDimensions()));
    } else {
        for (GridCells::const_iterator cell = cells.begin();
             cell != cells.end(); ++cell) {
            if (!coordinatesInInterval(cell.getCoordinates(),
                                       patches[cell.getRefLevel()]))
                patches[cell.getRefLevel()].push_back
                (Interval<3>(cell.getCoordinates(),
                             findPatch(cells, cell.getCoordinates(),
                                       cell.getRefLevel(), patches)));
        }
    }
    return patches;
}

/** \brief Finds and returns a new AMR patch (see function findAMRPatches)
 *         that begin at given coordinates.
 *
 * Searches only for patches that are at given refinement level and
 * not in patches given in argument "readyPatches". Assumes that
 * given begin coordinates are "free", so at least this cell
 * belongs in new patch.
 */
Coordinates<3> findPatch(const GridCells& cells,
                         const Coordinates<3>& begin,
                         unsigned int refLevel,
                         const vector<vector<Interval<3> > >& readyPatches)
{
    std::size_t dist = 1;  // current size of patch (length of edge of cube)
    while(true) {
        bool cont = true;
        // Stop if coordinates are outside of cells
        if (!coordinatesInInterval(begin + Coordinates<3>(dist, dist, dist),
                                   Interval<3>(Coordinates<3>(0, 0, 0),
                                               cells.getDimensions() *
                                               twoPower(refLevel))))
            break;
        for (Dimension dim = 0; dim < 3; ++dim) {
            // Take new surface that is tried to be added to patch
            Interval<3> surface = cubeSurface(begin, dist, CoordinateID(dim));
            GridCells slice = cells.slice(surface, refLevel);
            /* Check that new side (surface) of cube doesn't intersect
             * any already found patch or other refinement level. */
            if (std::find_if(readyPatches[refLevel].begin(),
                             readyPatches[refLevel].end(),
                             IntervalsIntersect<3>(surface)) !=
                readyPatches[refLevel].end() ||
                std::find_if(slice.allBegin(), slice.allEnd(), EmptyCell()) !=
                slice.allEnd()) {
                cont = false;
                break;
            }
        }
        if (!cont)
            break;
        ++dist;
    }
    return Coordinates<3>(begin[0] + dist, begin[1] + dist, begin[2] + dist);
}

/** \brief Calls given function and puts result values (that it assumes to
 *         be in vector) to corresponding indexes in vector "components".
 *
 * Passes argument "popId" as second argument to function that is called.
 * At first call, puts result values to first indexes in arrays in
 * vector "components". Increases this index by one at every call.
 */
template <class T, class Func> struct ArrayWriter {
    ArrayWriter(Func f, const vector<int>& popId, vector<T*>& components)
        : m_function(f), m_popId(popId), m_components(components),
          m_idx(0), m_retVals(components.size()) { }
    template <class Arg>
    bool operator()(Arg& argument) {
        m_retVals = m_function(argument, m_popId);
        for (unsigned int i = 0; i < m_components.size(); ++i)
            m_components[i][m_idx] = m_retVals[i];
        ++m_idx;
        return true;
    }
private:
    Func m_function;
    const vector<int>& m_popId;
    vector<T*>& m_components;
    std::size_t m_idx;
    vector<T> m_retVals;
};

/** \brief Convenience function that returns new ArrayWriter with given
 *         constructor parameters.
 *
 * With this function, template parameters of ArrayWriter don't need to
 * be written, they are deduced automatically.
 */
template <class T, class Func>
ArrayWriter<T, Func> getArrayWriter(Func function, const vector<int>& popId,
                                    vector<T*>& components)
{
    return ArrayWriter<T, Func>(function, popId, components);
}

/** \brief Returns Particle that contains data from given particle.
 */
struct ParticleDataGetter {
    CommonDataSourceData::Particle operator()
    (const TLinkedParticle& particle) const {
        CommonDataSourceData::Particle p = {
            SequenceHandle<real>
            (new ArraySequence<real>(particle.x, particle.y, particle.z))
        };
        return p;
    }
};

/** \brief Returns cell vector variable that contains data from variable
 *         at given index "varIndex".
 */
struct CellVectorVariableGenerator {
    struct CellVariableData;
    CellVectorVariableGenerator(const CellVariableData& varData)
        : m_varData(varData) { }
    CommonDataSourceData::VectorVariable operator()
    (unsigned int varIndex) const {
        typedef CommonDataSourceData::VectorVariable VectorVar;
        CellVariableData::VectorVariable var =
            m_varData.vectorVariables[varIndex];
        vector<vector<real*> > varValues(var.nComponents);
        for (unsigned int refLevel = 0; refLevel < m_varData.patches.size();
             ++refLevel)
            for (unsigned int patchIdx = 0;
                 patchIdx < m_varData.patches[refLevel].size(); ++patchIdx) {
                vector<real*> patchComponents;
                for (Dimension compIdx = 0;
                     compIdx < var.nComponents; ++compIdx) {
                    real* componentValues =
                        new real[volume(m_varData.patches[refLevel][patchIdx])];
                    varValues[compIdx].push_back(componentValues);
                    patchComponents.push_back(componentValues);
                }
                GridCells patchCells = m_varData.cells.slice
                                       (m_varData.patches[refLevel][patchIdx], refLevel);
                std::for_each(patchCells.begin(), patchCells.end(),
                              getArrayWriter(var.formula,
                                             m_varData.popId[var.popIdIdx],
                                             patchComponents));
            }
        VectorVar::Component* components =
            new VectorVar::Component[var.nComponents];
        for (Dimension compIdx = 0; compIdx < var.nComponents; ++compIdx) {
            components[compIdx].name = var.name + static_cast<char>('x' + compIdx);
            components[compIdx].values = SequenceHandle<real*>
                                         (new AutoADelArraySequence<real*>(varValues[compIdx]));
        }
        VectorVar v;
        v.name = var.name;
        v.components = SequenceHandle<VectorVar::Component>
                       (new ArraySequence<VectorVar::Component>(components, var.nComponents));
        return v;
    }

    /** \brief Struct containing data from which data of variables
     *         is determined.
     */
    struct CellVariableData {
        struct VectorVariable;
        GridCells cells;
        vector<vector<Interval<3> > > patches;
        vector<vector<int> > popId;
        vector<VectorVariable> vectorVariables;
        struct VectorVariable {
            string name;
            unsigned int nComponents;
            Function<vector<real> (const GridPriv::Tcell& cell,
                                   const vector<int>& popId)> formula;
            int popIdIdx;
            VectorVariable(const string& name, unsigned int nComponents,
                           const Function<vector<real>
                           (const GridPriv::Tcell& cell,
                            const vector<int>& popId)>& formula,
                           int popIdIdx)
                : name(name), nComponents(nComponents), formula(formula),
                  popIdIdx(popIdIdx) { }
        };
        CellVariableData (const GridCells& cells,
                          const vector<vector<Interval<3> > >& patches,
                          const vector<vector<int> >& popId,
                          const vector<VectorVariable>& vectorVariables)
            : cells(cells), patches(patches), popId(popId),
              vectorVariables(vectorVariables) { }
    };

private:
    const CellVariableData m_varData;
};

/** \brief Returns particle vector variable that contains data from variable
 *         at given index "varIndex".
 */
struct ParticleVectorVariableGenerator {
    struct ParticleVariableData;
    ParticleVectorVariableGenerator(const ParticleVariableData& varData)
        : m_varData(varData) { }
    CommonDataSourceData::VectorVariable operator()(unsigned int varIndex) const {
        typedef CommonDataSourceData::VectorVariable VectorVar;
        ParticleVariableData::VectorVariable var =
            m_varData.vectorVariables[varIndex];
        VectorVar v;
        v.name = var.name;
        v.components = SequenceHandle<VectorVar::Component>
                       (new ArraySequence<VectorVar::Component>
                        (new VectorVar::Component[var.nComponents], var.nComponents));
        vector<real*> compVals;
        for (Dimension compIdx = 0; compIdx < var.nComponents; ++compIdx) {
            v.components[compIdx].name = var.name +
                                         static_cast<char>('x' + compIdx);
            v.components[compIdx].values = SequenceHandle<real*>
                                           (new AutoADelArraySequence<real*>(new real*[1], 1));
            v.components[compIdx].values[0] = new real[m_varData.nParticles];
            compVals.push_back(v.components[compIdx].values[0]);
        }
        ArrayWriter<real,  Function<vector<real>
        (const TLinkedParticle& particle,
         const vector<int>& popId)> > valueWriter =
             getArrayWriter(var.formula, vector<int>(0), compVals);
        for (GridCells::const_iterator cell = m_varData.cells.begin();
             cell != m_varData.cells.end(); ++cell)
            cell->plist.pass(valueWriter);
        return v;
    }

    /** \brief Struct containing data from which data of variables
     *         is determined.
     */
    struct ParticleVariableData {
        struct VectorVariable;
        GridCells cells;
        std::size_t nParticles;
        vector<VectorVariable> vectorVariables;
        struct VectorVariable {
            string name;
            unsigned int nComponents;
            Function<vector<real> (const TLinkedParticle& particle,
                                   const vector<int>& popId)> formula;
            VectorVariable(const string& name, unsigned int nComponents,
                           const Function<vector<real>
                           (const TLinkedParticle& particle,
                            const vector<int>& popId)>& formula)
                : name(name), nComponents(nComponents), formula(formula) { }
        };
        ParticleVariableData(const GridCells& cells, std::size_t nParticles,
                             const vector<VectorVariable>& vectorVariables)
            : cells(cells), nParticles(nParticles), vectorVariables(vectorVariables) {
        }
    };
private:
    ParticleVariableData m_varData;
};

/** \brief Currently unimplemented function that returns data of
 *         cell at given index "cellIndex".
 */
struct VisCellGenerator {
    VisData::Patch::Cell operator()(unsigned int cellIndex) const {
        VisData::Patch::Cell c;
        return c;
    }
};

/** \brief Converts vector variables from VariableDefinitions
 *         to suitable form (CellVariableData::VectorVariable)
 *         for CellVectorVariableGenerator.
 */
vector<CellVectorVariableGenerator::CellVariableData::VectorVariable>
getCellVectorVariableData(const VariableDefinitions& defs)
{
    typedef CellVectorVariableGenerator::CellVariableData::VectorVariable
    CellVar;
    vector<string> hcFilePrefix;
    vector< vector<int> > popId;
    Population::getHcFileConfigs(hcFilePrefix, popId);
    vector<CellVar> vars;
    for (vector<VariableDefinitions::CellVectorVariable>::const_iterator var
         = defs.cellVectorVariables.begin();
         var != defs.cellVectorVariables.end(); ++var)
        if (!var->isPopulationVariable)
            vars.push_back(CellVar(var->name, var->nComponents, var->formula, 0));
        else
            for (unsigned int popu = 0; popu < hcFilePrefix.size(); ++popu) {
                vars.push_back(CellVar(var->name + "_" + hcFilePrefix[popu],
                                       var->nComponents, var->formula, popu));
            }
    return vars;
}

/** \brief Converts vector variables from VariableDefinitions
 *         to suitable form (ParticleVariableData::VectorVariable)
 *         for CellVectorVariableGenerator.
 */
vector<ParticleVectorVariableGenerator::ParticleVariableData::VectorVariable>
getParticleVectorVariableData(const VariableDefinitions& defs)
{
    typedef ParticleVectorVariableGenerator::ParticleVariableData::
    VectorVariable ParticleVar;
    vector<ParticleVar> vars;
    for (vector<VariableDefinitions::ParticleVectorVariable>::const_iterator
         var = defs.particleVectorVariables.begin();
         var != defs.particleVectorVariables.end(); ++var)
        vars.push_back(ParticleVar(var->name, var->nComponents, var->formula));
    return vars;
}

}

/** \brief When operator() is called, returns result when function
 *         given to constructor is applied to particle at index
 *         "particleIdx".
 *
 * Particles are indexed at arbitrary order.
 * Caches location of particle when operator() is called, so
 * it's fast to get particle at following index.
 * Template parameter Ret is return type of the function and Func
 * is type of the function.
 */
template <class Ret, class Func>
struct ParticleMapper {
    struct ParticleData;
    /** \brief Creates new ParticleMapper.
     *
     * \param pData Determines which particles are "iterated" through.
     * \param f Function that is applied to particle.
     */
    ParticleMapper(const ParticleData& pData,
                   Func f)
        : m_pData(new ParticleData(pData)), m_func(f), prevIdx(0),
          currentCell(m_pData->cells.begin()) {
        findFirstParticleFromCurrentCell();
    }
    Ret operator()(std::size_t particleIdx) const {
        if (particleIdx < prevIdx) {
            currentCell = m_pData->cells.begin();
            findFirstParticleFromCurrentCell();
            prevIdx = 0;
        }
        while(prevIdx != particleIdx)
            gotoNextParticle();
        return m_func(*prevParticle);
    }
    struct ParticleData {
        const GridCells cells;
    };
private:
    SharedPtr<ParticleData> m_pData;
    Func m_func;
    mutable std::size_t prevIdx;
    mutable GridCells::const_iterator currentCell;
    mutable const TLinkedParticle* prevParticle;
    void gotoNextParticle() const {
        prevParticle = prevParticle->next;
        while (prevParticle == 0) {
            ++currentCell;
            prevParticle = currentCell->plist.first;
        }
        ++ prevIdx;
    }
    void findFirstParticleFromCurrentCell() const {
        prevParticle = currentCell->plist.first;
        while (prevParticle == 0) {
            ++currentCell;
            prevParticle = currentCell->plist.first;
        }
    }
};

//! "Pimpl" implementation class of SimulationVisDataSourceImpl
class SimulationVisDataSourceImplPrivate
{
public:
    SimulationVisDataSourceImplPrivate(const Params& p, const Tgrid& g);
    //! Params instance from which simulation data is taken
    const Params& m_params;
    //! Grid instance from which simulation data is taken
    const Tgrid& m_grid;
    GridCells m_cells;                        //! Cells whose data is visualized
    VariableDefinitions m_varDefinitions;     //! Variables that are visualized
    ExpressionDefinitions m_exprDefinitions;  //! Expressions that are visualized
    vector<vector<Interval<3> > > m_patches;  //! AMR patches in grid
    VisData getVisData();
    /** \brief Defines what kind of variables and expressions are visualized
     */
    void defineData();
    //! \brief Returns patches in suitable form for VisData
    ConstSequenceHandle<VisData::Patch> getPatches();
    //! \brief Returns particles in suitable form for VisData
    ConstSequenceHandle<VisData::Particle> getParticles();
    //! \brief Returns cell vector variables in suitable form for VisData
    ConstSequenceHandle<VisData::VectorVariable> getCellVectorVariables();
    //! \brief Returns particle vector variables in suitable form for VisData
    ConstSequenceHandle<VisData::VectorVariable> getParticleVectorVariables();
    //! \brief Returns expressions in suitable form for VisData
    ConstSequenceHandle<VisData::Expression> getExpressions();
};

SimulationVisDataSourceImplPrivate::SimulationVisDataSourceImplPrivate
(const Params& p, const Tgrid& g)
    : m_params(p), m_grid(g),
      m_cells(GridCells(g.cells, Coordinates<3>
#ifndef VTK_SHOW_GHOST_CELLS
                        (Params::nx+2, Params::ny+2, Params::nz+2)).slice(Interval<3>(Coordinates<3>(1, 1, 1), Coordinates<3>(Params::nx+1, Params::ny+1, Params::nz+1))))
#else
                        (Params::nx+2, Params::ny+2, Params::nz+2)))
#endif
{
    m_patches = findAMRPatches(m_cells);
    defineData();
}

VisData SimulationVisDataSourceImplPrivate::getVisData()
{
    VisData data;
    data.time = Params::t;
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
#ifndef VTK_SHOW_GHOST_CELLS
    data.startCoordinates = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.x_1 + m_grid.bgdx, m_grid.y_1 + m_grid.bgdx,m_grid.z_1 + m_grid.bgdx));
#else
    data.startCoordinates = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.x_1, m_grid.y_1, m_grid.z_1));
#endif
    data.refRatio1CellSize = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.bgdx, m_grid.bgdx, m_grid.bgdx));
#else // spherical
#ifndef VTK_SHOW_GHOST_CELLS
    // 1. don't include ghost cells
    data.startCoordinates = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.x_1 + m_grid.bgdx, m_grid.y_1 + m_grid.sph_bgdy, m_grid.z_1 + m_grid.sph_bgdz));
    //! 1a. don't include ghost cells. First and last layers in x(r) direction are excluded
    //(new ArraySequence<real>(m_grid.x_1 + 2*m_grid.bgdx, m_grid.y_1 + m_grid.sph_bgdy, m_grid.z_1 + m_grid.sph_bgdz));
#else
    //! 2. include ghost cells
    data.startCoordinates = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.x_1, m_grid.y_1, m_grid.z_1));
#endif
    data.refRatio1CellSize = ConstSequenceHandle<real>(new ArraySequence<real>(m_grid.bgdx, m_grid.sph_bgdy, m_grid.sph_bgdz));
#endif
    data.patches = getPatches();
    //data.particles = getParticles();
    data.cellVectorVariables = getCellVectorVariables();
    data.particleVectorVariables = getParticleVectorVariables();
    data.expressions = getExpressions();
    return data;
}

void SimulationVisDataSourceImplPrivate::defineData()
{
    typedef VariableDefinitions::CellVectorVariable CellVectVar;
    //typedef VariableDefinitions::ParticleVectorVariable PVectVar;
    //typedef ExpressionDefinitions::Expression Expr;
    vector<CellVectVar>& cellVectVars = m_varDefinitions.cellVectorVariables;
    //vector<PVectVar>& pVectVars = m_varDefinitions.particleVectorVariables;
    //vector<Expr>& exprs = m_exprDefinitions.expressions;
    // define cell-variables
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    cellVectVars.push_back(CellVectVar("B1", 3, false, false,CellFormulas::B1Formula()));
    if(Params::constantMagneticFieldProfile.size() > 0) {
        cellVectVars.push_back(CellVectVar("B0", 3, false, false,CellFormulas::B0Formula()));
    }
    cellVectVars.push_back(CellVectVar("v", 3, true, true,CellFormulas::vFormula()));
    cellVectVars.push_back(CellVectVar("n", 1, true, true,CellFormulas::nFormula()));
    cellVectVars.push_back(CellVectVar("vtot", 3, true, false,CellFormulas::vFormula()));
    cellVectVars.push_back(CellVectVar("ntot", 1, true, false,CellFormulas::nFormula()));
    if(Params::resistivityFunction.isDefined() == true) {
        cellVectVars.push_back(CellVectVar("Eta*J", 3, false, false,CellFormulas::EResistiveFormula()));
    }
    if(Params::electronPressure==true) {
        cellVectVars.push_back(CellVectVar("-nablaPe", 3,false,false,CellFormulas::EPolarizationFormula()));
    }
    cellVectVars.push_back(CellVectVar("Ue", 3, false, false,CellFormulas::UeFormula()));
    cellVectVars.push_back(CellVectVar("MacroParticles", 1, true, false,CellFormulas::MacroParticlesFormula()));
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    cellVectVars.push_back(CellVectVar("SaveParticleOrbitFlag", 1, false, false,CellFormulas::SaveParticleFlagFormula()));
#endif
#endif
    cellVectVars.push_back(CellVectVar("T", 1, true, true,CellFormulas::TFormula()));
    cellVectVars.push_back(CellVectVar("Ttot", 1, true, false,CellFormulas::TFormula()));
    cellVectVars.push_back(CellVectVar("-UexB", 3, false, false,CellFormulas::EConvectiveFormula()));
    // Before the averaged values can be saved Tgrid::end_average() needs to be called.
    if(Params::averaging == true) {
        cellVectVars.push_back(CellVectVar("B1_average", 3, false, false,CellFormulas::B1AveFormula()));
        cellVectVars.push_back(CellVectVar("ntot_average", 1, true, false,CellFormulas::nTotAveFormula()));
#ifdef SAVE_POPULATION_AVERAGES
        cellVectVars.push_back(CellVectVar("n_ave", 1, true, true,CellFormulas::nAveFormula()));
        cellVectVars.push_back(CellVectVar("v_ave", 3, true, true,CellFormulas::vAveFormula()));
#endif
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        cellVectVars.push_back(CellVectVar("sph_B1_average", 3, false, false,CellFormulas::sph_B1AveFormula()));
#endif
    }
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    cellVectVars.push_back(CellVectVar("sph_B1", 3, false, false,CellFormulas::sph_B1Formula()));
    if (Params::constantMagneticFieldProfile.size() > 0)  {
        cellVectVars.push_back(CellVectVar("sph_B0", 3, false, false,CellFormulas::sph_B0Formula()));
        cellVectVars.push_back(CellVectVar("sph_B_tot", 3, false, false,CellFormulas::sph_B_totFormula()));
    }
    cellVectVars.push_back(CellVectVar("sph_dBdt", 3, false, false,CellFormulas::sph_dBdt_Formula()));
    cellVectVars.push_back(CellVectVar("sph_E", 3, false, false,CellFormulas::sph_EFormula()));
    cellVectVars.push_back(CellVectVar("sph_v", 3, true, true,CellFormulas::sph_vFormula()));
    cellVectVars.push_back(CellVectVar("sph_n", 1, true, true,CellFormulas::sph_nFormula()));
    cellVectVars.push_back(CellVectVar("sph_rho", 1, true, true,CellFormulas::sph_rhoFormula()));
    cellVectVars.push_back(CellVectVar("sph_rho_q", 1, true, true,CellFormulas::sph_rhoqFormula()));
    cellVectVars.push_back(CellVectVar("sph_rho_q_tot", 1, true, true,CellFormulas::sph_rhoq_totFormula()));
    cellVectVars.push_back(CellVectVar("sph_vtot", 3, true, false,CellFormulas::sph_vFormula()));
    cellVectVars.push_back(CellVectVar("sph_ntot", 1, true, false,CellFormulas::sph_nFormula()));
    if (Params::electronPressure==true) {
        cellVectVars.push_back(CellVectVar("sph_-nablaPe", 3,false,false,CellFormulas::sph_EPolarizationFormula()));
    }
    cellVectVars.push_back(CellVectVar("sph_Ue", 3, false, false,CellFormulas::sph_UeFormula()));
    cellVectVars.push_back(CellVectVar("sph_J", 3, false, false,CellFormulas::sph_JFormula()));
    cellVectVars.push_back(CellVectVar("sph_Ji", 3, false, false,CellFormulas::sph_JiFormula()));
#endif
    // define particle-variables
    //pVectVars.push_back(PVectVar("popID", 1, false,ParticleFormulas::popIDFormula()));
    // define expressions
    //exprs.push_back(Expr("rhov", true, "$rho * $v"));
}

ConstSequenceHandle<VisData::Patch>
SimulationVisDataSourceImplPrivate::getPatches()
{
    vector<VisData::Patch> patches;
    for (unsigned int refLevel = 0; refLevel < m_patches.size(); ++refLevel)
        for (vector<Interval<3> >::iterator patch = m_patches[refLevel].begin();
             patch != m_patches[refLevel].end(); ++patch) {
            /*cerr << "in: getPatches():loop, m_patches[refLevel].size:"
             << m_patches[refLevel].size() << " , refLevel: " << refLevel
             << endl;*/
            VisData::Patch p;
            p.refinementRatio = twoPower(refLevel);
            p.startCoordinates = SequenceHandle<std::size_t>
                                 (new ArraySequence<std::size_t>
                                  (patch->getBegin()[0], patch->getBegin()[1], patch->getBegin()[2]));
            p.endCoordinates = SequenceHandle<std::size_t>
                               (new ArraySequence<std::size_t>
                                (patch->getEnd()[0] - 1,
                                 patch->getEnd()[1] - 1,
                                 patch->getEnd()[2] - 1));
            p.cells = ConstSequenceHandle<VisData::Patch::Cell>
                      (new CachedFunctionSequence<VisData::Patch::Cell, VisCellGenerator>
                       (VisCellGenerator(), volume(*patch)));
            patches.push_back(p);
        }
    return ConstSequenceHandle<VisData::Patch>
           (new ArraySequence<VisData::Patch>(patches));
}


ConstSequenceHandle<VisData::Particle>
SimulationVisDataSourceImplPrivate::getParticles()
{
    ParticleMapper<VisData::Particle, ParticleDataGetter>::ParticleData pData
        = { m_cells };
    ConstSequenceHandle<VisData::Particle> particles =
        ConstSequenceHandle<VisData::Particle>
        ((new FunctionSequence<VisData::Particle,
          ParticleMapper<VisData::Particle, ParticleDataGetter> >
          (ParticleMapper<VisData::Particle, ParticleDataGetter>
           (pData, ParticleDataGetter()), m_grid.n_particles)));
    return particles;
}

ConstSequenceHandle<VisData::VectorVariable>
SimulationVisDataSourceImplPrivate::getCellVectorVariables()
{
    vector<string> hcFilePrefix;
    vector< vector<int> > popId;
    Population::getHcFileConfigs(hcFilePrefix, popId);
    CellVectorVariableGenerator::CellVariableData data
    (m_cells, m_patches, popId, getCellVectorVariableData(m_varDefinitions));
    return ConstSequenceHandle<VisData::VectorVariable>
           (new FunctionSequence<VisData::VectorVariable,
            CellVectorVariableGenerator>
            (CellVectorVariableGenerator(data), data.vectorVariables.size()));
}

ConstSequenceHandle<VisData::VectorVariable>
SimulationVisDataSourceImplPrivate::getParticleVectorVariables()
{
    ParticleVectorVariableGenerator::ParticleVariableData data
    (m_cells, m_grid.n_particles,
     getParticleVectorVariableData(m_varDefinitions));
    return ConstSequenceHandle<VisData::VectorVariable>
           (new CachedFunctionSequence<VisData::VectorVariable,
            ParticleVectorVariableGenerator>(ParticleVectorVariableGenerator(data),
                    data.vectorVariables.size()));
}

ConstSequenceHandle<VisData::Expression>
SimulationVisDataSourceImplPrivate::getExpressions()
{
    vector<string> hcFilePrefix;
    vector< vector<int> > popId;
    Population::getHcFileConfigs(hcFilePrefix, popId);
    vector<VisData::Expression> exprs;
    for (vector<ExpressionDefinitions::Expression>::iterator expr
         = m_exprDefinitions.expressions.begin();
         expr != m_exprDefinitions.expressions.end(); ++expr)
        if (!expr->isPopulationExpression) {
            VisData::Expression e = {expr->name, expr->formula};
            exprs.push_back(e);
        } else {
            for (vector<string>::iterator prefix = hcFilePrefix.begin();
                 prefix != hcFilePrefix.end(); ++prefix) {
                VisData::Expression e = {*prefix + "_" + expr->name,
                                         replace(expr->formula, "$", *prefix + "_")
                                        };
                exprs.push_back(e);
            }
        }
    return ConstSequenceHandle<VisData::Expression>
           (new ArraySequence<VisData::Expression>(exprs));
}

SimulationVisDataSourceImpl::SimulationVisDataSourceImpl(const Params& p,
        const Tgrid& g)
    : pimpl(new SimulationVisDataSourceImplPrivate(p, g)) { }

VisData SimulationVisDataSourceImpl::getVisData() const
{
    return pimpl->getVisData();
}

