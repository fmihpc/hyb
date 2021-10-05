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

#ifndef GRID_H
#define GRID_H

#include <fstream>
#include <cstring>
#include "definitions.h"
#include "particle.h"
#include "atmosphere.h"
#include "refinement.h"
#include "resistivity.h"
#include "forbidsplitjoin.h"
#include "backgroundcharge.h"
#include "magneticfield.h"

//! Magnetic field log
struct MagneticLog {
    real avgBx;
    real avgBy;
    real avgBz;
    real avgB;
    real avgDivB;
    real energyB;
    real maxDivB;
    real maxB;
    gridreal posMaxB[3];
    real maxDxDivBperB;
    MagneticLog();
};

//! Counters for field quantities
struct FieldCounter {
    real cutRateE;
    real cutRateRhoQ;
    real cutRateUe;
    real cutRateMaxVw;
    int resetTimestep;
    FieldCounter();
    void reset();
    void finalize();
};

//! Simulation grid
class Tgrid
{
public:
#ifndef USE_SPHERICAL_COORDINATE_SYSTEM
    //! Grid cell face average quantities
    enum TFaceDataSelect {FACEDATA_B=0, FACEDATA_J=1, FACEDATA_BSTAR=2, FACEDATA_MINUSDB=3, FACEDATA_AVEB=4};
    //! Grid cell node quantities
    enum TNodeDataSelect {NODEDATA_B=0, NODEDATA_UE=1, NODEDATA_E=2, NODEDATA_J=3, NODEDATA_Ji=4, NODEDATA_ne=5};
    //! Grid cell volume average quantities
    enum TCellDataSelect {CELLDATA_Ji=0, CELLDATA_UE=1, CELLDATA_J=2, CELLDATA_B=3, CELLDATA_TEMP1=4, CELLDATA_TEMP2=5};
    enum {NFACEDATA=5};
    enum {NCELLDATA=6};
    enum {NNODEDATA=6};
#else
    //! (SPHERICAL) Grid cell face average quantities
    enum TFaceDataSelect {FACEDATA_B=0, FACEDATA_J=1, FACEDATA_BSTAR=2, FACEDATA_MINUSDB=3, FACEDATA_AVEB=4, FACEDATA_UE=5};
    //! (SPHERICAL) Grid cell node quantities
    enum TNodeDataSelect {NODEDATA_B=0, NODEDATA_UE=1, NODEDATA_E=2, NODEDATA_J=3, NODEDATA_Ji=4, NODEDATA_ne=5, NODEDATA_TEMP2=6};
    //! (SPHERICAL) Grid cell volume average quantities
    enum TCellDataSelect {CELLDATA_Ji=0, CELLDATA_UE=1, CELLDATA_J=2, CELLDATA_B=3, CELLDATA_TEMP1=4, CELLDATA_TEMP2=5, CELLDATA_E=6, CELLDATA_J_TEMP=7, CELLDATA_B_TEMP1=8, CELLDATA_B_TEMP2=9};
    enum {NFACEDATA=6};
    enum {NNODEDATA=7};
    enum {NCELLDATA=10};
#endif
    //! Get nth bit of word (n=0 is leftmost)
    static bool GetBit(int word, int n) {
        return (word & (1 << n)) != 0;
    }
    //! Return word with nth bit set to b (b=0,1)
    static int SetBit(int word, int n, bool b) {
        if (b)
            return word | (1 << n);
        else
            return word & ~(1 << n);
    }
    static FieldCounter fieldCounter; //!< Field counter
    friend class GridPrivatesHelper;
private:
#define PUBLIC_TOBJECT
    struct Tcell; //!< Grid cell
    struct Trefintf; //!< Grid cell refinement interface
    struct Tface; //!< Grid cell face
    struct Tnode; //!< Grid cell node
    typedef Tnode *TNodePtr; //! Grid node pointer
    typedef Tface *TFacePtr; //! Grid face pointer
    typedef Tcell *TCellPtr; //! Grid cell pointer
    //! Grid cell
    struct Tcell PUBLIC_TOBJECT {
        bool haschildren; //!< Has the cell child cells?
        bool refine_it; //!< Refine this cell
        bool recoarsen_it; //!< Recoarsen this cell
        bool forbid_psplit; //!< Forbid split&join
        int refstatus;
        bool isrefined_face(int dim, int dir) const {
            return GetBit(refstatus,2*dim+dir);
        }
        void setrefined_face(int dim, int dir, bool flag) {
            refstatus = SetBit(refstatus,2*dim+dir,flag);
        }
        bool anyrefined_face() const {
            return refstatus != 0;
        }
        int level;
#ifdef SAVE_PARTICLES_ALONG_ORBIT
        bool save_particles; //!< If to save particles (along orbit)
#endif
        TCellPtr parent; //!< Pointer to parent cell, or 0 for level=0 cells
        TCellPtr neighbour[3][2]; //!< Pointers to direct neighbouring cells (without refinement)
        int flatind; //!< Flatindex of the cell (root cell), or childorder (0..7) for non-root cell
        int running_index; //!< Used when saving to file only, need not even be initialized
        //! Def: cell is root cell iff parent==0.
        union {
            TFacePtr face[3][2]; //!< Leaf cell: dim=0,1,2 (x/y/z), direction=0,1 (left/right)
            Trefintf *refintf[3][2]; //!< Leaf cell: dim=0,1,2 (x/y/z), direction=0,1 (left/right)
        };
        TCellPtr child[2][2][2]; //!< Nonleaf cell: Pointers to children (x/y/z=0/1)
        datareal nc; //!< Number density of the physical particles inside the cell [#/m^3]
        datareal rho_q; //!< Charge density inside the cell [C/m^3]
        datareal rho_q_bg; //!< Background charge density inside the cell [C/m^3]
        datareal ave_nc; //!< Temporally averaged density inside the cell [#/m^3]
#ifdef SAVE_POPULATION_AVERAGES
        std::vector<datareal> pop_ave_n; //!< Temporal population average: density
        std::vector<datareal> pop_ave_vx; //!< Temporal population average: vx
        std::vector<datareal> pop_ave_vy; //!< Temporal population average: vy
        std::vector<datareal> pop_ave_vz; //!< Temporal population average: vz
#endif
#ifdef SAVE_PARTICLE_CELL_SPECTRA
        std::vector< std::vector<datareal> > spectra; //!< Particle cell energy spectra
#endif
        datareal celldata[NCELLDATA][3]; //!< Cell data
        gridreal centroid[3]; //!< Centroid coordinates of the cell
        gridreal r2; //!< Square of the distance to the box origin [m^2]
        gridreal size; //!< Side length of the cell [m]
        gridreal invsize; //!< 1.0/size [1/m]
        TParticleList plist; //!< List of macroparticles residing in this cell
        real faceave(int dim, int dir, TFaceDataSelect s) const;
        void childave(TCellDataSelect cs, real result[3]) const;
        real childave_rhoq() const;
        real childave_nc() const;
        void FC_recursive(TFaceDataSelect fs, TCellDataSelect cs);
        void NC_smoothing_recursive();
        void NC_recursive(TNodeDataSelect ns,TCellDataSelect cs);
        void CN_recursive(TCellDataSelect cs, TNodeDataSelect ns);
        void CN_smoothing_recursive();
        void CN_rhoq_recursive();
        void CN_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
        void zero_rhoq_nc_Vq_recursive();
        void calc_ue_recursive(void);
        void calc_node_E_recursive(void);
        void calc_cell_E_recursive(void);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        void sph_calc_cell_E_recursive(void);
#endif
        void FaceCurl_recursive(TNodeDataSelect nsB, TFaceDataSelect fsj, int d, real factor);
        void NF_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d);
        void NF_rhoq_recursive(int d);
        void FacePropagate_recursive(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt);
        void set_B_recursive(void (*)(const gridreal[3], datareal[3]));
        void set_bgRhoQ_recursive(BackgroundChargeDensityProfile func);
        void calc_facediv_recursive(TFaceDataSelect fs, MagneticLog& result) const;
        void CalcGradient_rhoq_recursive();
        int Ncells_recursive() const;
        int Nfaces() const;
        int Nparticles_recursive() const;
        void enum_children_recursive();
        void writeMHD_children_recursive(std::ostream& o,const int filetype,std::vector<int> popId) const;
        void writeMHD(std::ostream& o,const int filetype,std::vector<int> popId) const;
        void writeDBUG_children_recursive(std::ostream& o) const;
        void writeDBUG(std::ostream& o) const;
        void writeEXTRA_children_recursive(std::ostream& o,ScalarField* s) const;
        void writeEXTRA(std::ostream& o,ScalarField* s) const;
#ifdef SAVE_PARTICLE_CELL_SPECTRA
        void writeSPECTRA_children_recursive(std::ostream& o,std::vector<int> popId) const;
        void writeSPECTRA(std::ostream& o,std::vector<int> popId) const;
#endif
        int mark_refinement_recursive(GridRefinementProfile refFunc);
        int mark_recoarsening_recursive(gridreal (*mindx)(const gridreal[]));
        void refine_recursive(Tgrid& g);
        void recoarsen_recursive(Tgrid& g);
        bool refine(Tgrid& g);
        bool recoarsen(Tgrid& g);
        template <class Func> int particle_pass_recursive(Func& op, bool relocate);
        int particle_pass_recursive(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), bool relocate);
        template <class Func> void cellPassRecursive(Func& op);
        void split_and_join_recursive(int& nsplit, int& njoined);
        int forbid_split_and_join_recursive(ForbidSplitAndJoinProfile forb);
        void begin_average_recursive();
        void end_average_recursive(real inv_ave_ntimes);
        void set_resistivity_recursive(ResistivityProfile res);
        void prepare_PDF_recursive(ScalarField* pdffunc, Tgrid*);
        void generate_random_point(gridreal r[3]);
        void cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, std::vector<int> popId);
#ifdef SAVE_PARTICLES_ALONG_ORBIT
        void particles_write_recursive();
#endif
        void CN_ne_recursive();
        void calc_node_j_recursive();
        void calc_node_ue_recursive();
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        int sph_CN_boundary_flag; //!< (SPHERICAL) 
        gridreal sph_sizey; //!< (SPHERICAL) 
        gridreal sph_sizez; //!< (SPHERICAL) 
        gridreal sph_invsizey; //!< (SPHERICAL) 
        gridreal sph_invsizez; //!< (SPHERICAL) 
        gridreal sph_centroid[3]; //!< (SPHERICAL) 
        gridreal sph_coor[3]; //!< (SPHERICAL) 
        gridreal sph_coor_next[3]; //!< (SPHERICAL) 
        gridreal sph_size[3]; //!< (SPHERICAL) 
        gridreal sph_dl[3]; //!< (SPHERICAL) 
        gridreal sph_diag; //!< (SPHERICAL) 
        gridreal sph_dS[3]; //!< (SPHERICAL) 
        gridreal sph_dS_centroid[3]; //!< (SPHERICAL) 
        gridreal sph_dS_next[3]; //!< (SPHERICAL) 
        gridreal sph_dV; //!< (SPHERICAL) Cell volume
        void sph_FC_recursive(TFaceDataSelect fs, TCellDataSelect cs, int FC_boundary_flag);
        void FC_copy_recursive(TFaceDataSelect fs, TCellDataSelect cs);
        void sph_NC_smoothing_recursive();
        void sph_NC_recursive(TNodeDataSelect ns,TCellDataSelect cs);
        void NC_copy_recursive(TNodeDataSelect ns,TCellDataSelect cs);
        void sph_NC_copy_recursive(TNodeDataSelect ns,TCellDataSelect cs);
        void sph_CN_recursive(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CNb_recursive(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag);
        void sph_CN_vol_recursive(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CN_Ji_recursive(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CNb_Ji_recursive(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag);
        void C_2_recursive(TCellDataSelect cs);
        void sph_Cell_BC_recursive(TCellDataSelect cs);
        void sph_CN_smoothing_recursive();
        void sph_CNb_smoothing_recursive(int sph_CN_boundary_flag);
        void sph_CN_rhoq_recursive();
        void sph_CNb_rhoq_recursive(int sph_CN_boundary_flag);
        void sph_CN_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
        void sph_CNb_donor_recursive(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt, int sph_CN_boundary_flag);
        void sph_calc_ue_recursive(void);
        void sph_calc_ue_app_recursive(void);
        void sph_calc_node_E_recursive(void);
        void sph_Node_BC_recursive(TNodeDataSelect ns);
        void sph_calc_node_E_app_recursive(void);
        void sph_FaceCurl_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor);
        void sph_FaceCurl_2_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor);
        void sph_NF_recursive(TNodeDataSelect ns, TFaceDataSelect fs, int d);
        void sph_NF_rhoq_recursive(int d);
        void sph_calc_facediv_recursive(TFaceDataSelect fs, TCellDataSelect cs, MagneticLog& result) const;
        void sph_CalcGradient_rhoq_recursive();
        void sph_generate_random_point(gridreal r[3]);
        void sph_cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, std::vector<int> popId);
        void sph_CN_ne_recursive();
        void sph_CNb_ne_recursive(int sph_CN_boundary_flag);
        void sph_calc_node_j_recursive();
        void sph_calc_node_j_2_recursive();
        void sph_calc_node_ue_recursive();
#endif
    private:
        void take27neighbours(TCellPtr celltab[3][3][3]);
    }; //! Grid cell
    //! Grid cell refinement interface, used when grid cell size changes
    struct Trefintf PUBLIC_TOBJECT {
        TFacePtr face[4]; //!< Pointers to the four faces in (-y,-z),(+y,-z),(-y,+z),(+y,+z) order (first-y-then-z order)
    };
    //! Grid cell face
    struct Tface PUBLIC_TOBJECT {
        TNodePtr node[4]; //!< Pointers to the four corner nodes in (-y,-z),(+y,-z),(+y,+z),(-y,+z) order (cyclic order)
        TNodePtr sidenode[4]; //!< sidenode[0] is between node[0] and node[1], sidenode[1] between node[1] and node[2], etc.
        datareal facedata[NFACEDATA]; //!< B and j (FACEDATA_B, FACEDATA_J)
        Tface() {
            memset(facedata,0,sizeof(datareal)*NFACEDATA);
            node[0] = node[1] = node[2] = node[3] = 0;
            sidenode[0] = sidenode[1] = sidenode[2] = sidenode[3] = 0;
        }
        void Curl1(TNodeDataSelect nsB, TFaceDataSelect fsj, int d, gridreal dx, real factor);
        void NF1(TNodeDataSelect ns, TFaceDataSelect fs, int d);
        void NF_rhoq1();
        void Propagate1(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt);
        //! Set magnetic field (B1) on a face
        void set_B1(const gridreal r[3], void (*f)(const gridreal[3], datareal[3]), int d) {
            datareal b[3];
            (*f)(r, b);
            facedata[FACEDATA_B] = b[d];
        }
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        void sph_Curl1(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor);
        void sph_Curl1_2(TNodeDataSelect ns, TFaceDataSelect fs, int d, real factor);
        void sph_NF1(TNodeDataSelect ns, TFaceDataSelect fs, int d);
        void sph_NF_rhoq1(int d);
#endif
    }; // Grid cell face
    //! Grid cell node
    struct Tnode PUBLIC_TOBJECT {
        TCellPtr cell[2][2][2]; //!< Each node touches 8 cells, pointers to them (if fewer, the rest are null)
        gridreal centroid[3]; //!< Node coordinates (exact)
        gridreal r2; //!< Square of the distance to the box origin [m^2]
        gridreal r0; //!< Distance to the box origin [m^2]
        datareal nodedata[NNODEDATA][3]; //< Cells are ordered -+x, -+y, -+z, for example +x,-y,-z cell is [1][0][0]
        datareal nn; //!< Density at the node
        datareal eta; //!< Resistivity at the node
        //! Constructor
        Tnode() {
            memset(&nodedata[0][0],0,sizeof(datareal)*NNODEDATA*3);
            nn = 0;
            eta = 0;
        }
        void CN1(TCellDataSelect cs, TNodeDataSelect ns);
        void CN_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
        void CN1_density();
        void CN1_smoothing();
        void CN1_rhoq();
        void set_resistivity(ResistivityProfile res);
        void calc_E1(void);
        void update_cell_pointers(gridreal size, Tgrid& g);
        real jcomponent(int edir, real dx);
        void calc_j(real dx);
        void CN1_ne();
        void calc_ue1();
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        gridreal sph_centroid[3];
        void sph_CN1(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CNb1(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag);
        void sph_CN1_vol(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CN1_Ji(TCellDataSelect cs, TNodeDataSelect ns);
        void sph_CNb1_Ji(TCellDataSelect cs, TNodeDataSelect ns, int sph_CN_boundary_flag);
        void sph_CN_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
        void sph_CNb_donor1(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt, int sph_CN_boundary_flag);
        void sph_CN1_smoothing();
        void sph_CNb1_smoothing(int sph_CN_boundary_flag);
        void sph_CN1_rhoq();
        void sph_CNb1_rhoq(int sph_CN_boundary_flag);
        void sph_calc_E1(void);
        void sph_Node_BC1(TNodeDataSelect ns);
        void sph_calc_E1_app(void);
        real sph_jcomponent(int edir);
        real sph_jcomponent_2(int edir);
        void sph_calc_j();
        void sph_calc_j_2();
        void sph_CN1_ne();
        void sph_CNb1_ne(int sph_CN_boundary_flag);
        void sph_calc_ue1();
#endif
    }; // Grid cell node
    //! TPtrHash indexes void* pointers using gridreal triples.
    class TPtrHash
    {
    private:
        enum {PTR_HASHTABLE_SIZE = 39};
        struct Thashnode {
            void *data;
            int iq[3];
            Thashnode *next;
        };
        typedef Thashnode *Thashnodeptr;
        Tgrid *gptr; //!< Needed to access Tgrid::intcoords
        Thashnode **pool;
        int HashFunction(const gridreal[3], int iq[3]) const;
    public:
        TPtrHash(Tgrid& g) : gptr(&g) {
            pool = new Thashnodeptr [PTR_HASHTABLE_SIZE];
            memset(pool,0,sizeof(Thashnodeptr)*PTR_HASHTABLE_SIZE);
        }
        void *add(const gridreal r[3], void *value);
        void *add_unique(const gridreal r[3], void *value);
        //! Convenient abbreviation call
        void *add_unique(Tgrid::TNodePtr nodeptr) {
            return add_unique(nodeptr->centroid,(void*)nodeptr);
        }
        void *read(const gridreal r[3]) const;
        void remove(const gridreal r[3]);
        void delete_all_as_nodeptr();
        ~TPtrHash();
    };
    //! TBoxDef defines a cubical cell: its lower-left corner and size
    struct TBoxDef {
        gridreal lowx,lowy,lowz;
        gridreal size;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
        gridreal sph_sizey, sph_sizez;
#endif
    };
    //! Probability density function (PDF) object associated with Tgrid
    struct TPDFTable {
        int n; //!< Number of interior cells, length of tables
        shortreal *cpdf; //!< Cumulative PDF, monotonically increasing from 0 to 1, vector of length n
        TCellPtr *cellptrs; //!< Pointer to interior cell for each CPDF entry, vector of length n (computed only for the FIRST PDF allocated (pdftables[0]) because it is the same for subsequent ones)
    };

    // ---------------- Private data of Tgrid: ------------------
    
    int nx,ny,nz; //!< Basegrid size including ghost cells
    gridreal x_1, y_1, z_1; //!< Basegrid (including ghosts) lowest corner point
    gridreal bgdx,invbgdx; //!< Grid spacing (isotropic) and its inverse (invbgdx=1/bgdx)
    real inv_unit; //!< 1/(smallest representable unit wrt. gridreal "epsilon")
    TCellPtr *cells;
    TCellPtr saved_cellptr; //!< Routines which get r[3] as input saves the found cell here (avoids unnecessary findcell() call)
    const static char *celldata_names[NCELLDATA];
    static int cell_running_index; //!< Running cell index
    static TPtrHash *hp;
    int n_particles; //!< Number of macro particles
    int ave_ntimes; //!< Temporal averaging counter
    TCellPtr previous_found_cell;
    enum {MAX_PDFTABLES = 100};
    TPDFTable pdftables[MAX_PDFTABLES];
    int n_pdftables; //!< Length of entries in pdftables, initially 0
    int pdftable_iterator; //!< Used only by prepare_PDF, prepare_PDF_recursive
    
    // ---------------- Private functions of Tgrid: --------------

    //! Compute flat index of a cell
    int flatindex(int i, int j, int k) const {
        return (i*ny + j)*nz + k;
    }
    //! Decompose a flat index of a cell
    void decompose(int ijk, int& i, int& j, int& k) const {
        k = ijk % nz;
        const int ij = ijk / nz;
        j = ij % ny;
        i = ij / ny;
    }
    //! "Round" coordinate point to unique integer triple
    void intcoords(const gridreal r[3], int iq[3]) const {
        iq[0] = int((r[0] - x_1)*inv_unit + 0.5);
        iq[1] = int((r[1] - y_1)*inv_unit + 0.5);
        iq[2] = int((r[2] - z_1)*inv_unit + 0.5);
    }
    //! Inverse operation of intcoords
    void realcoords(const int iq[3], gridreal r[3]) const {
        r[0] = x_1 + iq[0]/inv_unit;
        r[1] = y_1 + iq[1]/inv_unit;
        r[2] = z_1 + iq[2]/inv_unit;
    }
    void update_neighbours(Tcell *c);
    //! Get cell size
    gridreal size(const Tcell *c) const {
        return c->size;
    }
    //! Get 1/(cell size)
    gridreal invsize(const Tcell *c) const {
        return c->invsize;
    }
    //! Get cell volume
    gridreal volume(const Tcell *c) const {
        const gridreal sz = size(c);
        return sz*sz*sz;
    }
    TCellPtr moveto(TCellPtr c, int dim, bool movetoright) const {
        return c->neighbour[dim][movetoright];
    }
    gridreal accumulate_PIC_recursive(const TBoxDef& cloudbox, Tcell *c, gridreal invvol, const shortreal v[3], real w, int popid);
    static gridreal intersection_volume(const TBoxDef& boxA, const TBoxDef& boxB);
    static gridreal intersection_volume_samesize_nochecks(const shortreal rA[3], const gridreal rB[3], gridreal size);
    void copy_celldata(int cTo, int cFrom, TCellDataSelect cs);
    void copy_rhoq(int cTo, int cFrom);
    void copy_smoothing(int cTo, int cFrom);
    void finalize_accum_recursive(Tcell *c);
    struct writeParticle;
    struct writeMagneticField;
    struct readMagneticField;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    gridreal sph_theta_1, sph_phi_1; //!< (SPHERICAL)
    gridreal sph_bgdy, sph_bgdz, sph_bgdtheta, sph_bgdphi, sph_invbgdy, sph_invbgdz; //!< (SPHERICAL)
    //! (SPHERICAL) Get size y
    gridreal sph_sizey(const Tcell *c) const {
        return c->sph_sizey;
    }
    //! (SPHERICAL) Get size z
    gridreal sph_sizez(const Tcell *c) const {
        return c->sph_sizez;
    }
    //! (SPHERICAL) Get volume
    gridreal sph_dV(const Tcell *c) const {
        return c->sph_dV;
    }
    //! (SPHERICAL) Get volume (check negative value)
    gridreal sph_volume(const Tcell *c) {
        gridreal dV = sph_dV(c);
        if (dV < 0.0) dV = -dV;
        return dV;
    }
    gridreal sph_accumulate_PIC_recursive(const TBoxDef& cloudbox, Tcell *c, gridreal invvol, const shortreal v[3], real w, int popid);
    static gridreal sph_intersection_volume(const TBoxDef& boxA, const TBoxDef& boxB);
    void copy_celldata_temp(int cTo, int cFrom, TCellDataSelect cs);
    void copy_celldata_SCS(int cTo, int cFrom, TCellDataSelect cs);
    void copy_celldata_ave(int cTo, int cFrom, TCellDataSelect cs);
    void sph_copy_rhoq(int cTo, int cFrom);
    void sph_copy_smoothing(int cTo, int cFrom);
    void sph_copy_smoothing_0(int cTo);
    void sph_finalize_accum_recursive(Tcell *c);
#endif
    // Tgrid
public:
    TCellPtr findcell(const shortreal r[3], gridreal* lowercorner=0);
    TCellPtr findcell_to_maxlevel(const shortreal r[3], int maxlevel) const;
    TParticleList *find_plist(const TLinkedParticle& p);
    Tgrid();
    ~Tgrid();
    void init(int nx1, int ny1, int nz1, gridreal x1, gridreal y1, gridreal z1, gridreal bgdx1);
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    void sph_init(int nx1, int ny1, int nz1, gridreal x1, gridreal y1, gridreal z1, gridreal theta1, gridreal phi1, gridreal bgdx1, gridreal bgdy1, gridreal bgdz1);
#endif
    void getbox(gridreal& x1, gridreal& x2, gridreal& y1, gridreal& y2, gridreal& z1, gridreal& z2) const {
        x1 = x_1+bgdx;
        y1 = y_1+bgdx;
        z1 = z_1+bgdx;
        x2 = x1 + (nx-2)*bgdx;
        y2 = y1 + (ny-2)*bgdx;
        z2 = z1 + (nz-2)*bgdx;
    }
    gridreal getbasedx() const {
        return bgdx;
    }
    int Ncells_with_ghosts() const;
    int Ncells_without_ghosts() const;
    int Nbasegrid_cells_without_ghosts() const {
        return (nx-2)*(ny-2)*(nz-2);
    }
    int Nbasegrid_cells_with_ghosts() const {
        return nx*ny*nz;
    }
    int Ncells() const {
        return Ncells_without_ghosts();
    }
    int Nbasegrid_cells() const {
        return Nbasegrid_cells_without_ghosts();
    }
    real BasegridSpacing() const {
        return bgdx;
    }
    void faceintpol(const shortreal r[3], TFaceDataSelect s, real result[3]);
    void cellintpol(const shortreal r[3], TCellDataSelect s, real result[3]);
    void cellintpol(TCellDataSelect s, real result[3]) const;
    void cellintpol_fluid(const shortreal r[3], real& n, real& vx, real& vy, real& vz, real& P, std::vector<int> popId);
    void cellintpol_fluid(real& n, real& vx, real& vy, real& vz, real& P, std::vector<int> popId);
    void FC(TFaceDataSelect fs, TCellDataSelect cs);
    void NC_smoothing();
    void NC(TNodeDataSelect ns, TCellDataSelect cs);
    void CN(TCellDataSelect cs, TNodeDataSelect ns);
    void CN_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
    void CN_smoothing();
    void CN_rhoq();
    void set_resistivity(ResistivityProfile res);
    void FaceCurl(TNodeDataSelect ns, TFaceDataSelect fs, real factor);
    void NF(TNodeDataSelect ns, TFaceDataSelect fs);
    void NF_rhoq();
    void zero_rhoq_nc_Vq();
    void accumulate_PIC(const shortreal r[3], const shortreal v[3], real w, int popid);
    void finalize_accum();
    void calc_ue(void);
    void Neumann(TCellDataSelect cs);
    void Neumann_rhoq();
    void Neumann_smoothing();
    void smoothing();
    void smoothing_E();
    void set_B(void (*)(const gridreal[3], datareal[3]));
    void set_bgRhoQ(BackgroundChargeDensityProfile func);
    void boundarypass(int dim, bool toRight, void (*)(datareal cdata[NCELLDATA][3], int));
    void calc_node_E(void);
    void calc_cell_E(void);
    void FacePropagate(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt);
#ifdef SAVE_PARTICLES_ALONG_ORBIT
    void set_save_particles_orbit(const char *fn);
    void particles_write();
#endif
    bool hcwrite_MHD(const char *fn,std::string ascbin,std::string hctype,std::vector<int> popId);
    bool hcwrite_DBUG(const char *fn);
    bool hcwrite_EXTRA(std::string fileName,ScalarField* s,std::string);
#ifdef SAVE_PARTICLE_CELL_SPECTRA
    bool hcwrite_SPECTRA(const char *fn,std::string ascbin,std::vector<int> popId);
#endif
    void dumpState(std::ostream& os);
    void readState(std::istream& is);
    void Refine(GridRefinementProfile refFunc);
    void recoarsen(gridreal (*mindx)(const gridreal[3]));
    void calc_facediv(TFaceDataSelect fs, MagneticLog& result) const;
    void CalcGradient_rhoq();
    int approx_bytes_per_cell() {
        return sizeof(Tcell) + sizeof(Tnode) + 3*sizeof(Tface);
    }
    static size_t bytes_allocated() {
        return 0;
    }
    void addparticle(shortreal x, shortreal y, shortreal z,
                     shortreal vx,shortreal vy,shortreal vz,
                     shortreal w, int popid, bool inject=true);
    template <class Func> int particle_pass(Func op, bool relocate=false);
    int particle_pass(bool (*op)(TLinkedParticle& p, ParticlePassArgs a), bool relocate=false);
    template <class Func> void cellPass(Func op);
    /** \brief Call operator for all particles in the grid
     *
     * If op returns false, delete the particle afterwards,
     * return number of deletions _with_ relocation: assume op
     * may have changed x,y,z, and assign particle to new
     * cell's plist if needed.
     */
    int particle_pass_with_relocation(bool (*op)(TLinkedParticle& p)) {
        return particle_pass(op,true);
    }
    int Nparticles() const;
    void split_and_join(int& nsplit, int& njoined);
    int forbid_split_and_join(ForbidSplitAndJoinProfile forb);
    void begin_average();
    bool end_average();
    void prepare_PDF(ScalarField* pdffunc, TPDF_ID& pdfid, real& cumsumvalue);
    void generate_random_point(const TPDF_ID& pdfid, gridreal r[3]);
    void boundary_faces(TFaceDataSelect cs);
    void CN_ne();
    void calc_node_j();
    void calc_node_ue();
    //! So that Tcell::refine(g) can access g->findcell()
    friend class Tcell;
    //! HashFunction needs to access Tgrid::intcoords
    friend class TPtrHash;
    friend class SimulationVisDataSourceImplPrivate;
#ifdef USE_SPHERICAL_COORDINATE_SYSTEM
    void sph_faceintpol(gridreal r[3], TFaceDataSelect s, real result[3]);
    void sph_FC(TFaceDataSelect fs, TCellDataSelect cs);
    void FC_copy(TFaceDataSelect fs, TCellDataSelect cs);
    void sph_NC_smoothing();
    void sph_NC(TNodeDataSelect ns, TCellDataSelect cs);
    void NC_copy(TNodeDataSelect ns,TCellDataSelect cs);
    void sph_NC_copy(TNodeDataSelect ns,TCellDataSelect cs);
    void sph_CN(TCellDataSelect cs, TNodeDataSelect ns);
    void sph_CNb(TCellDataSelect cs, TNodeDataSelect ns);
    void sph_CN_vol(TCellDataSelect cs, TNodeDataSelect ns);
    void sph_CN_Ji(TCellDataSelect cs, TNodeDataSelect ns);
    void sph_CNb_Ji(TCellDataSelect cs, TNodeDataSelect ns);
    void sph_YCave(TCellDataSelect cs);
    void sph_YNave(TNodeDataSelect ns);
    void sph_YNave_Ji(TNodeDataSelect ns);
    void C_2(TCellDataSelect cs);
    void sph_Cell_BC(TCellDataSelect cs);
    void sph_CN_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
    void sph_CNb_donor(TCellDataSelect cs, TNodeDataSelect ns, TNodeDataSelect uns, real dt);
    void sph_CN_smoothing();
    void sph_CNb_smoothing();
    void sph_CN_rhoq();
    void sph_CNb_rhoq();
    void sph_CNF_map(TCellDataSelect cs, TNodeDataSelect ns, TFaceDataSelect fs);
    void sph_FaceCurl(TNodeDataSelect ns, TFaceDataSelect fs, real factor);
    void sph_FaceCurl_2(TNodeDataSelect ns, TFaceDataSelect fs, real factor);
    void sph_NF(TNodeDataSelect ns, TFaceDataSelect fs);
    void sph_NF_rhoq();
    void sph_accumulate_PIC(shortreal r[3], const shortreal v[3], real w, int popid);
    void sph_finalize_accum();
    void sph_calc_ue(void);
    void sph_calc_ue_app(void);
    void sph_calc_node_E(void);
    void sph_Neumann_periodic(TCellDataSelect cs);
    void sph_Neumann(TCellDataSelect cs);
    void sph_Neumann_0(TCellDataSelect cs);
    void sph_Neumann_SCS(TCellDataSelect cs);
    void sph_Neumann_rhoq();
    void sph_Neumann_rhoq_0();
    void sph_Neumann_smoothing();
    void sph_Neumann_smoothing_0();
    void sph_smoothing();
    void sph_smoothing_E();
    void sph_boundarypass(int dim, bool toRight, void (*)(datareal cdata[NCELLDATA][3], gridreal sph_centroid[3], int));
    void sph_Node_BC(TNodeDataSelect ns);
    void sph_calc_node_E_app(void);
    void sph_calc_cell_E(void);
    void sph_calc_facediv(TFaceDataSelect fs, TCellDataSelect cs, MagneticLog& result) const;
    void sph_CalcGradient_rhoq();
    void sph_addparticle(shortreal x, shortreal y, shortreal z, shortreal vx,shortreal vy,shortreal vz, shortreal w, int popid);
    void sph_boundary_faces(TFaceDataSelect cs);
    void sph_CN_ne();
    void sph_CNb_ne();
    void sph_calc_node_j();
    void sph_calc_node_j_2();
    void sph_calc_node_ue();
#endif
}; // Tgrid

//! Intersection volume of two rectangular boxes
inline gridreal Tgrid::intersection_volume(const TBoxDef& boxA, const TBoxDef& boxB)
{
    gridreal V;
    V = min2(boxA.lowx+boxA.size, boxB.lowx+boxB.size) - max2(boxA.lowx, boxB.lowx);
    if (V <= 0) return 0;
    V*= min2(boxA.lowy+boxA.size, boxB.lowy+boxB.size) - max2(boxA.lowy, boxB.lowy);
    if (V <= 0) return 0;
    V*= min2(boxA.lowz+boxA.size, boxB.lowz+boxB.size) - max2(boxA.lowz, boxB.lowz);
    if (V <= 0) return 0;
    return V;
}

/** \brief Intersection volume of two rectangular boxes
 *
 * Assumes:
 * (1) that the boxes do intersect, i.e. that the volume is never zero
 * (2) that the boxes are of the same size
 */
inline gridreal Tgrid::intersection_volume_samesize_nochecks(const shortreal rA[3], const gridreal rB[3], gridreal size)
{
    gridreal V;
    V = size - fabs(rA[0]-rB[0]);
    V*= size - fabs(rA[1]-rB[1]);
    V*= size - fabs(rA[2]-rB[2]);
    return V;
}

//! Face propagate: Bnew = Bold - curl(E)*dt
inline void Tgrid::Tface::Propagate1(TFaceDataSelect Bold, TFaceDataSelect Bnew, real dt)
{
    facedata[Bnew] = facedata[Bold] - dt*facedata[FACEDATA_MINUSDB];
}

#ifdef USE_SPHERICAL_COORDINATE_SYSTEM

//! (SPHERICAL) Spherical version of "intersection_volume". Here we include the fact that the sizes in different directions may be different.
inline gridreal Tgrid::sph_intersection_volume(const TBoxDef& boxA, const TBoxDef& boxB)
{
    gridreal V;
    V = min2(boxA.lowx+boxA.size, boxB.lowx+boxB.size) - max2(boxA.lowx, boxB.lowx);
    if (V <= 0) return 0;
    V*= min2(boxA.lowy+boxA.sph_sizey, boxB.lowy+boxB.sph_sizey) - max2(boxA.lowy, boxB.lowy);
    if (V <= 0) return 0;
    V*= min2(boxA.lowz+boxA.sph_sizez, boxB.lowz+boxB.sph_sizez) - max2(boxA.lowz, boxB.lowz);
    if (V <= 0) return 0;
    return V;
}
#endif

#endif

