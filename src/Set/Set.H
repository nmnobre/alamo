#ifndef SET_SET_
#define SET_SET_

#include <iomanip>

#include "Util/Util.H"

#include "Set/Base.H"
/// \brief A collection of data types and symmetry-reduced data structures

namespace Set
{
//class Field : public amrex::MultiFab
//{

//};
template <class T>
class Field : public amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>
{
public:
    Field() {} 
    Field(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost) 
    {
        Define(a_levs,a_grids,a_dmap,a_ncomp,a_nghost);
    } 
    Field(int size) : amrex::Vector<std::unique_ptr<amrex::FabArray<amrex::BaseFab<T>>>>(size) {}
    void Define(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_levs == a_grids.size()));
        Util::Assert(INFO,TEST(a_levs == a_dmap.size()));
        this->resize(a_levs);
        for (int lev = 0; lev < a_levs; lev++)
        {
            (*this)[lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(a_grids[lev],a_dmap[lev],a_ncomp,a_nghost));
        }
    }
    void Define(int a_lev, const amrex::BoxArray & a_grid, const amrex::DistributionMapping & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_lev < this->size()));
        (*this)[a_lev].reset(new amrex::FabArray<amrex::BaseFab<T>>(a_grid,a_dmap,a_ncomp,a_nghost));
    }
    int finest_level = 0;
};

template <>
class Field<Set::Scalar> : public amrex::Vector<std::unique_ptr<amrex::MultiFab>>
{
public:
    Field() {} 
    Field(int size) : amrex::Vector<std::unique_ptr<amrex::MultiFab>>(size) {}
    void Define(int a_levs, const amrex::Vector<amrex::BoxArray> & a_grids, const amrex::Vector<amrex::DistributionMapping> & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_levs == a_grids.size()));
        Util::Assert(INFO,TEST(a_levs == a_dmap.size()));
        this->resize(a_levs);
        for (int lev = 0; lev < a_levs; lev++)
        {
            (*this)[lev].reset(new amrex::MultiFab(a_grids[lev],a_dmap[lev],a_ncomp,a_nghost));
        }
    }
    void Define(int a_lev, const amrex::BoxArray & a_grid, const amrex::DistributionMapping & a_dmap, int a_ncomp, int a_nghost)
    {
        Util::Assert(INFO,TEST(a_lev < this->size()));
        (*this)[a_lev].reset(new amrex::MultiFab(a_grid,a_dmap,a_ncomp,a_nghost));
    }
    int finest_level = 0;
};
}

namespace Util
{
Set::Scalar Random();
}

namespace Set
{
namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
static const Set::Scalar INVSQ2 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105351942519376716382078636750692311545614851246241802792536860632206074854996791570661133296375279637789997525057639103028574;
static const Set::Scalar INVSQ3 = 0.57735026918962576450914878050195745564760175127012687601860232648397767230293334569371539558574952522520871380513556767665664836499965082627055183736479121617603107730076852735599160670036155830775501;

    AMREX_FORCE_INLINE
    Set::Scalar PermutationTensor(int i, int j, int k)
    {
        if (i==0 && j==1 && k==2) return 1.0;
        else if (i==2 && j==0 && k==1) return 1.0;
        else if (i==1 && j==2 && k==0) return 1.0;
        else if (i==2 && j==1 && k==0) return -1.0;
        else if (i==0 && j==2 && k==1) return -1.0;
        else if (i==1 && j==0 && k==2) return -1.0;
        else return 0.0;
    }
}
}

#include "Set/Matrix3.H"
#include "Set/Matrix4.H"


#endif
