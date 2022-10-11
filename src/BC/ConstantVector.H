#ifndef BC_VECTOR_H_
#define BC_VECTOR_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_Array.H>

#include "Set/Set.H"
#include "BC/BC.H"
#include "Numeric/Interpolator/Linear.H"

namespace BC
{
class ConstantVector
    : public BC<Set::Vector>
{

    #if AMREX_SPACEDIM==2
    enum Face{
        XLO, YLO, XHI, YHI, 
        INT
    };
    #elif AMREX_SPACEDIM==3
    enum Face{
        XLO, YLO, ZLO, XHI, YHI, ZHI, // 6
        INT
    };
    #endif

public:
    //Constant (amrex::Vector<amrex::Geometry> &_geom) : geom(_geom) {};
    ConstantVector() {} ;
    ConstantVector(IO::ParmParse &pp, std::string name) 
    {pp.queryclass(name,*this);} ;

    virtual ~ConstantVector() {};

    virtual void FillBoundary (amrex::BaseFab<Set::Vector> &in, const amrex::Box &box,
                    int ngrow, int dcomp, int ncomp, amrex::Real time,
                    Orientation face = Orientation::All,
                    const amrex::Mask *mask = nullptr) override;

    using BC::FillBoundary;
    
    amrex::BCRec GetBCRec() override;
    virtual amrex::Array<int,AMREX_SPACEDIM> IsPeriodic() override;
    virtual amrex::Periodicity Periodicity () const override;
    virtual amrex::Periodicity Periodicity (const amrex::Box& b) override;



    template<class T>
    const amrex::Array<amrex::Array<T,AMREX_SPACEDIM>,2> GetBCTypes()
    {
        return {{{AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])},
                {AMREX_D_DECL((T)m_bc_type[Face::XLO][0],(T)m_bc_type[Face::YLO][0],(T)m_bc_type[Face::ZLO][0])}}};
    }


private:
    #if AMREX_SPACEDIM==2
    static const int m_nfaces = 4;
    #elif AMREX_SPACEDIM==3
    static const int m_nfaces = 6;
    #endif

    unsigned int m_ncomp = AMREX_SPACEDIM;

    //int bc_lo[BL_SPACEDIM];
    //int bc_hi[BL_SPACEDIM];
    //amrex::Vector<amrex::Real> AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3);
    //amrex::Vector<amrex::Real> AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3);

    std::array<std::vector<int>, m_nfaces> m_bc_type; 
    std::array<Set::Vector,      m_nfaces> m_bc_val; 

public:
    static void Parse(ConstantVector & value, IO::ParmParse & pp)
    {
        std::map<std::string, int> bcmap;
            bcmap["BOGUS_BC"]         = BOGUS_BC;
            bcmap["INT_DIR"]          = INT_DIR;
            bcmap["REFLECT_ODD"]      = REFLECT_ODD;
            bcmap["INT_DIR"]          = INT_DIR;
            bcmap["REFLECT_EVEN"]     = REFLECT_EVEN;
            bcmap["FOEXTRAP"]         = FOEXTRAP;
            bcmap["EXT_DIR"]          = EXT_DIR;
            bcmap["HOEXTRAP"]         = HOEXTRAP;
            bcmap["Interior"]         = Interior;
            bcmap["Inflow"]           = Inflow;
            bcmap["Outflow"]          = Outflow;
            bcmap["Symmetry"]         = Symmetry;
            bcmap["SlipWall"]         = SlipWall;
            bcmap["NoSlipWall"]       = NoSlipWall;
        // From <AMReX_LO_BCTYPES.H>
        bcmap["interior"]         = (int)amrex::LinOpBCType::interior;
        bcmap["Dirichlet"]        = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["dirichlet"]        = (int)amrex::LinOpBCType::Dirichlet;
        bcmap["Neumann"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["NEUMANN"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["neumann"]          = (int)amrex::LinOpBCType::Neumann;
        bcmap["reflect_odd"]      = (int)amrex::LinOpBCType::reflect_odd;
        bcmap["Marshak"]          = (int)amrex::LinOpBCType::Marshak;
        bcmap["SanchezPomraning"] = (int)amrex::LinOpBCType::SanchezPomraning;
        bcmap["inflow"]           = (int)amrex::LinOpBCType::inflow;
        bcmap["Periodic"]         = (int)amrex::LinOpBCType::Periodic;
        bcmap["periodic"]         = (int)amrex::LinOpBCType::Periodic;


        
        value.m_bc_type[Face::XLO].clear(); value.m_bc_val[Face::XLO];
        value.m_bc_type[Face::XHI].clear(); value.m_bc_val[Face::XHI];
        value.m_bc_type[Face::YLO].clear(); value.m_bc_val[Face::YLO];
        value.m_bc_type[Face::YHI].clear(); value.m_bc_val[Face::YHI];
        #if AMREX_SPACEDIM == 3
        value.m_bc_type[Face::ZLO].clear(); value.m_bc_val[Face::ZLO];
        value.m_bc_type[Face::ZHI].clear(); value.m_bc_val[Face::ZHI];
        #endif
        
        // TYPES
        
        std::vector<std::string> str;
        pp.queryarr("type.xlo",str); 
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::XLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.xhi",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::XHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::XHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.ylo",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::YLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.yhi",str);
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::YHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::YHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        pp.queryarr("type.zlo",str);
        #if AMREX_SPACEDIM==3
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZLO].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::ZLO].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        #endif
        pp.queryarr("type.zhi",str);
        #if AMREX_SPACEDIM==3
            for (unsigned int i = 0; i < str.size(); i++) if (!bcmap.count(str[i])) Util::Abort(INFO,"Invalid BC: ", str[i]);
            if (str.size() == value.m_ncomp) for (unsigned int i = 0; i < value.m_ncomp; i++) value.m_bc_type[Face::ZHI].push_back(bcmap[str[i]]); 
            else if (str.size() == 1) value.m_bc_type[Face::ZHI].resize(value.m_ncomp,bcmap[str[0]]);
            else Util::Abort(INFO,"Incorrect number of ",pp.prefix()," BC type args: expected ", value.m_ncomp, " or 1 but got ", str.size());
        #endif

        // VALS
        std::vector<Set::Scalar> val;
        pp.queryarr("val.xlo",value.m_bc_val[Face::XLO]); 
        pp.queryarr("val.xhi",value.m_bc_val[Face::XHI]); 
        pp.queryarr("val.ylo",value.m_bc_val[Face::YLO]); 
        pp.queryarr("val.yhi",value.m_bc_val[Face::YHI]);
        pp.contains("val.zlo"); 
        pp.contains("val.zhi"); 
        #if AMREX_SPACEDIM==3
        pp.queryarr("val.zlo",value.m_bc_val[Face::ZLO]); 
        pp.queryarr("val.zhi",value.m_bc_val[Face::ZHI]); 
        #endif
    }

};
}
#endif
