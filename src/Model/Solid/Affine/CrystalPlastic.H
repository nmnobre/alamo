#ifndef MODEL_SOLID_PLASTIC_CRYSTALPLASTIC_H_
#define MODEL_SOLID_PLASTIC_CRYSTALPLASTIC_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Cubic.H"

/*
    Note: everything is hardcoded for FCC right now. 
    We will deal with templating BCC and HCP later.
*/

namespace Model
{
namespace Solid
{
namespace Affine
{
class CrystalPlastic : public Cubic
{
public:
    CrystalPlastic(){};
    CrystalPlastic(Solid<Set::Sym::MajorMinor> base) : Cubic(base) {};
    CrystalPlastic(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Set::Scalar a_phi1, Set::Scalar a_Phi, Set::Scalar a_phi2)
    {
        Eigen::Matrix3d m;
    	m =     Eigen::AngleAxisd(a_phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(a_Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(a_phi1, Eigen::Vector3d::UnitX());
        Define(a_C11, a_C12, a_C44, m);
    }
    CrystalPlastic(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m)
    {
        Define(a_C11, a_C12, a_C44, m);
    }

    void Define(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Set::Scalar a_phi1, Set::Scalar a_Phi, Set::Scalar a_phi2)
    {
        Eigen::Matrix3d m;
    	m =     Eigen::AngleAxisd(a_phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(a_Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(a_phi1, Eigen::Vector3d::UnitX());
        Define(a_C11, a_C12, a_C44, m);
    }

    void Define(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m)
    {
        Affine::Cubic::Define(a_C11, a_C12, a_C44, m, Set::Matrix::Zero());
        
        // Set with default values
        DefineSlipSystem(16.0, 10.0, 0.001);
    }
    
    void Define (Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m, std::array<Set::Scalar, 12> &a_g_alpha, std::array<Set::Scalar, 12> &a_exp, std::array<Set::Scalar, 12> &a_gam0dot)
    {
        Affine::Cubic::Define(a_C11, a_C12, a_C44, m, Set::Matrix::Zero());
        DefineSlipSystem(a_g_alpha, a_exp, a_gam0dot);
    }

    void DefineSlipSystem(Set::Scalar a_g_alpha, Set::Scalar a_exp, Set::Scalar a_gam0got)
    {
        std::array<Set::Scalar, 12> f_g_alpha, f_m_exp, f_gam0dot;
        for(int i = 0; i<12; i++)
        {
            f_g_alpha[i] = a_g_alpha;
            f_m_exp[i] = a_exp;
            f_gam0dot[i] = a_gam0got;
        }
        DefineSlipSystem(f_g_alpha, f_m_exp, f_gam0dot);
    }

    void DefineSlipSystem(std::array<Set::Scalar, 12> &a_g_alpha, std::array<Set::Scalar, 12> &a_exp, std::array<Set::Scalar, 12> &a_gam0dot)
    {
        Set::Matrix R = Set::Matrix::Zero();
        AMREX_D_TERM(
            R(0,0) = m(0,0);
            ,
            R(0,1) = m(0,1); R(1,0) = m(1,0); R(1,1) = m(1,1);
            ,
            R(0,2) = m(0,2); R(1,2) = m(1,2); R(2,0) = m(2,0); R(2,1) = m(2,1); R(2,2) = m(2,2);
        );
        Set::Scalar f_dir = AMREX_D_PICK(1.0,1.0,INVSQ2);
        Set::Scalar f_nor = AMREX_D_PICK(1.0, INVSQ2, INVSQ3);
        
        plane[0] = R*Set::Vector(AMREX_D_DECL(f_nor, f_nor, f_nor));
        plane[1] = R*Set::Vector(AMREX_D_DECL(f_nor, f_nor, f_nor));
        plane[2] = R*Set::Vector(AMREX_D_DECL(f_nor, f_nor, f_nor));
        dir[0] = R*Set::Vector(AMREX_D_DECL(0.0, -f_dir, f_dir));
        dir[1] = R*Set::Vector(AMREX_D_DECL(f_dir, 0.0, -f_dir));
        dir[2] = R*Set::Vector(AMREX_D_DECL(-f_dir, f_dir, 0.0));

        plane[3] = R*Set::Vector(AMREX_D_DECL(-f_nor, -f_nor, f_nor));
        plane[4] = R*Set::Vector(AMREX_D_DECL(-f_nor, -f_nor, f_nor));
        plane[5] = R*Set::Vector(AMREX_D_DECL(-f_nor, -f_nor, f_nor));
        dir[3] = R*Set::Vector(AMREX_D_DECL(0.0, f_dir, f_dir));
        dir[4] = R*Set::Vector(AMREX_D_DECL(-f_dir, 0.0, -f_dir));
        dir[5] = R*Set::Vector(AMREX_D_DECL(f_dir, -f_dir, 0.0));

        plane[6] = R*Set::Vector(AMREX_D_DECL(-f_nor, f_nor, f_nor));
        plane[7] = R*Set::Vector(AMREX_D_DECL(-f_nor, f_nor, f_nor));
        plane[8] = R*Set::Vector(AMREX_D_DECL(-f_nor, f_nor, f_nor));
        dir[6] = R*Set::Vector(AMREX_D_DECL(0.0, -f_dir, f_dir));
        dir[7] = R*Set::Vector(AMREX_D_DECL(-f_dir, 0.0, -f_dir));
        dir[8] = R*Set::Vector(AMREX_D_DECL(f_dir, f_dir, 0.0));

        plane[9] = R*Set::Vector(AMREX_D_DECL(f_nor, -f_nor, f_nor));
        plane[10] = R*Set::Vector(AMREX_D_DECL(f_nor, -f_nor, f_nor));
        plane[11] = R*Set::Vector(AMREX_D_DECL(f_nor, -f_nor, f_nor));
        dir[9] = R*Set::Vector(AMREX_D_DECL(0.0, f_dir, f_dir));
        dir[10] = R*Set::Vector(AMREX_D_DECL(f_dir, 0.0, -f_dir));
        dir[11] = R*Set::Vector(AMREX_D_DECL(-f_dir, -f_dir, 0.0));

        for (int i = 0; i < 12; i++)
        {
            g_alpha[i] = a_g_alpha[i];
            m_exp[i] = a_exp[i];
            gam0dot[i] = a_gam0dot[i];
        }
    
        prev = Set::Matrix::Zero();
        curr = Set::Matrix::Zero();
    }

    Set::Scalar YieldSurface()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<12; i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = dir[i](j)*plane[i](k);
            yield += g_alpha[i]*R;
        }
        return yield.norm();
    }

    Set::Scalar PlasticEnergy()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<12; i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = dir[i](j)*plane[i](k);
            yield += g_alpha[i]*R;
        }
        return (yield.transpose()*curr).trace();
    }

    // Hardening law
    // \dot{g}^{\alpha} = \sum_{\beta=1}^N h_{\alpha\beta} |\dot{\gamma}^\beta|
    // h_{\alpha\beta} = q_{\alpha\beta} h(\bar{\gamma})
    // h(\bar{\gamma}) = h_s + (h_s - h_0) \sech^2 (\frac{h_0-h_s}{g_s-g_0} \bar{\gamma})
    // \bar{\gamma} = \int_0^t \sum_{\alpha=1}^N |\dot{\gamma}^\alpha| dt
    // q_{\alpha\beta} = 1 (if \alpha,\beta coplanar), q (else)

    
    void Advance(const Set::Matrix sigma, const Set::Matrix /*strain*/, const Set::Matrix3 /*epsgrad*/, Set::Scalar dt)
    {
        curr = prev;
        Set::Matrix R = Set::Matrix::Zero();

        std::array<Set::Scalar, 12> gamdot;

        for (int i = 0; i < 12; i++)
        {
            gamdot[i] = 0.0;
            Set::Scalar tau = dir[i].dot(sigma*plane[i]);
            Set::Scalar strength = g_alpha[i];
            if(tau/strength < 1.0) continue;
            gamdot[i] = gam0dot[i]*std::pow((tau/strength), m_exp[i]);
        }

        for (int i = 0; i < 12; i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = dir[i](j)*plane[i](k);

            curr += dt*gamdot[i]*R;
            gam_bar += dt*std::abs(gamdot[i]);
        }

        if (hardening)
        {
            Set::Scalar sech_gambar = 1.0/(std::cosh(((h0-hs)/(gs-g0))*gam_bar));
            Set::Scalar h_gambar = hs + (hs-h0)*sech_gambar*sech_gambar;
            for (int i = 0; i < 12; i++)
            {
                Set::Scalar gdotalpha = 0;
                for (int j = 0; j < 12; j++)
                {
                    Set::Scalar h_ab = ((plane[i] == plane[j]) ? 1.0 : q) * h_gambar;
                    gdotalpha += h_ab*std::abs(gamdot[j]);
                }
                g_alpha[i] += dt*gdotalpha;
            }
        }
        Cubic::SetF0(curr);
    }

    void SetPlasticStrains(Set::Matrix &a_state)
    {
        prev = a_state;
        curr = a_state;
    }

    Set::Matrix GetPlasticState()
    {
        return curr;
    }

    void SetHardening(bool flag) { hardening = flag; }

private: // this is to fix an issue with the intel compiler and constexpr
  constexpr static Set::Scalar INVSQ2 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105351942519376716382078636750692311545614851246241802792536860632206074854996791570661133296375279637789997525057639103028574;
  constexpr static Set::Scalar INVSQ3 = 0.57735026918962576450914878050195745564760175127012687601860232648397767230293334569371539558574952522520871380513556767665664836499965082627055183736479121617603107730076852735599160670036155830775501;
public:

    // For now FCC is hardcoded in here.
    std::array<Set::Vector, 12> plane;
    std::array<Set::Vector, 12> dir;
    std::array<Set::Scalar, 12> m_exp;
    std::array<Set::Scalar, 12> gam0dot;
    std::array<Set::Scalar, 12> g_alpha;

    bool hardening = false;
    Set::Scalar g0, gs, h0, hs, q, gam_bar;
    Set::Matrix curr, prev;

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.Define(0.0,0.0,0.0,0.0,0.0,0.0);
        ret.DefineSlipSystem(0.,0.,0.);
        ret.hardening = false;
        ret.g0 = 0.0;
        ret.gs = 0.0;
        ret.h0 = 0.0;
        ret.hs = 0.0;
        ret.q = 0.0;
        ret.gam_bar = 0.0;
        return ret;
    }

    static void Parse(CrystalPlastic & value, IO::ParmParse & pp)
    {
        Affine::Cubic::Parse(value,pp);
        
        amrex::Vector<Set::Scalar> input_g_alpha;
        std::array<Set::Scalar,12> f_g_alpha;
        
        pp.query("hardening.on", value.hardening);
        if (value.hardening)
        {
            pp.query("hardening.g0",value.g0);
            pp.query("hardening.gs",value.gs);
            pp.query("hardening.h0",value.h0);
            pp.query("hardening.hs",value.hs);
            pp.query("hardening.q",value.q);
            value.gam_bar = 0.0;
            for (int i = 0; i< 12; i++) f_g_alpha[i] = value.g0;
        }
        else
        {
            if(pp.contains("g_alpha"))
            {
                pp.queryarr("g_alpha", input_g_alpha);
                if(input_g_alpha.size() == 1)
                {
                    for (int i = 0; i< 12; i++) 
                        f_g_alpha[i] = input_g_alpha[0];
                }
                else if(input_g_alpha.size() <= 0 || input_g_alpha.size() != 12) Util::Abort(INFO, "Invalid number of values in g_alpha");
                else
                {
                    for (int i = 0; i< 12; i++) 
                        f_g_alpha[i] = input_g_alpha[i];
                }
            }
            else
                for (int i = 0; i< 12; i++) f_g_alpha[i] = 16.0;
        }

        amrex::Vector<Set::Scalar> input_exp;
        std::array<Set::Scalar, 12> f_g_exp;
        
        if(pp.contains("exp"))
        {
            pp.queryarr("exp",input_exp);
            if(input_exp.size() == 1)
            {
                for (int i = 0; i< 12; i++) 
                    f_g_exp[i] = input_exp[0];
            }
            else if(input_exp.size() <= 0 || input_exp.size() != 12) Util::Abort(INFO, "Invalid number of values in exp");
            else
            {
                for (int i = 0; i< 12; i++) 
                    f_g_exp[i] = input_exp[i];
            }
        }
        else
            for (int i = 0; i< 12; i++) f_g_exp[i] = 10.0;

        amrex::Vector<Set::Scalar> input_gam0dot;
        std::array<Set::Scalar,12> f_gam0dot;
        if(pp.contains("gam0dot"))
        {
            pp.queryarr("gam0dot",input_gam0dot);
            if(input_gam0dot.size() == 1)
            {
                for (int i = 0; i< 12; i++) 
                    f_gam0dot[i] = input_gam0dot[0];
            }
            else if(input_gam0dot.size() <= 0 || input_gam0dot.size() != 12) Util::Abort(INFO, "Invalid number of values in gam0dot");
            else
            {
                for (int i = 0; i< 12; i++) 
                    f_gam0dot[i] = input_gam0dot[i];
            }
        }
        else
            for (int i = 0; i< 12; i++) f_gam0dot[i] = 0.001;

        value.DefineSlipSystem(f_g_alpha,f_g_exp,f_gam0dot);
    }
    AMREX_FORCE_INLINE
	void operator += (const CrystalPlastic &rhs)
    {
        curr += rhs.curr;
        prev += rhs.prev;
        for (int i = 0; i < 12; i++)
        {
            m_exp[i] += rhs.m_exp[i];
            g_alpha[i] += rhs.g_alpha[i];
            gam0dot[i] += rhs.gam0dot[i];
        }
        
        g0 += rhs.g0; 
        gs += rhs.gs;
        h0 += rhs.h0;
        hs += rhs.hs;
        q += rhs.q;
        gam_bar += rhs.gam_bar;
        ddw += rhs.ddw;
        F0  += rhs.F0;
    }
    AMREX_FORCE_INLINE
    CrystalPlastic operator * (Set::Scalar alpha)
    {
        CrystalPlastic ret;
        ret.curr = curr*alpha;
        ret.prev = prev*alpha;

        ret.DefineSlipSystem(0.,0.,0.);

        for (int i = 0; i < 12; i++)
        {
            ret.m_exp[i] = alpha*m_exp[i];
            ret.g_alpha[i] = alpha*g_alpha[i];
            ret.gam0dot[i] = alpha*gam0dot[i];
        }

        ret.g0 = g0*alpha; 
        ret.gs = gs*alpha; 
        ret.h0 = h0*alpha; 
        ret.hs = hs*alpha;
        ret.q = q*alpha;
        ret.gam_bar = gam_bar*alpha;
        ret.hardening = hardening;
        ret.ddw = ddw*alpha;
        ret.F0  = alpha*F0;
    
        return ret;
    }
    friend CrystalPlastic operator * (const Set::Scalar alpha, const CrystalPlastic b);
    friend CrystalPlastic operator + (const CrystalPlastic a, const CrystalPlastic b);
    friend CrystalPlastic operator - (const CrystalPlastic a, const CrystalPlastic b);
};
AMREX_FORCE_INLINE
CrystalPlastic operator * (const Set::Scalar alpha, const CrystalPlastic b)
{
    CrystalPlastic ret;
    ret.curr = b.curr*alpha;
    ret.prev = b.prev*alpha;
    ret.DefineSlipSystem(0.,0.,0.);

    for (int i = 0; i < 12; i++)
    {
        ret.m_exp[i] = alpha*b.m_exp[i];
        ret.g_alpha[i] = alpha*b.g_alpha[i];
        ret.gam0dot[i] = alpha*b.gam0dot[i];
    }
    ret.g0 = b.g0*alpha; 
    ret.gs = b.gs*alpha; 
    ret.h0 = b.h0*alpha; 
    ret.hs = b.hs*alpha;
    ret.q = b.q*alpha;
    ret.gam_bar = b.gam_bar*alpha;
    ret.hardening = b.hardening;

    ret.ddw = b.ddw*alpha;
    ret.F0  = alpha*b.F0;
    
    return ret;
}

AMREX_FORCE_INLINE
CrystalPlastic operator + (const CrystalPlastic a, const CrystalPlastic b)
{
    CrystalPlastic ret;
    ret.curr = a.curr + b.curr;
    ret.prev = a.prev + b.prev;

    ret.DefineSlipSystem(0.,0.,0.);

    for (int i = 0; i < 12; i++)
    {
        ret.m_exp[i] = a.m_exp[i] + b.m_exp[i];
        ret.g_alpha[i] = a.g_alpha[i] + b.g_alpha[i];
        ret.gam0dot[i] = a.gam0dot[i] + b.gam0dot[i] ;
    }

    ret.g0 = a.g0 + b.g0; 
    ret.gs = a.gs + b.gs; 
    ret.h0 = a.h0 + b.h0; 
    ret.hs = a.hs + b.hs;
    ret.q = a.q + b.q;
    ret.gam_bar = a.gam_bar + b.gam_bar;
    ret.hardening = a.hardening && b.hardening;

    ret.ddw = a.ddw + b.ddw;
    ret.F0  = a.F0 + b.F0;
    
    return ret;
}

AMREX_FORCE_INLINE
CrystalPlastic operator - (const CrystalPlastic a, const CrystalPlastic b)
{
    CrystalPlastic ret;
    ret.curr = a.curr - b.curr;
    ret.prev = a.prev - b.prev;
    ret.DefineSlipSystem(0.,0.,0.);

    for (int i = 0; i < 12; i++)
    {
        ret.m_exp[i] = a.m_exp[i] - b.m_exp[i];
        ret.g_alpha[i] = a.g_alpha[i] - b.g_alpha[i];
        ret.gam0dot[i] = a.gam0dot[i] - b.gam0dot[i] ;
    }

    ret.g0 = a.g0 - b.g0; 
    ret.gs = a.gs - b.gs; 
    ret.h0 = a.h0 - b.h0; 
    ret.hs = a.hs - b.hs;
    ret.q = a.q - b.q;
    ret.gam_bar = a.gam_bar - b.gam_bar;
    ret.hardening = a.hardening && b.hardening;

    ret.ddw = a.ddw - b.ddw;
    ret.F0  = a.F0 - b.F0;
    
    return ret;
}

}
}
}

#endif
