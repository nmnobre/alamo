#ifndef MODEL_SOLID_PLASTIC_CRYSTALPLASTIC_H_
#define MODEL_SOLID_PLASTIC_CRYSTALPLASTIC_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Cubic.H"

#if AMREX_SPACEDIM==3
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
        for (int i = 0; i < 12; i++)
        {
            m_g_alpha[i] = a_g_alpha[i];
            // m_g_alpha0[i] = a_g_alpha[i];
            m_exp[i] = a_exp[i];
            m_gam0dot[i] = a_gam0dot[i];
        }
    
        F0 = Set::Matrix::Zero();
    }

    Set::Scalar YieldSurface()
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        Set::Matrix_3D yield = Set::Matrix_3D::Zero(), R = Set::Matrix_3D::Zero();
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                for (int k=0; k<3; k++)
                {
                    for (int l=0; l<3; l++) 
                        R(k,l) = (R_crystal*plane[i])(k)*(R_crystal*dir[i][j])(l);
                }
                yield += m_g_alpha[3*i + j]*R;
            }
        }
        return yield.norm();
    }

    Set::Scalar PlasticEnergy()
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        Set::Matrix_3D yield = Set::Matrix_3D::Zero(), R = Set::Matrix_3D::Zero();
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                for (int k=0; k<3; k++)
                {
                    for (int l=0; l<3; l++)
                        R(k,l) = (R_crystal*plane[i])(k)*(R_crystal*dir[i][j])(l);
                }
                yield += m_g_alpha[3*i + j]*R;
            }
        }
        return (yield.transpose()*F0).trace();
    }

    void DegradeYieldSurface(Set::Scalar eta)
    {
        if (std::isnan(eta) || std::isinf(eta)) Util::Abort(INFO, "Nans or Infs in eta. ", eta);
        if (eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta. ", eta);
        // for (int i=0; i<12; i++)
        // {
        //     m_g_alpha[i] = m_g_alpha0[i]/eta;
        //     if (hardening)
        //     {
        //         m_g0 = m_g00/eta; m_gs = m_gs0/eta;
        //         m_h0 = m_h00/eta; m_hs = m_hs0/eta;
        //     }
        // }
    }

    Set::Scalar OriginalYieldSurface()
    {
        YieldSurface();
        // Set::Matrix_3D R_crystal;
        // R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        // Set::Matrix_3D yield = Set::Matrix_3D::Zero(), R = Set::Matrix_3D::Zero();
        // for (int i=0; i<4; i++)
        // {
        //     for (int j=0; j<3; j++)
        //     {
        //         for (int k=0; k<3; k++)
        //         {
        //             for (int l=0; l<3; l++) 
        //                 R(k,l) = (R_crystal*plane[i])(k)*(R_crystal*dir[i][j])(l);
        //         }
        //         yield += m_g_alpha0[3*i + j]*R;
        //     }
        // }
        // return yield.norm();
    }

    Set::Scalar OriginalPlasticEnergy()
    {
        PlasticEnergy();
        // Set::Matrix_3D R_crystal;
        // R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        // Set::Matrix_3D yield = Set::Matrix_3D::Zero(), R = Set::Matrix_3D::Zero();
        // for (int i=0; i<4; i++)
        // {
        //     for (int j=0; j<3; j++)
        //     {
        //         for (int k=0; k<3; k++)
        //         {
        //             for (int l=0; l<3; l++)
        //                 R(k,l) = (R_crystal*plane[i])(k)*(R_crystal*dir[i][j])(l);
        //         }
        //         yield += m_g_alpha0[3*i + j]*R;
        //     }
        // }
        // return (yield.transpose()*F0).trace();
    }

    // Hardening law
    // \dot{g}^{\alpha} = \sum_{\beta=1}^N h_{\alpha\beta} |\dot{\gamma}^\beta|
    // h_{\alpha\beta} = qh_{\alpha\beta} h(\bar{\gamma})
    // h(\bar{\gamma}) = h_s + (h_s - h_0) \sech^2 (\frac{h_0-h_s}{g_s-g_0} \bar{\gamma})
    // \bar{\gamma} = \int_0^t \sum_{\alpha=1}^N |\dot{\gamma}^\alpha| dt
    // qh_{\alpha\beta} = 1 (if \alpha,\beta coplanar), qh (else)
    void Advance(Set::Scalar dt, Set::Matrix /*eps*/, Set::Matrix sig) override
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        std::array<Set::Scalar, 12> gamdot;
        
        for (int i=0; i< 4; i++)
        {
            Set::Scalar tau = 0.0;
            for (int j=0; j<3; j++)
            {
                gamdot[3*i + j] = 0.0;
                
                Set::Scalar strength = m_g_alpha[3*i + j];
                tau = (R_crystal*dir[i][j]).dot(sig*(R_crystal*plane[i]));
                if(fabs(tau) < fabs(strength)) continue;

                Set::Scalar incr_abs = m_gam0dot[3*i + j]*std::pow(fabs(tau/strength), m_exp[3*i + j]);
                gamdot[3*i + j] = tau>0 ? +incr_abs : -incr_abs;

                F0 += dt * gamdot[3*i + j] * ((R_crystal*dir[i][j]) * (R_crystal*plane[i]).transpose());
                m_gam_bar += dt*std::abs(gamdot[3*i + j]);

                if(std::isnan(F0.norm()) || std::isnan(m_gam_bar))
                {
                    Util::Message(INFO, "tau = ", tau, ". gamdot [", 3*i+j, "] = ", gamdot[3*i + j], ". gambar = ", m_gam_bar, ". incr_abs = ", incr_abs);
                    Util::Message(INFO, "gam0dot[",3*i+j,"] = ",m_gam0dot[3*i+j], ". strength = ", strength, ". m_exp[",3*i+j,"] = ",m_exp[3*i+j]);
                    Util::Abort(INFO, "Nans detected");
                }
            }
        }

        if (hardening)
        {
            Set::Scalar sech_gambar = 1.0/(std::cosh(((m_h0-m_hs)/(m_gs-m_g0))*m_gam_bar));
            Set::Scalar h_gambar = m_hs + (m_hs-m_h0)*sech_gambar*sech_gambar;

            for (int i = 0; i < 12; i++)
            {
                Set::Scalar gdotalpha = 0;
                for (int j = 0; j < 12; j++)
                {
                    Set::Scalar h_ab = ((plane[i/3] == plane[j/3]) ? 1.0 : m_qh) * h_gambar;
                    gdotalpha += h_ab*std::abs(gamdot[j]);
                }
                m_g_alpha[i] += dt*gdotalpha;
            }
        }
    }

    void SetHardening(bool flag) { hardening = flag; }

public:

    // For now FCC is hardcoded in here.
    inline static const std::array<Set::Vector_3D, 4> plane = 
    {
        Set::Vector_3D(1,1,1).normalized(),
        Set::Vector_3D(-1,-1,1).normalized(),
        Set::Vector_3D(-1,1,1).normalized(),
        Set::Vector_3D(1,-1,1).normalized()
    };
    inline static const std::array<std::array<Set::Vector_3D,3>,4> dir = 
    {
        Set::Vector_3D(0,-1,1).normalized(), Set::Vector_3D(1,0,-1).normalized(), Set::Vector_3D(-1,1,0).normalized(),
        Set::Vector_3D(0,1,1).normalized(), Set::Vector_3D(-1,0,-1).normalized(), Set::Vector_3D(1,-1,0).normalized(),
        Set::Vector_3D(0,-1,1).normalized(), Set::Vector_3D(-1,0,-1).normalized(), Set::Vector_3D(1,1,0).normalized(),
        Set::Vector_3D(0,1,1).normalized(), Set::Vector_3D(1,0,-1).normalized(), Set::Vector_3D(-1,-1,0).normalized()
    };
    //std::array<Set::Vector_3D, 12> dir;
    std::array<Set::Scalar, 12> m_exp;
    std::array<Set::Scalar, 12> m_gam0dot;
    std::array<Set::Scalar, 12> m_g_alpha;
    // std::array<Set::Scalar, 12> m_g_alpha0;

    Set::Scalar m_gam_bar;
    Set::Scalar m_g0, m_gs, m_h0, m_hs, m_qh;
    // Set::Scalar m_g00, m_gs0, m_h00, m_hs0;

    static bool hardening;

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;

        ret.q.coeffs() *= 0;
        ret.C11 = 0.0; ret.C12 = 0.0; ret.C44 = 0.0;
        // ret.C110 = 0.0; ret.C120 = 0.0; ret.C440 = 0.0;
        ret.F0 = Set::Matrix::Zero();

        ret.DefineSlipSystem(0.,0.,0.);
        ret.hardening = false;
        
        ret.m_g0 = 0.0; //ret.m_g00 = 0.0;
        ret.m_gs = 0.0; //ret.m_gs0 = 0.0;
        ret.m_h0 = 0.0; //ret.m_h00 = 0.0;
        ret.m_hs = 0.0; //ret.m_hs0 = 0.0;
        ret.m_qh = 0.0;
        ret.m_gam_bar = 0.0;
        return ret;
    }

    static void Parse(CrystalPlastic & value, IO::ParmParse & pp)
    {
        Cubic::Parse(value,pp);
        
        amrex::Vector<Set::Scalar> input_g_alpha;
        std::array<Set::Scalar,12> f_g_alpha;
        
        pp.query("hardening.on", value.hardening);
        if (value.hardening)
        {
            pp.query("hardening.g0",value.m_g0);
            pp.query("hardening.gs",value.m_gs);
            pp.query("hardening.h0",value.m_h0);
            pp.query("hardening.hs",value.m_hs);
            pp.query("hardening.qh",value.m_qh);
            
            // value.m_g00 = value.m_g0;
            // value.m_gs0 = value.m_gs;
            // value.m_h00 = value.m_h0;
            // value.m_hs0 = value.m_hs;

            value.m_gam_bar = 0.0;
            for (int i = 0; i< 12; i++) f_g_alpha[i] = value.m_g0;
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
                for (int i = 0; i< 12; i++) f_g_alpha[i] = 0.016;
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

    #define OP_CLASS CrystalPlastic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) X(m_g0) X(m_gs) X(m_h0) X(m_hs) X(m_qh) X(m_gam_bar) // X(C110) X(C120) X(C440)
    #define OP_LOOP_NUM 12
    #define OP_LOOP_VARS X(m_exp) X(m_g_alpha) X(m_gam0dot) // X(m_g_alpha0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

bool CrystalPlastic::hardening = false;

}
}
}

template<>
inline int Set::Field<Model::Solid::Affine::CrystalPlastic>::NComp() const 
{
    return AMREX_SPACEDIM*AMREX_SPACEDIM;
}

template<>
inline std::string Set::Field<Model::Solid::Affine::CrystalPlastic>::Name(int i) const 
{
    if (i==0) return name + ".epsp_xx";
    if (i==1) return name + ".epsp_xy";
    if (i==2) return name + AMREX_D_PICK("",".epsp_yx",".epsp_xz");
    if (i==3) return name + AMREX_D_PICK("",".epsp_yy",".epsp_yx");
    if (i==4) return name + ".epsp_yy";
    if (i==5) return name + ".epsp_yz";
    if (i==6) return name + ".epsp_zx";
    if (i==7) return name + ".epsp_zy";
    if (i==8) return name + ".epsp_zz";
    return name;
}

template<>
inline void Set::Field<Model::Solid::Affine::CrystalPlastic>::Copy(int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok())
        {
            amrex::Array4<const Model::Solid::Affine::CrystalPlastic> const & src = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const & dst = a_dst.array(mfi);
            for (int n = 0; n < AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp + n) = src(i,j,k).F0(n/AMREX_SPACEDIM,n%AMREX_SPACEDIM);
                });
            }
        }
    }    
}

#endif
#endif
