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
        Set::Matrix_3D R;
        R = q.normalized().toRotationMatrix();

        Set::Scalar f_dir = Set::Constant::INVSQ2;
        Set::Scalar f_nor = Set::Constant::INVSQ3;

        for (int i = 0; i < 12; i++)
        {
            g_alpha[i] = a_g_alpha[i];
            m_exp[i] = a_exp[i];
            gam0dot[i] = a_gam0dot[i];
        }
    
        F0 = Set::Matrix::Zero();
    }

    Set::Scalar YieldSurface()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                for (int k=0; k<3; k++)
                {
                    for (int l=0; l<3; l++)
                        R(k,l) = plane[i](k)*dir[i][j](l);
                }
                yield += g_alpha[3*i + j]*R;
            }
        }
        // for(int i=0; i<12; i++)
        //{
            //for (int j = 0; j < AMREX_SPACEDIM; j++) // TODO
            //    for (int k = 0; k < AMREX_SPACEDIM; k++) // TODO
            //        R(j,k) = dir[i](j)*plane[i](k); // TODO
            //yield += g_alpha[i]*R;
        //}
        return yield.norm();
    }

    Set::Scalar PlasticEnergy()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                for (int k=0; k<3; k++)
                {
                    for (int l=0; l<3; l++)
                        R(k,l) = plane[i](k)*dir[i][j](l);
                }
                yield += g_alpha[3*i + j]*R;
            }
        }
        //for(int i=0; i<12; i++)
        //{
            //for (int j = 0; j < AMREX_SPACEDIM; j++) // TODO
            //    for (int k = 0; k < AMREX_SPACEDIM; k++) // TODO
            //        R(j,k) = dir[i](j)*plane[i](k); // TODO
            //yield += g_alpha[i]*R;
        //}
        return (yield.transpose()*F0).trace();
    }

    // Hardening law
    // \dot{g}^{\alpha} = \sum_{\beta=1}^N h_{\alpha\beta} |\dot{\gamma}^\beta|
    // h_{\alpha\beta} = qh_{\alpha\beta} h(\bar{\gamma})
    // h(\bar{\gamma}) = h_s + (h_s - h_0) \sech^2 (\frac{h_0-h_s}{g_s-g_0} \bar{\gamma})
    // \bar{\gamma} = \int_0^t \sum_{\alpha=1}^N |\dot{\gamma}^\alpha| dt
    // qh_{\alpha\beta} = 1 (if \alpha,\beta coplanar), qh (else)
    void Advance(const Set::Matrix sigma, const Set::Matrix /*strain*/, const Set::Matrix3 /*epsgrad*/, Set::Scalar dt)
    {
        std::array<Set::Scalar, 12> gamdot;
        //std::array<Set::Scalar, 12> temp_g_alpha;

        for (int i=0; i< 4; i++)
        {
            Set::Scalar tau = 0.0;
            for (int j=0; j<3; j++)
            {
                gamdot[3*i + j] = 0.0;
                
                Set::Scalar strength = g_alpha[3*i + 4];
                tau = dir[i][j].dot(sigma*plane[i]);
                if(fabs(tau) < fabs(strength)) continue;

                Set::Scalar incr_abs = gam0dot[3*i + j]*std::pow(fabs(tau/strength), m_exp[3*i + j]);
                gamdot[3*i + j] = tau>0 ? +incr_abs : -incr_abs;

                F0 += dt * gamdot[3*i + j] * (dir[i][j] * plane[i].transpose());
                gam_bar += dt*std::abs(gamdot[3*i + j]);
            }
        }

        // for (int i = 0; i < 12; i++)
        // {
        //     gamdot[i] = 0.0;
        //     Set::Scalar tau = 0.0;
        //     Util::Message(INFO,"sigma=",sigma);
        //     Util::Message(INFO,"plane=",plane[i]);

        //     Set::Scalar strength = g_alpha[i];
        //     temp_g_alpha[i] = g_alpha[i];

        //     // if (fabs(strength) < 1.e-6) Util::Abort(INFO, "g_alpha[",i,"] = ", g_alpha[i]);

        //     Util::Message(INFO,tau," ",strength);
        //     if(fabs(tau) < fabs(strength)) continue;

        //     Set::Scalar incr_abs = gam0dot[i]*std::pow(fabs(tau/strength), m_exp[i]);
        //     gamdot[i] = tau>0 ? +incr_abs : -incr_abs;

        //     // if(gamdot[i] > 1.e10) 
        //     // {
        //     //     Util::Message(INFO, "strength = ", strength);
        //     //     Util::Message(INFO, "tau = ", tau);
        //     //     Util::Message(INFO, "gam0dot[",i,"] = ",gam0dot[i]);
        //     //     Util::Message(INFO, "m_exp[",i,"] = ", m_exp[i]);
        //     //     Util::Abort(INFO, "incr_abs = ", incr_abs);
        //     // }
        //     // Util::Message(INFO, "tau = ", tau, ". strength = ", strength, ". gamdot[",i,"] = ",gamdot[i]);
        // }

        // for (int i = 0; i < 12; i++)
        // {
        //     //F0 += dt * gamdot[i] * (dir[i] * plane[i].transpose()); // TODO
        //     gam_bar += dt*std::abs(gamdot[i]);
        // }
        // Util::Message(INFO,"gam_bar = ", gam_bar);

        if (hardening)
        {
            Set::Scalar sech_gambar = 1.0/(std::cosh(((h0-hs)/(gs-g0))*gam_bar));
            Set::Scalar h_gambar = hs + (hs-h0)*sech_gambar*sech_gambar;

            for (int i = 0; i < 12; i++)
            {
                Set::Scalar gdotalpha = 0;
                for (int j = 0; j < 12; j++)
                {
                    Set::Scalar h_ab = ((plane[i/3] == plane[j/3]) ? 1.0 : qh) * h_gambar;
                    gdotalpha += h_ab*std::abs(gamdot[j]);
                }
                g_alpha[i] += dt*gdotalpha;
            }
        }

        // if (F0.norm()>1.e6) //debugging
        // {
        //     Set::Scalar tau = 0.0;
        //     for (int i = 0; i< 12; i++)
        //     {
        //         //for (int j = 0; j < 12;  j++) tau += dir[i].dot(sigma*plane[j]); // TODO
        //         Util::Message(INFO, "g_alpha[",i,"] = ", g_alpha[i]);
        //         Util::Message(INFO, "temp_g_alpha[",i,"] = ", temp_g_alpha[i]);
        //         Util::Message(INFO, "gamdot[",i,"] = ", gamdot[i]);
        //         Util::Message(INFO, "gam0dot[",i,"] = ", gam0dot[i]);
        //         Util::Message(INFO, "m_exp[",i,"] = ", m_exp[i]);
        //     }
        //     Util::Message(INFO, "tau = ", tau);
        //     Util::Message(INFO, "gambar = ", gam_bar);
        //     Util::Abort(INFO);
        // }
    }

    void SetPlasticStrains(Set::Matrix &a_state)
    {
        F0 = a_state;
    }

    Set::Matrix GetPlasticState()
    {
        return F0;
    }

    void SetHardening(bool flag) { hardening = flag; }

public:

    // For now FCC is hardcoded in here.
    inline static const std::array<Set::Vector_3D, 4> plane = 
    {
        Set::Vector_3D(1,1,1).normalized(),
        Set::Vector_3D(1,1,1).normalized(),
        Set::Vector_3D(1,1,1).normalized(),
        Set::Vector_3D(1,1,1).normalized()
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
    std::array<Set::Scalar, 12> gam0dot;
    std::array<Set::Scalar, 12> g_alpha;

    Set::Scalar g0, gs, h0, hs, qh, gam_bar;

    static bool hardening;

public:
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        ret.q.coeffs() *= 0;
        ret.C11 = 0.0; ret.C12 = 0.0; ret.C44 = 0.0;
        ret.F0 = Set::Matrix::Zero();
        ret.DefineSlipSystem(0.,0.,0.);
        ret.hardening = false;
        ret.g0 = 0.0;
        ret.gs = 0.0;
        ret.h0 = 0.0;
        ret.hs = 0.0;
        ret.qh = 0.0;
        ret.gam_bar = 0.0;
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
            pp.query("hardening.g0",value.g0);
            pp.query("hardening.gs",value.gs);
            pp.query("hardening.h0",value.h0);
            pp.query("hardening.hs",value.hs);
            pp.query("hardening.qh",value.qh);
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
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) X(g0) X(gs) X(h0) X(hs) X(qh) X(gam_bar) 
    #define OP_LOOP_NUM 12
    #define OP_LOOP_VARS X(m_exp) X(g_alpha) X(gam0dot)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

bool CrystalPlastic::hardening = false;

//static constexpr std::array<Set::Vector_3D, 4> CrystalPlastic::plane = 
//{
//    Set::Vector_3D(1,1,1).normalized(),
//    Set::Vector_3D(1,1,1).normalized(),
//    Set::Vector_3D(1,1,1).normalized(),
//    Set::Vector_3D(1,1,1).normalized()
//}



}
}
}

#endif
