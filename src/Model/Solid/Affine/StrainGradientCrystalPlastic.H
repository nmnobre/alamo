#ifndef MODEL_SOLID_PLASTIC_STRAINGRADIENTCRYSTALPLASTIC_H_
#define MODEL_SOLID_PLASTIC_STRAINGRADIENTCRYSTALPLASTIC_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Affine/CrystalPlastic.H"

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
class StrainGradientCrystalPlastic : public CrystalPlastic
{
public:
    StrainGradientCrystalPlastic(){};

    void DefineSlipGradientSystem()
    {
        DefineSlipGradientSystem(0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
    }

    void DefineSlipGradientSystem(Set::Scalar a_Ke, Set::Scalar a_b, Set::Scalar a_de, Set::Scalar a_ds, Set::Scalar a_tau0, Set::Scalar a_mu, 
                                    Set::Scalar a_c0, Set::Scalar a_w1, Set::Scalar a_w2, Set::Scalar a_h, Set::Scalar a_q1, Set::Scalar a_q2,
                                    Set::Scalar a_rho0, Set::Scalar a_F0, Set::Scalar a_T)
    {
        m_Ke = a_Ke;
        m_b = a_b;
        m_de = a_de;
        m_ds = a_ds;
        m_tau0 = a_tau0;
        m_mu = a_mu;
        m_c0 = a_c0;
        m_w1 = a_w1;
        m_w2 = a_w2;
        m_h = a_h;
        m_q1 = a_q1;
        m_q2 = a_q2;
        m_rho0 = a_rho0;
        m_F0 = a_F0;
        m_T = a_T;

        for (int i=0; i<12; i++)
        {
            rho_se[i] = a_rho0/24.0;
            rho_ss[i] = a_rho0/24.0;
            rho_gen[i] = 0.0;
            rho_get[i] = 0.0;
            rho_gsm[i] = 0.0;
            gamalpha[i] = 0.0; 
            gamdotalpha[i] = 0.0;
        }
    }

    void ComputePlasticStrainRate(  Set::Matrix sigma, 
                                    std::array<Set::Vector, 12> &gamdotgrad, 
                                    std::array<Set::Vector, 12> &gamgrad,
                                    Set::Scalar dt)
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        Set::Scalar total_rho = 0.0;
        for (int i = 0; i<12; i++) 
            total_rho += rho_se[i] + rho_ss[i] + std::abs(rho_gen[i]) + std::abs(rho_get[i]) + std::abs(rho_gsm[i]);
        
        // Let's first compute all dislocation densities
        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                int alpha = 3*i + j;
                Set::Vector t_alpha = R_crystal*dir[i][j].cross(R_crystal*plane[i]);
                Set::Scalar rho_se_dot = ((m_Ke/m_b)*std::sqrt(total_rho) - m_de*rho_se[alpha])*std::abs(gamdotalpha[alpha]);
                Set::Scalar rho_ss_dot = ((1. - (Set::Constant::Pi*m_ds*m_ds*rho_ss[alpha]/4.0))*(2.0*m_Ke/m_b)*std::sqrt(total_rho) - m_ds*rho_ss[alpha])*std::abs(gamdotalpha[alpha]);
                Set::Vector curlTerm = Set::Vector::Zero();
                for (int bi=0; bi<4; bi++)
                {
                    for (int bj=0; bj<3; bj++)
                    {
                        int beta = 3*bi + bj;
                        for (int p=0; p<3; p++)
                            for (int q=0; q<3; q++)
                                for (int r=0; r<3; r++)
                                    curlTerm(r) = curlTerm(r) + (gamdotgrad[alpha](q)*gamalpha[beta] + gamdotalpha[alpha]*gamgrad[beta](q))*((R_crystal*plane[bi])(p))
                                                    *((R_crystal*plane[i]).dot(R_crystal*plane[bi]))
                                                    *Set::Constant::PermutationTensor(p,q,r);
                    }
                }
                
                Set::Scalar rho_gsm_dot = (1.0/m_b)*(R_crystal*dir[i][j]).dot(curlTerm);
                Set::Scalar rho_get_dot = (1.0/m_b)*(t_alpha).dot(curlTerm);
                Set::Scalar rho_gen_dot = (1.0/m_b)*(R_crystal*plane[i]).dot(curlTerm);

                rho_se[alpha] += rho_se_dot*dt;
                rho_ss[alpha] += rho_ss_dot*dt;
                rho_gsm[alpha] += rho_gsm_dot*dt;
                rho_get[alpha] += rho_get_dot*dt;
                rho_gen[alpha] += rho_gen_dot*dt;
            }
        }

        total_rho = 0.0;
        for (int k = 0; k<12; k++) total_rho += rho_se[k] + rho_ss[k] + std::abs(rho_gen[k]) + std::abs(rho_get[k]) + std::abs(rho_gsm[k]);

        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                int alpha = 3*i+j;
                Set::Scalar tau = (R_crystal*dir[i][j]).dot(sigma*(R_crystal*plane[i]));
                Set::Scalar tau_prime = m_tau0 + m_mu*m_b*m_c0*std::sqrt(total_rho);

                Set::Scalar anisotropic_rho = 0.0;
                for (int k=0; k<12; k++) anisotropic_rho += (m_w1 + (k == alpha ? 1.0 : 0.0)*(1.0-m_w2))*m_h*(rho_se[k] + rho_ss[k] + std::abs(rho_gen[k]) + std::abs(rho_get[k]) + std::abs(rho_gsm[k]));
                Set::Scalar S_alpha = m_mu*m_b*std::sqrt(anisotropic_rho);
                
                Set::Scalar temp1 = (std::abs(tau) - S_alpha) / tau_prime;
                temp1 = temp1 > 0.0 ? std::pow(temp1,m_q1) : 0.0;
                Set::Scalar temp2 = 1. - temp1;
                temp2 = temp2 > 0.0 ? std::pow(temp2,m_q2) : 0.0;
                gamdotalpha[alpha] = gam0dot[alpha]*std::exp(-m_F0*temp2/(m_T))*(tau>0.0 ? 1.0 : (tau < 0.0 ? -1.0 : 0.0));
            }
        }
        
    }

    void Advance(const Set::Matrix /*sigma*/, const Set::Matrix /*strain*/, const Set::Matrix3 /*epsgrad*/, Set::Scalar dt)
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        //std::array<Set::Scalar, 12> gamdot;
        
        for (int i=0; i< 4; i++)
        {
            for (int j=0; j<3; j++)
            {
                F0 += dt * gamdotalpha[3*i + j] * ((R_crystal*dir[i][j]) * (R_crystal*plane[i]).transpose());
                gamalpha[3*i + j] += dt*gamdotalpha[3*i + j];
            }
        }

    }

public:

    std::array<Set::Scalar,12> gamalpha;
    std::array<Set::Scalar,12> gamdotalpha;
    std::array<Set::Scalar,12> rho_gsm;
    std::array<Set::Scalar,12> rho_get;
    std::array<Set::Scalar,12> rho_gen;
    std::array<Set::Scalar,12> rho_se;
    std::array<Set::Scalar,12> rho_ss;

    Set::Scalar m_Ke;
    Set::Scalar m_b;
    Set::Scalar m_de;
    Set::Scalar m_ds;
    Set::Scalar m_tau0;
    Set::Scalar m_mu;
    Set::Scalar m_c0;
    Set::Scalar m_rho0;
    Set::Scalar m_w1, m_w2, m_h;
    Set::Scalar m_q1, m_q2;
    Set::Scalar m_F0, m_T;

    
public:
    static StrainGradientCrystalPlastic Zero()
    {
        StrainGradientCrystalPlastic ret;
        ret.q.coeffs() *= 0;
        ret.C11 = 0.0; ret.C12 = 0.0; ret.C44 = 0.0;
        ret.F0 = Set::Matrix::Zero();
        ret.DefineSlipSystem(0.,0.,0.);
        ret.DefineSlipGradientSystem();
        ret.hardening = false;
        ret.g0 = 0.0;
        ret.gs = 0.0;
        ret.h0 = 0.0;
        ret.hs = 0.0;
        ret.qh = 0.0;
        ret.gam_bar = 0.0;
        return ret;
    }

    static void Parse(StrainGradientCrystalPlastic & value, IO::ParmParse & pp)
    {
        CrystalPlastic::Parse(value,pp);

        Set::Scalar a_Ke = 35.0e3;
        Set::Scalar a_b = 0.248; //nm
        Set::Scalar a_de = 5.e0;
        Set::Scalar a_ds = 25.0e0;
        Set::Scalar a_tau0 = 3.40e-1; //GPa
        Set::Scalar a_mu = value.C44;
        Set::Scalar a_c0 = 8.7e-5;
        Set::Scalar a_w1 = 1.0, a_w2 = 1.0, a_h = 0.084*0.084;
        Set::Scalar a_q1 = 0.5, a_q2 = 1.25e2;
        Set::Scalar a_rho0 = 6.e12;
        Set::Scalar a_F0 = 3.114e4; // divided by kB
        Set::Scalar a_T = 298.;

        pp.query("Ke", a_Ke);
        pp.query("b", a_b);
        pp.query("de", a_de);
        pp.query("ds", a_ds);
        pp.query("tau0", a_tau0);
        pp.query("mu", a_mu);
        pp.query("c0", a_c0);
        pp.query("w1", a_w1);
        pp.query("w2", a_w2);
        pp.query("h", a_h);
        pp.query("q1", a_q1);
        pp.query("q2", a_q2);
        pp.query("rho0", a_rho0);
        pp.query("m_F0", a_F0);
        pp.query("T", a_T);

        value.DefineSlipGradientSystem(a_Ke, a_b, a_de, a_ds, a_tau0, a_mu, 
                                    a_c0, a_w1, a_w2, a_h, a_q1, a_q2, a_rho0, a_F0, a_T);
    }

    #define OP_CLASS StrainGradientCrystalPlastic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) X(g0) X(gs) X(h0) X(hs) X(qh) X(gam_bar) X(m_Ke) X(m_b) X(m_de) X(m_ds) X(m_tau0) X(m_mu) X(m_c0) X(m_w1) X(m_w2) X(m_h) X(m_q1) X(m_q2) X(m_rho0) X(m_F0) X(m_T)
    #define OP_LOOP_NUM 12
    #define OP_LOOP_VARS X(m_exp) X(g_alpha) X(gam0dot) X(gamalpha) X(gamdotalpha) X(rho_se) X(rho_ss) X(rho_gen) X(rho_gsm) X(rho_get)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"


}
}
}

#endif
