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

    void DefineSlipGradientSystem(Set::Scalar a_Le, Set::Scalar a_Ld, Set::Scalar a_visc=0.01)
    {
        std::array<Set::Scalar,12> temp_Le, temp_Ld;
        for (int i=0; i<12; i++)
        {
            temp_Le[i] = a_Le;
            temp_Ld[i] = a_Ld;
        }
        DefineSlipGradientSystem(temp_Le, temp_Ld, a_visc);
    }

    void DefineSlipGradientSystem(std::array<Set::Scalar,12> &a_Le, std::array<Set::Scalar,12> &a_Ld, Set::Scalar a_visc=0.01)
    {
        visc = a_visc;
        for (int i=0; i<12; i++)
        {
            Le[i] = a_Le[i];
            Ld[i] = a_Ld[i];
            gamalpha[i] = 0.0; 
            gamdotalpha[i] = 0.0;
        }
    }

    void ComputePlasticStrainRate(  const Set::Matrix sigma, 
                                    const std::array<Set::Vector, 12> &gamdotgrad, 
                                    const std::array<Set::Vector, 12> &gamgrad,
                                    const std::array<Set::Matrix, 12> &gamdotgradgrad,
                                    const std::array<Set::Matrix, 12> &gamgradgrad)
    {
        Set::Scalar tol = 1.e-10;
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        for (int i=0; i<4; i++)
        {
            for (int j=0; j<3; j++)
            {
                Set::Scalar tau = (R_crystal*dir[i][j]).dot(sigma*(R_crystal*plane[i]));

                Set::Scalar gamedot = std::sqrt(gamdotalpha[3*i+j]*gamdotalpha[3*i+j] + 
                                                (Ld[3*i+j]*gamdotgrad[3*i+j].dot(R_crystal*dir[i][j]))*
                                                (Ld[3*i+j]*gamdotgrad[3*i+j].dot(R_crystal*dir[i][j])));
                
                Set::Vector gamedotgrad = Set::Vector::Zero();
                if (std::abs(gamedot) > tol)
                {
                    gamedotgrad =   (1.0/gamedot)*
                                    ((gamdotalpha[3*i+j]*gamdotgrad[3*i+j]) +
                                    (Ld[3*i+j]*Ld[3*i+j]*(gamdotgrad[3*i+j].dot(R_crystal*dir[i][j]))*(gamdotgradgrad[3*i+j]*R_crystal*dir[i][j])));
                }

                Set::Scalar taue = g_alpha[3*i+j]*std::pow((gamedot/gam0dot[3*i+j]),m_exp[3*i+j]);
                Set::Vector tauegrad = g_alpha[3*i+j]*m_exp[3*i+j]*gamedotgrad*std::pow(gamedot,m_exp[3*i+j]-1)*std::pow(gam0dot[3*i+j],m_exp[3*i+j]);

                Set::Scalar qbar = taue*gamdotalpha[3*i+j]/gamedot;
                //Set::Scalar ziebar = taue*(Ld[3*i+j]*Ld[3*i+j]*gamdotgrad[3*i+j].dot(R_crystal*dir[i][j]))/gamedot;
                //Set::Scalar zietilde = C44*Le[3*i+j]*Le[3*i+j]*gamgrad[3*i+j].dot(R_crystal*dir[i][j]);

                Set::Scalar grad_zietilde_s = C44*Le[3*i+j]*Le[3*i+j]*(gamgradgrad[3*i+j]*(R_crystal*dir[i][j])).dot(R_crystal*dir[i][j]);
                Set::Scalar grad_ziebar_s = (tauegrad.dot(R_crystal*dir[i][j]))*(Ld[3*i+j]*Ld[3*i+j]*gamdotgrad[3*i+j].dot(R_crystal*dir[i][j]))/gamedot +
                                            taue*(Ld[3*i+j]*Ld[3*i+j]*(gamdotgradgrad[3*i+j]*R_crystal*dir[i][j]).dot(R_crystal*dir[i][j]))/gamedot - 
                                            taue*(Ld[3*i+j]*Ld[3*i+j]*gamdotgrad[3*i+j].dot(R_crystal*dir[i][j]))*(gamedotgrad.dot(R_crystal*dir[i][j]))/(gamedot*gamedot);

                gamdotalpha[3*i+j] -= visc*(qbar - tau - grad_zietilde_s - grad_ziebar_s);
            }
        }
    }

    // Hardening law
    // \dot{g}^{\alpha} = \sum_{\beta=1}^N h_{\alpha\beta} |\dot{\gamma}^\beta|
    // h_{\alpha\beta} = qh_{\alpha\beta} h(\bar{\gamma})
    // h(\bar{\gamma}) = h_s + (h_s - h_0) \sech^2 (\frac{h_0-h_s}{g_s-g_0} \bar{\gamma})
    // \bar{\gamma} = \int_0^t \sum_{\alpha=1}^N |\dot{\gamma}^\alpha| dt
    // qh_{\alpha\beta} = 1 (if \alpha,\beta coplanar), qh (else)
    void Advance(const Set::Matrix /*sigma*/, const Set::Matrix /*strain*/, const Set::Matrix3 /*epsgrad*/, Set::Scalar dt)
    {
        Set::Matrix_3D R_crystal;
        R_crystal = Model::Solid::Linear::Cubic::q.normalized().toRotationMatrix();

        std::array<Set::Scalar, 12> gamdot;
        
        for (int i=0; i< 4; i++)
        {
            for (int j=0; j<3; j++)
            {
                F0 += dt * gamdotalpha[3*i + j] * ((R_crystal*dir[i][j]) * (R_crystal*plane[i]).transpose());
                gamalpha[3*i+j] += dt*std::abs(gamdot[3*i + j]);
            }
        }

        //if (hardening)
    }

public:

    Set::Scalar visc;
    std::array<Set::Scalar,12> Le;
    std::array<Set::Scalar,12> Ld;
    std::array<Set::Scalar,12> gamalpha;
    std::array<Set::Scalar,12> gamdotalpha;
    
public:
    static StrainGradientCrystalPlastic Zero()
    {
        StrainGradientCrystalPlastic ret;
        ret.q.coeffs() *= 0;
        ret.C11 = 0.0; ret.C12 = 0.0; ret.C44 = 0.0;
        ret.F0 = Set::Matrix::Zero();
        ret.DefineSlipSystem(0.,0.,0.);
        ret.DefineSlipGradientSystem(0.,0.,0.);
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

        amrex::Vector<Set::Scalar> input_Ld, input_Le;
        std::array<Set::Scalar, 12> f_Le, f_Ld;
        if (pp.contains("Le"))
        {
            pp.queryarr("Le",input_Le);
            if(input_Le.size() == 1)
            {
                for (int i = 0; i< 12; i++) f_Le[i] = input_Le[0];
            }
            else if(input_Le.size() <= 0 || input_Le.size() != 12) Util::Abort(INFO, "Invalid number of values in Le");
            else
            {
                for (int i = 0; i< 12; i++) f_Le[i] = input_Le[i];
            }
        }
        else
            for (int i = 0; i< 12; i++) f_Le[i] = 0.1;
        
        if (pp.contains("Ld"))
        {
            pp.queryarr("Ld",input_Ld);
            if(input_Ld.size() == 1)
            {
                for (int i = 0; i< 12; i++) f_Ld[i] = input_Ld[0];
            }
            else if(input_Ld.size() <= 0 || input_Ld.size() != 12) Util::Abort(INFO, "Invalid number of values in Ld");
            else
            {
                for (int i = 0; i< 12; i++) f_Ld[i] = input_Ld[i];
            }
        }
        else
            for (int i = 0; i< 12; i++) f_Ld[i] = 0.1;
        
        Set::Scalar input_visc;
        pp.query("visc",input_visc);
        value.DefineSlipGradientSystem(f_Le,f_Ld,input_visc);
    }

    #define OP_CLASS StrainGradientCrystalPlastic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) X(g0) X(gs) X(h0) X(hs) X(qh) X(gam_bar) X(visc)
    #define OP_LOOP_NUM 12
    #define OP_LOOP_VARS X(m_exp) X(g_alpha) X(gam0dot) X(gamalpha) X(gamdotalpha) X(Le) X(Ld)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"


}
}
}

#endif
