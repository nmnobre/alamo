#ifndef MODEL_SOLID_PLASTIC_J2_H_
#define MODEL_SOLID_PLASTIC_J2_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class J2Plastic : public Isotropic
{
public:
    J2Plastic() {};
    //J2Plastic(Solid<Set::Sym::Isotropic> base) : Isotropic(base) {};

    // J2Plastic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta=1.0)
    // {
    //     Define(a_mu, a_lambda, a_yield, a_hardening, a_theta);
    // }

    void DefineSlipSystem(Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta, Set::Scalar a_mu)
    {
        theta = a_theta;
        yield_strength = a_yield;
        hardening_modulus = a_hardening;
        m_mu = a_mu;
        
        F0 = Set::Matrix::Zero(); 
        beta = Set::Matrix::Zero();
    }

    Set::Scalar YieldSurface()
    {
        return yield_strength + F0.norm()*hardening_modulus;
    }

    Set::Scalar PlasticEnergy()
    {
        return (yield_strength*F0.norm() + 0.5*hardening_modulus*F0.norm()*F0.norm());
    }
    
    void Advance(Set::Matrix sigma, Set::Matrix strain, Set::Matrix3 /*eps_grad*/, Set::Scalar /*dt*/)
    {
        Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));
        Set::Matrix sigdev = sigma - (1.0/((double)AMREX_SPACEDIM))*sigma.trace()*Set::Matrix::Identity();
        Set::Matrix epsdev = strain - strain.trace()*Set::Matrix::Identity();

        Set::Matrix zeta_trial = sigdev  - beta + 2.0*m_mu*epsdev;

        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*F0.norm());
        if( f_trial <= 0.0)
        {
            // curr = prev; 
            // beta = beta_prev;
            // curr.alpha = prev.alpha; curr.beta = prev.beta; 
            return;
        }
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(m_mu)*(1.0 + (hardening_modulus/(3.0*m_mu))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        //curr.alpha = prev.alpha + SQ2O3*dGamma;
        beta += SQ2O3*dH*n_new;
        F0 += dGamma*n_new;

        if(std::isnan(F0.norm()))
        {
            Util::Message(INFO, "f_trial = ", f_trial);
            Util::Abort(INFO, "Nans detected");
        }
        //Isotropic::SetF0(curr);
    }

    void SetPlasticStrains()
    {
        F0 = Set::Matrix::Zero();
    }

    void SetPlasticStrains(Set::Matrix &a_state)
    {
        F0 = a_state;
    }

    Set::Matrix GetPlasticState()
    {
        return F0;
    }

public:
    Set::Matrix beta;
    Set::Scalar theta;                  // isotropic and kinematic hardening parameter
    Set::Scalar yield_strength;          // yield strength
    Set::Scalar hardening_modulus;      // hardening modulus
    Set::Scalar m_mu;

public:
    static J2Plastic Zero()
    {
        J2Plastic ret;
        ret.Define(0.,0.);
        ret.DefineSlipSystem(0.,0.,0.,0.);
        return ret;
    }
    static J2Plastic Random()
    {
        J2Plastic ret;
        Set::Scalar mu = Util::Random();
        ret.Define(mu,Util::Random());
        ret.DefineSlipSystem(Util::Random(),Util::Random(), Util::Random(),mu);
        ret.F0 = Set::Matrix::Random();
        ret.beta = Set::Matrix::Random();
        return ret;
    }
    static void Parse(J2Plastic & value, IO::ParmParse & pp)
    {
        Isotropic::Parse(value,pp);
        Set::Scalar mu, lambda;
        if (pp.contains("lambda") && pp.contains("mu"))
        {
            pp.query("lambda",lambda);
            pp.query("mu",mu);
        }
        else if (pp.contains("E") && pp.contains("nu"))
        {
            Set::Scalar E, nu;
            pp.query("E",E);
            pp.query("nu",nu);
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        Set::Scalar yield = 1.0, hardening = 1.0, theta_tmp = 1.0;

        pp.query("yield", yield);
        pp.query("hardening", hardening);
        pp.query("theta", theta_tmp);

        value.DefineSlipSystem(yield, hardening, theta_tmp, mu);
    }
    #define OP_CLASS J2Plastic
    #define OP_VARS X(ddw) X(F0) X(beta) X(theta) X(yield_strength) X(hardening_modulus) X(m_mu)
    #include "Model/Solid/InClassOperators.H"

};
#include "Model/Solid/ExtClassOperators.H"


}
}
}
#endif