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
    
    void DefineSlipSystem(Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta)
    {
        theta = a_theta;
        yield_strength = a_yield; yield_strength0 = a_yield;
        hardening_modulus = a_hardening; hardening_modulus0 = a_hardening;
        
        F0 = Set::Matrix::Zero(); 
        beta = Set::Matrix::Zero();
    }

    void DegradeYieldSurface(const Set::Scalar eta)
    {
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
        if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
        if(std::isinf(eta)) Util::Abort(INFO, "Infs in eta");
        yield_strength = yield_strength0*(1. - eta);
        hardening_modulus = hardening_modulus0*(1.-eta);
    }

    void DegradeYieldSurface(const amrex::Vector<Set::Scalar> &eta)
    {
        if (eta.size() == 0) return;
        else if (eta.size() == 1)
        {
            if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
            yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[0]);
        }
        else if (eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[1]);
        }
        else return;
    }

    Set::Scalar YieldSurface()
    {
        return yield_strength + F0.norm()*hardening_modulus;
    }

    Set::Scalar OriginalYieldSurface()
    {
        return yield_strength0 + F0.norm()*hardening_modulus0;
    }

    Set::Scalar PlasticEnergy()
    {
        return (yield_strength*F0.norm() + 0.5*hardening_modulus*F0.norm()*F0.norm());
    }
    
    Set::Scalar OriginalPlasticEnergy()
    {
        return (yield_strength0*F0.norm() + 0.5*hardening_modulus0*F0.norm()*F0.norm());
    }

    Set::Scalar OriginalPlasticEnergy(Set::Scalar alpha)
    {
        return (yield_strength0*alpha + 0.5*hardening_modulus0*alpha*alpha);
    }
    
    void Advance(Set::Matrix sigma, Set::Matrix strain, Set::Matrix3 /*eps_grad*/, Set::Scalar /*dt*/)
    {
        Set::Scalar SQ2O3 = sqrt(1.0 - 1.0/((double)AMREX_SPACEDIM));
        Set::Matrix sigdev = sigma - (1.0/((double)AMREX_SPACEDIM))*sigma.trace()*Set::Matrix::Identity();
        Set::Matrix epsdev = strain - strain.trace()*Set::Matrix::Identity();

        Set::Matrix zeta_trial = sigdev  - beta + 2.0*m_mu*epsdev;

        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*F0.norm());
        if( f_trial <= 0.0) return;
        
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(m_mu)*(1.0 + (hardening_modulus/(3.0*m_mu))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        beta += SQ2O3*dH*n_new;
        F0 += dGamma*n_new;

        if(std::isnan(F0.norm()))
        {
            Util::Message(INFO, "f_trial = ", f_trial);
            Util::Abort(INFO, "Nans detected");
        }
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
    Set::Scalar yield_strength, yield_strength0;          // yield strength
    Set::Scalar hardening_modulus, hardening_modulus0;  // hardening modulus

public:
    static J2Plastic Zero()
    {
        J2Plastic ret;
        ret.Define(0.,0.);
        ret.DefineSlipSystem(0.,0.,0.);
        return ret;
    }
    static J2Plastic Random()
    {
        J2Plastic ret;
        ret.Define(Util::Random(),Util::Random());
        ret.DefineSlipSystem(Util::Random(),Util::Random(), Util::Random());
        ret.F0 = Set::Matrix::Random();
        ret.beta = Set::Matrix::Random();
        return ret;
    }
    static void Parse(J2Plastic & value, IO::ParmParse & pp)
    {
        Isotropic::Parse(value,pp);
        Set::Scalar yield = 1.0, hardening = 1.0, theta_tmp = 1.0;

        pp.query("yield", yield);
        pp.query("hardening", hardening);
        pp.query("theta", theta_tmp);

        value.DefineSlipSystem(yield, hardening, theta_tmp);
    }
    #define OP_CLASS J2Plastic
    #define OP_VARS X(m_lambda) X(m_lambda0) X(m_mu) X(m_mu0) X(F0) X(beta) X(theta) X(yield_strength) X(hardening_modulus) X(yield_strength0) X(hardening_modulus0)
    #include "Model/Solid/InClassOperators.H"

};
#include "Model/Solid/ExtClassOperators.H"


}
}
}
#endif
