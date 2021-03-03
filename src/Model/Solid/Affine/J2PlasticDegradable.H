#ifndef MODEL_SOLID_PLASTIC_J2DEGRADABLE_H_
#define MODEL_SOLID_PLASTIC_J2DEGRADABLE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/IsotropicDegradable.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class J2PlasticDegradable : public IsotropicDegradable
{
public:
    J2PlasticDegradable() {};
    J2PlasticDegradable(Solid<Set::Sym::Isotropic> base) : IsotropicDegradable(base) {};

    J2PlasticDegradable(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta=1.0)
    {
        Define(a_mu, a_lambda, a_yield, a_hardening, a_theta);
    }

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening, Set::Scalar a_theta)
    {
        theta = a_theta;
        yield_strength = a_yield;
        yield_strength0 = a_yield;
        hardening_modulus = a_hardening;
        hardening_modulus0 = a_hardening;

        curr = Set::Matrix::Zero(); 
        prev = Set::Matrix::Zero();
        beta = Set::Matrix::Zero();
        beta_prev = Set::Matrix::Zero();

        IsotropicDegradable::Define(a_mu, a_lambda);
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
        return yield_strength + curr.norm()*hardening_modulus;
    }

    Set::Scalar OriginalYieldSurface()
    {
        return yield_strength0 + curr.norm()*hardening_modulus0;
    }

    Set::Scalar PlasticEnergy()
    {
        return (yield_strength*curr.norm() + 0.5*hardening_modulus*curr.norm()*curr.norm());
    }
    
    Set::Scalar OriginalPlasticEnergy()
    {
        return (yield_strength0*curr.norm() + 0.5*hardening_modulus0*curr.norm()*curr.norm());
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

        prev = curr;
        beta_prev = beta;
        Set::Matrix zeta_trial = sigdev  - beta_prev + 2.0*m_mu*epsdev;

        Set::Scalar f_trial = zeta_trial.norm() - SQ2O3*(yield_strength + theta*hardening_modulus*prev.norm());
        if( f_trial <= 0.0)
        {
            curr = prev; 
            beta = beta_prev;
            //curr.alpha = prev.alpha; curr.beta = prev.beta; 
            return;
        }
        Set::Matrix n_new = zeta_trial/zeta_trial.norm();
        Set::Scalar dGamma = f_trial/(2.0*(m_mu)*(1.0 + (hardening_modulus/(3.0*m_mu))));
        Set::Scalar dH = SQ2O3*(1.0-theta)*hardening_modulus*dGamma;

        //curr.alpha = prev.alpha + SQ2O3*dGamma;
        beta = beta_prev + SQ2O3*dH*n_new;
        curr = prev + dGamma*n_new;
        IsotropicDegradable::SetF0(curr);
    }

    void SetPlasticStrains()
    {
        prev = Set::Matrix::Zero();
        curr = Set::Matrix::Zero();
    }

    void SetPlasticStrains(Set::Matrix &a_state)
    {
        prev = a_state;
        curr = a_state;
        IsotropicDegradable::SetF0(curr);
    }

    Set::Matrix GetPlasticState()
    {
        return curr;
    }

public:
    Set::Matrix curr, prev;
    Set::Matrix beta, beta_prev;
    Set::Scalar theta;                  // isotropic and kinematic hardening parameter
    Set::Scalar yield_strength, yield_strength0;          // yield strength
    Set::Scalar hardening_modulus, hardening_modulus0;  // hardening modulus

public:

    static J2PlasticDegradable Zero()
    {
        J2PlasticDegradable ret;
        ret.Define(0.,0.,0.,0.,0.0);
        return ret;
    }
    static void Parse(J2PlasticDegradable & value, IO::ParmParse & pp)
    {
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

        value.Define(mu, lambda, yield, hardening, theta_tmp);
    }
    AMREX_FORCE_INLINE
	void operator += (const J2PlasticDegradable &rhs)
    {
        curr += rhs.curr;
        prev += rhs.prev;
        beta += rhs.beta;
        beta_prev += rhs.beta_prev;
        theta += rhs.theta;
        yield_strength0 += rhs.yield_strength0;
        yield_strength += rhs.yield_strength;
        hardening_modulus += rhs.hardening_modulus;
        hardening_modulus0 += rhs.hardening_modulus0;

        ddw += rhs.ddw;
        F0 += rhs.F0;
        _F0 += rhs._F0;
        mu0 += rhs.mu0;
		lambda0 += rhs.lambda0;
	}    
    AMREX_FORCE_INLINE
   	J2PlasticDegradable operator * (const Set::Scalar alpha) const
	{
        J2PlasticDegradable ret;
        ret.curr = alpha*curr;
        ret.prev = alpha*prev;
        ret.beta = alpha*beta;
        ret.beta_prev = alpha*beta_prev;
        ret.theta = alpha*theta;
        ret.yield_strength = alpha*yield_strength;
        ret.yield_strength0 = alpha*yield_strength0;
        ret.hardening_modulus0 = alpha*hardening_modulus0;
        ret.hardening_modulus = alpha*hardening_modulus;

        ret.ddw = ddw*alpha;
        ret.F0 = F0*alpha;
        ret._F0 = _F0*alpha;
        ret.mu0 = mu0*alpha;
		ret.lambda0 = lambda0*alpha;

        return ret;
	}    
    friend J2PlasticDegradable operator * (const Set::Scalar alpha, const J2PlasticDegradable b);
    friend J2PlasticDegradable operator + (const J2PlasticDegradable a, const J2PlasticDegradable b);
    friend J2PlasticDegradable operator - (const J2PlasticDegradable a, const J2PlasticDegradable b);

};
AMREX_FORCE_INLINE
J2PlasticDegradable operator * (const Set::Scalar alpha, const J2PlasticDegradable b)
{
    J2PlasticDegradable ret;
    ret.curr = alpha*b.curr;
    ret.prev = alpha*b.prev;
    ret.beta = alpha*b.beta;
    ret.beta_prev = alpha*b.beta_prev;
    ret.theta = alpha*b.theta;
    ret.yield_strength = alpha*b.yield_strength;
    ret.yield_strength0 = alpha*b.yield_strength0;
    ret.hardening_modulus0 = alpha*b.hardening_modulus0;
    ret.hardening_modulus = alpha*b.hardening_modulus;

    ret.ddw = b.ddw*alpha;
    ret.F0 = b.F0*alpha;
    ret._F0 = b._F0*alpha;
    ret.mu0 = b.mu0*alpha;
	ret.lambda0 = b.lambda0*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
J2PlasticDegradable operator + (const J2PlasticDegradable a, const J2PlasticDegradable b)
{
    J2PlasticDegradable ret;
    ret.curr = a.curr + b.curr;
    ret.prev = a.prev + b.prev;
    ret.beta = a.beta + b.beta;
    ret.beta_prev = a.beta_prev + b.beta_prev;
    ret.theta = a.theta + b.theta;
    ret.yield_strength0 = a.yield_strength0 + b.yield_strength0;
    ret.yield_strength = a.yield_strength + b.yield_strength;
    ret.hardening_modulus = a.hardening_modulus + b.hardening_modulus;
    ret.hardening_modulus0 = a.hardening_modulus0 + b.hardening_modulus0;

    ret.ddw = a.ddw + b.ddw;
    ret.F0 = a.F0 + b.F0;
    ret._F0 = a._F0 + b._F0;
	ret.mu0 = a.mu0 + b.mu0;
	ret.lambda0 = a.lambda0 + b.lambda0;
    return ret;
}    
AMREX_FORCE_INLINE
J2PlasticDegradable operator - (const J2PlasticDegradable a, const J2PlasticDegradable b)
{
    J2PlasticDegradable ret;
    ret.curr = a.curr - b.curr;
    ret.prev = a.prev - b.prev;
    ret.beta = a.beta - b.beta;
    ret.beta_prev = a.beta_prev - b.beta_prev;
    ret.theta = a.theta - b.theta;
    ret.yield_strength0 = a.yield_strength0 - b.yield_strength0;
    ret.yield_strength = a.yield_strength - b.yield_strength;
    ret.hardening_modulus = a.hardening_modulus - b.hardening_modulus;
    ret.hardening_modulus0 = a.hardening_modulus0 - b.hardening_modulus0;

    ret.ddw = a.ddw - b.ddw;
    ret.F0 = a.F0 - b.F0;
    ret._F0 = a._F0 - b._F0;
	ret.mu0 = a.mu0 - b.mu0;
	ret.lambda0 = a.lambda0 - b.lambda0;
    return ret;
} 


}
}
}
#endif