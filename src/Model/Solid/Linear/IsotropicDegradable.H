#ifndef MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLE_H_
#define MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLE_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Linear/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class IsotropicDegradable : public Isotropic
{
public:
    IsotropicDegradable() {};
	IsotropicDegradable(Isotropic base) : Isotropic(base) {};
	IsotropicDegradable(Set::Scalar a_mu, Set::Scalar a_lambda) 
    {
        Define(a_mu, a_mu, a_lambda, a_lambda);
    };

	void Define(Set::Scalar a_mu, Set::Scalar a_mu0, Set::Scalar a_lambda, Set::Scalar a_lambda0)
	{
		m_mu0 = a_mu0; m_lambda0 = a_lambda0;
		m_mu = a_mu; m_lambda = a_lambda;
        Isotropic::Define(a_mu,a_lambda);
	}

    virtual ~IsotropicDegradable() {};

    void DegradeModulus(const Set::Scalar eta)
	{
		if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta ", eta);
		if(std::isnan(m_lambda0) || std::isnan(m_mu0)) 
			Util::Abort(INFO, "Nans detected","ddw=", ddw, " mu0=", m_mu0, " lambda0=",m_lambda0); 

		m_mu = m_mu0*(1.-eta);
		m_lambda = m_lambda0*(1.-eta);
		Isotropic::Define( m_mu, m_lambda );
	}

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
		if(eta.size() == 0) return;
		else if(eta.size() == 1) DegradeModulus(eta[0]);
		else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
			m_mu = m_mu0*(1.-eta[0]);
			m_lambda = m_lambda0*(1.-eta[1]);
			Isotropic::Define( m_mu, m_lambda );
		}
		else
			return;
	}

public:
    Set::Scalar m_mu0 = NAN, m_lambda0 = NAN;
	Set::Scalar m_mu = NAN, m_lambda = NAN;

public:
	static void Parse(IsotropicDegradable & value, IO::ParmParse & pp)
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
		value.Define(mu,mu,lambda,lambda);
	}

	#define OP_CLASS IsotropicDegradable
    #define OP_VARS  X(ddw) X(m_mu0) X(m_lambda0) X(m_mu) X(m_lambda)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
