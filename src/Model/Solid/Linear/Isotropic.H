#ifndef MODEL_SOLID_LINEAR_ISOTROPIC_H_
#define MODEL_SOLID_LINEAR_ISOTROPIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Isotropic : public Solid<Set::Sym::Isotropic>
{
public:

    Isotropic() {};
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda) 
    {
        Define(a_mu,a_lambda);
    };
    virtual ~Isotropic() {};

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        m_lambda = a_lambda; m_lambda0 = a_lambda;
        m_mu = a_mu; m_mu0 = a_mu;
    }

    void DegradeModulus(const Set::Scalar eta)
    {
        if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta ", eta);
        if(std::isnan(m_lambda0) || std::isnan(m_mu0)) 
            Util::Abort(INFO, "Nans detected","ddw=", ddw, " mu0=", m_mu0, " lambda0=",m_lambda0); 

        m_mu = m_mu0*(1.-eta);
        m_lambda = m_lambda0*(1.-eta);
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
        }
        else
            return;
    }

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (DDW(gradu)*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return DDW(gradu)*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(m_lambda,m_mu);
        //return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << "lambda = " << m_lambda << ", mu = " << m_mu << std::endl;
        out << "lambda0 = " << m_lambda0 << ", mu0 = " << m_mu0 << std::endl;
    }

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    Set::Scalar m_lambda = NAN, m_lambda0 = NAN;
    Set::Scalar m_mu = NAN, m_mu0 = NAN;
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> ddw;

public:
    static Isotropic Random()
    {
        Isotropic ret;
        ret.Define(Util::Random(),Util::Random());
        return ret;
    }
    static void Parse(Isotropic & value, IO::ParmParse & pp)
    {
        Set::Scalar mu, lambda;
        if (pp.contains("lame") && pp.contains("shear"))
        {
            pp.query("lame",lambda);
            pp.query("shear",mu);
        }
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
        value.Define(mu,lambda);
    }

    #define OP_CLASS Isotropic
    #define OP_VARS X(m_lambda) X(m_lambda0) X(m_mu) X(m_mu0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"


}
}
}

#endif
