#ifndef MODEL_SOLID_AFFINE_ISOTROPICDEGRADABLE_H_
#define MODEL_SOLID_AFFINE_ISOTROPICDEGRADABLE_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/IsotropicDegradable.H"
//#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class IsotropicDegradable : public Linear::IsotropicDegradable
{
public:
    IsotropicDegradable() {};
    IsotropicDegradable(Linear::IsotropicDegradable base) : Linear::IsotropicDegradable(base) {};

    void Define(Set::Scalar a_mu, Set::Scalar a_mu0, Set::Scalar a_lambda, Set::Scalar a_lambda0)
    {
        Linear::IsotropicDegradable::Define(a_mu, a_mu0, a_lambda, a_lambda0);
    }

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return Linear::IsotropicDegradable::W(gradu - F0); //0.5*((F-F0).transpose() * (ddw*((F-F0)))).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return Linear::IsotropicDegradable::DW(gradu - F0); // ddw*(F-F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & gradu) const override
    {
        return Linear::IsotropicDegradable::DDW(gradu - F0); // ddw;
    }

    AMREX_FORCE_INLINE
    void SetF0(Set::Matrix a_F0) {F0 = a_F0; }

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;
    Set::Matrix F0 = Set::Matrix::Zero();

public:
    static IsotropicDegradable Zero()
    {
        IsotropicDegradable ret;
        ret.Define(0.,0.,0.,0.);
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }

    static IsotropicDegradable Random()
    {
        IsotropicDegradable ret;
        ret.Define(Util::Random(), Util::Random(), Util::Random(), Util::Random());
        ret.F0 = Set::Matrix::Random();
        return ret;
    }

    static void Parse(IsotropicDegradable & value, IO::ParmParse & pp)
    {
        Linear::IsotropicDegradable::Parse(value,pp);
        if (pp.contains("F0")) pp.queryarr("F0",value.F0);
    }

    virtual void Print(std::ostream &out) const override 
    {
        out << "ddw     = " << ddw << std::endl;
        out << "mu0     = " << m_mu0 << std::endl;
        out << "lambda0 = " << m_lambda0 << std::endl;
        out << "F0 = \n" << F0 << std::endl;
    }
    
    #define OP_CLASS IsotropicDegradable
    #define OP_VARS  X(ddw) X(m_mu0) X(m_lambda0) X(m_mu) X(m_lambda) X(F0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}
#endif
