#ifndef MODEL_SOLID_AFFINE_ISOTROPIC_H_
#define MODEL_SOLID_AFFINE_ISOTROPIC_H_

#include "AMReX.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Linear/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Isotropic : public Linear::Isotropic
{
public:

    Isotropic() {};
    Isotropic(Linear::Isotropic base) : Linear::Isotropic(base) {};

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return Linear::Isotropic::W(gradu - F0); //0.5*((F-F0).transpose() * (ddw*((F-F0)))).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return Linear::Isotropic::DW(gradu - F0); // ddw*(F-F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & gradu) const override
    {
        return Linear::Isotropic::DDW(gradu - F0); // ddw;
    }

    Set::Matrix F0 = Set::Matrix::Zero();
    static const KinematicVariable kinvar = KinematicVariable::gradu;

public:
    static Isotropic Zero()
    {
        Isotropic ret;
        ret.Define(0.,0.);
        ret.F0 = Set::Matrix::Zero();
        return ret;
    }
    static Isotropic Random()
    {
        Isotropic ret;
        ret.Define(Util::Random(), Util::Random());
        ret.F0 = Set::Matrix::Random();
        return ret;
    }
    static void Parse(Isotropic & value, IO::ParmParse & pp)
    {
        Linear::Isotropic::Parse(value,pp);
        if (pp.contains("F0")) pp.queryarr("F0",value.F0);
    }
    #define OP_CLASS Isotropic
    #define OP_VARS  X(m_lambda) X(m_lambda0) X(m_mu) X(m_mu0) X(F0)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

#endif
