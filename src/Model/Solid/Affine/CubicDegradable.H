#ifndef MODEL_SOLID_AFFINE_CUBIC_DEGRADABLE_H_
#define MODEL_SOLID_AFFINE_CUBIC_DEGRADABLE_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Affine/Affine.H"
#include "Model/Solid/Affine/Cubic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class CubicDegradable : public Cubic
{
public:

    CubicDegradable() {};
    CubicDegradable(Cubic base) : Cubic(base) {};

    void DegradeModulus(const Set::Scalar eta)
    {
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
        if(std::isnan(C110)) Util::Abort(INFO, "Nans in C110");
        if(std::isnan(C120)) Util::Abort(INFO, "Nans in C120");
        if(std::isnan(C440)) Util::Abort(INFO, "Nans in C440");

        C11 = C110*(1.-eta);
        C12 = C120*(1.-eta);
        C44 = C440*(1.-eta);
    }
    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
    {
        if(eta.size() == 0) return;
        else if (eta.size() == 1) DegradeModulus(eta[0]);
        else if(eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            C11 = C110*(1.-eta[0]);
            C12 = C120*(1.-eta[1]);
            C44 = C440*(1.-eta[1]);
        }
        else if (eta.size() == 3)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1]) || std::isnan(eta[2])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0 || eta[2] < 0.0 || eta[2] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            C11 = C110*(1.-eta[0]);
            C12 = C120*(1.-eta[1]);
            C44 = C440*(1.-eta[2]);
        }
        else
            return;
    }

public:
    Set::Scalar C110, C120, C440;

public:

    static CubicDegradable Zero()
    {
        CubicDegradable ret;
        ret.C11 = 0.0;
        ret.C12 = 0.0;
        ret.C44 = 0.0;
        ret.F0 = Set::Matrix::Zero();
        ret.q.coeffs() *= 0;
        ret.C110 = 0.0;
        ret.C120 = 0.0;
        ret.C440 = 0.0;
        return ret;
    }

    virtual void Print(std::ostream &out) const override 
    {
        out << "C11     = " << C11 << std::endl;
        out << "C110     = " << C110 << std::endl;
        out << "C12     = " << C12 << std::endl;
        out << "C120     = " << C120 << std::endl;
        out << "C44     = " << C44 << std::endl;
        out << "C440     = " << C440 << std::endl;
        out << "q      = " << q.coeffs() << std::endl();
        out << "F0 \n" << F0 << std::endl();
    }

    static CubicDegradable Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }

    static CubicDegradable Random(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44)
    {
        CubicDegradable ret = Cubic::Random(a_C11,a_C12,a_C44);
        ret.C110 = ret.C11;
        ret.C120 = ret.C12;
        ret.C440 = ret.C44;
        return ret;
    }

    // This class inherits from :ref:`Model::Solid::Affine::Cubic`.
    // It provides the ability to "degrade" while retaining information
    // about its original, pristine state
    static void Parse(CubicDegradable & value, IO::ParmParse & pp)
    {
        Cubic::Parse(value,pp);
        value.C110 = value.C11;
        value.C120 = value.C11;
        value.C440 = value.C44;
    }

    #define OP_CLASS CubicDegradable
    #define OP_VARS  X(C11) X(C12) X(C44) X(q.coeffs()) X(F0) X(C110) X(C120) X(C440)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}
#endif
