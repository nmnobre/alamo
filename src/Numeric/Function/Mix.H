#ifndef NUMERIC_FUNCTION_MIX_H
#define NUMERIC_FUNCTION_MIX_H
#include <AMReX.H>
#include <utility>

namespace Numeric
{
namespace Function
{

template <class T> class Mix
{
    Mix() {};
    T operator () (T rhs) const = 0;
};

template<>
class Mix<std::pair<Set::Scalar,Set::Scalar>>
{
    public:
    Mix() {}
    AMREX_FORCE_INLINE
    std::pair<Set::Scalar,Set::Scalar> operator () (std::pair<Set::Scalar,Set::Scalar> &x) const
    {
        std::pair<Set::Scalar,Set::Scalar> ret;
        Set::Scalar x1 = x.first, x2 = x.second;
        ret.first  = x1*x1;
        ret.second = x2*x2;
        Set::Scalar den = ret.first + ret.second;
        ret.first /= den;
        ret.second /= den;
        return ret;
    }

    Set::Scalar D1(Set::Scalar x1, Set::Scalar x2) const
    {
        Set::Scalar den = x1*x1 + x2*x2;
        return (x1*x2*x2) / (den*den);
    }
    Set::Scalar D2(Set::Scalar x1, Set::Scalar x2) const
    {
        Set::Scalar den = x1*x1 + x2*x2;
        return (x1*x1*x2) / (den*den);
    }
};
}
}







#endif