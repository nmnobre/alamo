#ifndef MODEL_SOLID_AFFINE_H_
#define MODEL_SOLID_AFFINE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
    enum LatticeType {FCC, BCC, HCP};   // standard Bravais lattices. More to be added later
//class Affine
//{

// struct PlasticState{
//     Set::Matrix epsp;   // Plastic strain
//     Set::Scalar alpha;  // Equivalent plastic strain
//     Set::Matrix beta;

//     static PlasticState Zero()
//     {
//         PlasticState ret;
//         ret.epsp = Set::Matrix::Zero();
//         ret.beta = Set::Matrix::Zero();
//         ret.alpha = 0.0;
//         return ret;
//     }

// //    AMREX_FORCE_INLINE
// //    void operator = (const PlasticState &rhs)
// //    {
// //        epsp = rhs.epsp;
// //        beta = rhs.beta;
// //        alpha = rhs.alpha;
// //    }

//     AMREX_FORCE_INLINE
//     void operator += (const PlasticState &rhs)
//     {
//         epsp += rhs.epsp;
//         alpha += rhs.alpha;
//         beta += rhs.beta;
//     }
//     AMREX_FORCE_INLINE
//     PlasticState operator * (const Set::Scalar _alpha) const
//     {
//         PlasticState ret;
//         ret.epsp = epsp*_alpha;
//         ret.alpha = alpha*_alpha;
//         ret.beta = beta*_alpha;
//         return ret;
//     }
//     friend PlasticState operator * (const Set::Scalar alpha, const PlasticState b);
//     friend PlasticState operator + (const PlasticState a, const PlasticState b);
//     friend PlasticState operator - (const PlasticState a, const PlasticState b);
// };

// AMREX_FORCE_INLINE
// PlasticState operator * (const Set::Scalar alpha, const PlasticState b)
// {
//     PlasticState ret;
//     ret.epsp = b.epsp*alpha;
//     ret.beta = b.beta*alpha;
//     ret.alpha = b.alpha*alpha;
//     return ret;
// }

// AMREX_FORCE_INLINE
// PlasticState operator + (const PlasticState a, const PlasticState b)
// {
//     PlasticState ret;
//     ret.epsp = a.epsp + b.epsp;
//     ret.beta = a.beta + b.beta;
//     ret.alpha = a.alpha + b.alpha;
//     return ret;
// }

// AMREX_FORCE_INLINE
// PlasticState operator - (const PlasticState a, const PlasticState b)
// {
//     PlasticState ret;
//     ret.epsp = a.epsp - b.epsp;
//     ret.beta = a.beta - b.beta;
//     ret.alpha = a.alpha - b.alpha;
//     return ret;
// }


// struct SlipSystem{
//     bool on;                ///< slip system active or not
//     Set::Vector s;          ///< slip plane normal
//     Set::Vector m;          ///< slip direction
//     Set::Scalar gam0dot;    ///< gamma0 dot 
//     Set::Scalar g_alpha;    ///< current strength
//     Set::Scalar m_exp;      ///< exponent

//     Set::Scalar g_alpha0;     ///< pristine material current strength

//     void Equal (const SlipSystem &rhs)
//     {
//         s = rhs.s;
//         m = rhs.m;
//         gam0dot = rhs.gam0dot;
//         g_alpha = rhs.g_alpha;
//         g_alpha0 = rhs.g_alpha0;
//         m_exp = rhs.m_exp;
//     }

// //    AMREX_FORCE_INLINE
// //    void operator = (const SlipSystem &rhs)
// //    {
// //        n = rhs.n;
// //        s = rhs.s;
// //        tcrss = rhs.tcrss;
// //        tcrss0 = rhs.tcrss0;
// //        s_exp = rhs.s_exp;
// //        gam0dot = rhs.gam0dot;
// //    }

//     AMREX_FORCE_INLINE
//     void operator += (const SlipSystem &rhs)
//     {
//         s += rhs.s;
//         m += rhs.m;
//         g_alpha += rhs.g_alpha;
//         g_alpha0 += rhs.g_alpha0;
//         m_exp += rhs.m_exp;
//         gam0dot += rhs.gam0dot;
//     }
//     AMREX_FORCE_INLINE
//     SlipSystem operator * (const Set::Scalar alpha) const
//     {
//         SlipSystem ret;
//         ret.s = s*alpha;
//         ret.m = m*alpha;
//         ret.g_alpha = g_alpha*alpha;
//         ret.g_alpha0 = g_alpha0*alpha;
//         ret.m_exp = m_exp*alpha;
//         ret.gam0dot = gam0dot*alpha;
        
//         return ret;
//     }
//     friend SlipSystem operator * (const Set::Scalar alpha, const SlipSystem b);
//     friend SlipSystem operator + (const SlipSystem a, const SlipSystem b);
//     friend SlipSystem operator - (const SlipSystem a, const SlipSystem b);
// };
// AMREX_FORCE_INLINE
// SlipSystem operator * (const Set::Scalar alpha, const SlipSystem b)
// {
//     SlipSystem ret;
//     ret.s = b.s*alpha;
//     ret.m = b.m*alpha;
//     ret.g_alpha = b.g_alpha*alpha;
//     ret.g_alpha0 = b.g_alpha0*alpha;
//     ret.m_exp = b.m_exp*alpha;
//     ret.gam0dot = b.gam0dot*alpha;
    
//     return ret;
// }
// AMREX_FORCE_INLINE
// SlipSystem operator + (const SlipSystem a, const SlipSystem b)
// {
//     SlipSystem ret;
//     ret.s = a.s + b.s;
//     ret.m = a.m + b.m;
//     ret.g_alpha = a.g_alpha + b.g_alpha;
//     ret.g_alpha0 = a.g_alpha0 + b.g_alpha0;
//     ret.m_exp = a.m_exp + b.m_exp;
//     ret.gam0dot = a.gam0dot + b.gam0dot;
//     return ret;
// }
// AMREX_FORCE_INLINE
// SlipSystem operator - (const SlipSystem a, const SlipSystem b)
// {
//     SlipSystem ret;
//     ret.s = a.s - b.s;
//     ret.m = a.m - b.m;
//     ret.g_alpha = a.g_alpha - b.g_alpha;
//     ret.g_alpha0 = a.g_alpha0 - b.g_alpha0;
//     ret.m_exp = a.m_exp - b.m_exp;
//     ret.gam0dot = a.gam0dot - b.gam0dot;
//     return ret;
// }
    
//};
}
}
}
#endif