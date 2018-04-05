#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Isotropic.H"


Operator::Elastic::PolyCrystal::Isotropic::Isotropic()
{
  lambda1 = 1.0;
  lambda2 = 1.0;
  mu1 = 1.0;
  mu2 = 2.0;
}

amrex::Real
Operator::Elastic::PolyCrystal::Isotropic::C(const int i, const int j, const int k, const int l,
					     const amrex::IntVect loc,
					     const int amrlev, const int mglev, const MFIter &mfi) const
{
  amrex::Real mu, lambda;

  const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);

  if (fabs(etafab(loc,0)) < 0.001 && fabs(etafab(loc,1)) < 0.001)
    amrex::Abort("etas are zero");

  mu = (mu1*etafab(loc,0) + mu2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  lambda = (lambda1*etafab(loc,0) + lambda2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));

  // TODO: This is a hack, needs to be fixed.
  if (mu != mu) mu = 0.5*(mu1+mu2);
  if (lambda != lambda) lambda = 0.5*(lambda1+lambda2);

  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  return ret;
}
