#ifndef MODEL_SOLID_LINEARELASTIC_MULTIWELL_H_
#define MODEL_SOLID_LINEARELASTIC_MULTIWELL_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
class Multiwell : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	AMREX_GPU_HOST_DEVICE 
	Multiwell() {} ;
	
	AMREX_GPU_HOST_DEVICE 
	Multiwell(Set::Scalar _mu, Set::Scalar _lambda, Set::Matrix _Fgb=Set::Matrix::Zero()) {Define(_mu,_lambda,_Fgb);} ;
	
	AMREX_GPU_HOST_DEVICE 
	void Define(Set::Scalar _mu, Set::Scalar _lambda,Set::Matrix _Fgb)
	{
		mu = _mu; lambda = _lambda; Fgb = _Fgb;
		for (int p = 0; p < AMREX_SPACEDIM; p++) gradFgb[p] = Set::Matrix::Zero();
	}

	AMREX_FORCE_INLINE
	AMREX_GPU_HOST_DEVICE virtual Set::Matrix operator () (Set::Matrix &gradu,bool a_homogeneous=true) const
	{
		Set::Matrix eps = gradu;
		if (!a_homogeneous) eps -= Fgb;
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};

	AMREX_FORCE_INLINE
	AMREX_GPU_HOST_DEVICE virtual Set::Vector operator () (Set::Matrix3 &gradgradu,bool a_homogeneous=true)
	{
		Set::Matrix3 gradeps = gradgradu;
		if (!a_homogeneous) for (int i = 0; i < AMREX_SPACEDIM; i++) gradeps[i] -= gradFgb[i];
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				{
					f(i) += mu*(gradeps(i,j,j) + gradeps(j,i,j))  + lambda*gradeps(j,j,i);
				}
		return f;
	}

	void Randomize()
	{
		mu = Util::Random();
		lambda = Util::Random();
	};

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda << std::endl;
	}

	AMREX_GPU_HOST_DEVICE 
	Multiwell operator + (const Multiwell &rhs) const
	{
		//gradFgb + gradFgb;
		Multiwell ret(mu+rhs.mu, lambda+rhs.lambda,Fgb + rhs.Fgb); ret.gradFgb = gradFgb; ret.gradFgb += rhs.gradFgb;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i] + rhs.gradFgb[i];
		return ret;
	}

	AMREX_GPU_HOST_DEVICE 
	Multiwell & operator = (const Multiwell &rhs)
	{
		mu=rhs.mu; lambda=rhs.lambda; Fgb = rhs.Fgb; gradFgb = rhs.gradFgb;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) gradFgb[i] = rhs.gradFgb[i];
		return *this;
	}

	AMREX_GPU_HOST_DEVICE 
	void operator += (const Multiwell &rhs)
	{
		mu+=rhs.mu; lambda+=rhs.lambda; Fgb += rhs.Fgb; gradFgb += rhs.gradFgb;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) gradFgb[i] += rhs.gradFgb[i];
	}

	AMREX_GPU_HOST_DEVICE 
	Multiwell operator - (const Multiwell &rhs) const
	{
		Multiwell ret(mu-rhs.mu, lambda-rhs.lambda,Fgb-rhs.Fgb); ret.gradFgb = gradFgb; ret.gradFgb -= rhs.gradFgb;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i] - rhs.gradFgb[i];
		return ret;
	}

	//Multiwell operator * (const Multiwell &rhs) const
	//{return Multiwell(mu*rhs.mu, lambda*rhs.lambda,Fgb);}

	//Multiwell operator / (const Multiwell &rhs) const
	//{return Multiwell(mu/rhs.mu, lambda/rhs.lambda);}

	AMREX_GPU_HOST_DEVICE 
	virtual Multiwell operator / (const Set::Scalar alpha) const
	{
		Multiwell ret(mu/alpha, lambda/alpha,Fgb/alpha); ret.gradFgb = gradFgb; ret.gradFgb /= alpha;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i]/alpha;
		return ret;
	}

	AMREX_GPU_HOST_DEVICE 
	virtual Multiwell operator * (const Set::Scalar alpha) const
	{
		Multiwell ret(mu*alpha, lambda*alpha,Fgb*alpha); ret.gradFgb = gradFgb; ret.gradFgb *= alpha;
		//for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i]*alpha;
		return ret;
	}
	
private:
	Set::Scalar mu; ///< First Lame parameter \f$\mu\f$
	Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
public:
	Set::Matrix Fgb = Set::Matrix::Zero();
	Set::Matrix3 gradFgb;
};

}
}
}

#endif


