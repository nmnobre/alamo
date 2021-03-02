#include "Hydro.H"
#include "IO/ParmParse.H"
#include "IC/Constant.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"

namespace Integrator
{
Hydro::Hydro() : Integrator()
{
  IO::ParmParse pp("physics"); 
  pp.query("gamma", gamma);

                {
			IO::ParmParse pp("bc");
			rhoBC = new BC::Constant(1);
			pp.queryclass("rho", *static_cast<BC::Constant *>(rhoBC));
			pBC = new BC::Constant(1);
			pp.queryclass("p", *static_cast<BC::Constant *>(pBC));
			eBC = new BC::Constant(1);
			pp.queryclass("e", *static_cast<BC::Constant *>(eBC));
		}


  RegisterNewFab(rho_mf,     rhoBC, 1, 0, "rho", true);
  RegisterNewFab(rho_old_mf, rhoBC, 1, 0, "rho_old", false);
  RegisterNodalFab(u_mf, 3, 0, "u", true);
  RegisterNodalFab(u_old_mf, 3, 0, "u_old", false);
  RegisterNewFab(p_mf,     pBC, 1, 0, "p", true);
  RegisterNewFab(p_old_mf, pBC, 1, 0, "p_old", false);
  RegisterNewFab(e_mf,     eBC, 1, 0, "e", true);
  RegisterNewFab(e_old_mf, eBC, 1, 0, "e_old", false);
}

void Hydro::Initialize(int lev)
{        
        rho_mf[lev]->setVal(1.0);
	rho_old_mf[lev]->setVal(1.0);

        p_mf[lev]->setVal(1.0);
	p_old_mf[lev]->setVal(1.0);

	e_mf[lev]->setVal(1.0);
	e_old_mf[lev]->setVal(1.0);

	u_mf[lev].get()->setVal(0.0);
	
}


void Hydro::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
  std::swap(rho_old_mf[lev], rho_mf[lev]);
  std::swap(u_old_mf[lev], u_mf[lev]);
  std::swap(p_old_mf[lev], p_mf[lev]);
  std::swap(e_old_mf[lev], e_mf[lev]);

  static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
				     dy(AMREX_D_DECL(0,1,0)),
				     dz(AMREX_D_DECL(0,0,1)));

  const amrex::Real *DX = geom[lev].CellSize();

  for (amrex::MFIter mfi(*rho_mf[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::Array4<Set::Scalar> const &rho		        = (*rho_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &rho_old   	= (*rho_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &u		        = (*u_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &u_old   	        = (*u_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &p		        = (*p_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &p_old		= (*p_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &e		        = (*e_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &e_old	        = (*e_old_mf[lev]).array(mfi);


		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
			{ 
				// e(i,j,k) value at i,j,k
				// e(i,j,k,n)  value at i,j,k for component n

				Set::Vector rho_grad = Numeric::Gradient(rho_old, i, j, k, 0, DX);
				Set::Vector p_grad = Numeric::Gradient(p_old, i, j, k, 0, DX);
				Set::Vector e_grad = Numeric::Gradient(e_old, i, j, k, 0, DX);

				if (j == bx.loVect()[1] || j == bx.hiVect()[1] || k == bx.loVect()[2] || k == bx.hiVect()[2])
				  {
				    Set::Vector u1_grad(1.0, 0.0, 0.0);
				    Set::Vector u2_grad = Set::Vector::Zero();
				    Set::Vector u3_grad = Set::Vector::Zero();

				    rho(i, j, k) = rho_old(i, j, k) - (u_old(i, j, k, 0)*(rho_grad(0)) + u_old(i, j, k, 1)*(rho_grad(1)) + u_old(i, j, k, 2)*(rho_grad(2)) + rho_old(i, j, k)*(u1_grad(0) + u2_grad(1) + u3_grad(2)))*dt;

				    if(i == bx.loVect()[0])
				      {u(i, j, k, 0) = 0.0;
				      u(i, j, k, 1) = 0.0;
				      u(i, j, k, 2) = 0.0;}

				    else if(i == bx.hiVect()[0])
				      {u(i, j, k, 0) = 1.0;
				      u(i, j, k, 1) = 0.0;
				      u(i, j, k, 2) = 0.0;}

				    else {
				      u(i, j, k, 0) = u_old(i, j, k, 0) - (u_old(i, j, k, 0)*(u_old(i, j, k, 0)*(u1_grad(0)) + u_old(i, j, k, 1)*(u1_grad(1)) + u_old(i, j, k, 2)*(u1_grad(2))) + 1/rho_old(i, j, k)*(p_grad(0)))*dt;
				
				      u(i, j, k, 1) = u_old(i, j, k, 1) - (u_old(i, j, k, 1)*(u_old(i, j, k, 0)*(u2_grad(0)) + u_old(i, j, k, 1)*(u2_grad(1)) + u_old(i, j, k, 2)*(u2_grad(2))) + 1/rho_old(i, j, k)*(p_grad(1)))*dt;
				
				      u(i, j, k, 2) = u_old(i, j, k, 2) - (u_old(i, j, k, 2)*(u_old(i, j, k, 0)*(u3_grad(0)) + u_old(i, j, k, 1)*(u3_grad(1)) + u_old(i, j, k, 2)*(u3_grad(2))) + 1/rho_old(i, j, k)*(p_grad(2)))*dt;
		       }
				
				    e(i,j,k) = e_old(i,j,k) - ((u_old(i, j, k, 0)*(e_grad(0)) + u_old(i, j, k, 1)*(e_grad(1)) + u_old(i, j, k, 2)*(e_grad(2))) + p_old(i, j, k)/rho_old(i, j, k)*(u1_grad(0) + u2_grad(1) + u3_grad(2)))*dt;
				
				    p(i, j, k) = (gamma - 1)*rho(i, j, k)*e(i, j, k);
				  }

				else
				  {
				    Set::Vector u1_grad = Numeric::Gradient(u_old, i, j, k, 0, DX);
				    Set::Vector u2_grad = Numeric::Gradient(u_old, i, j, k, 0, DX);
				    Set::Vector u3_grad = Numeric::Gradient(u_old, i, j, k, 0, DX);

				    rho(i, j, k) = rho_old(i, j, k) - (u_old(i, j, k, 0)*(rho_grad(0)) + u_old(i, j, k, 1)*(rho_grad(1)) + u_old(i, j, k, 2)*(rho_grad(2)) + rho_old(i, j, k)*(u1_grad(0) + u2_grad(1) + u3_grad(2)))*dt;

				    if(i == bx.loVect()[0])
				      {u(i, j, k, 0) = 0.0;
				      u(i, j, k, 1) = 0.0;
				      u(i, j, k, 2) = 0.0;}

				    else if(i == bx.hiVect()[0])
				      {u(i, j, k, 0) = 1.0;
				      u(i, j, k, 1) = 0.0;
				      u(i, j, k, 2) = 0.0;}

				    else {
				      u(i, j, k, 0) = u_old(i, j, k, 0) - (u_old(i, j, k, 0)*(u_old(i, j, k, 0)*(u1_grad(0)) + u_old(i, j, k, 1)*(u1_grad(1)) + u_old(i, j, k, 2)*(u1_grad(2))) + 1/rho_old(i, j, k)*(p_grad(0)))*dt;
				
				      u(i, j, k, 1) = u_old(i, j, k, 1) - (u_old(i, j, k, 1)*(u_old(i, j, k, 0)*(u2_grad(0)) + u_old(i, j, k, 1)*(u2_grad(1)) + u_old(i, j, k, 2)*(u2_grad(2))) + 1/rho_old(i, j, k)*(p_grad(1)))*dt;
				
				      u(i, j, k, 2) = u_old(i, j, k, 2) - (u_old(i, j, k, 2)*(u_old(i, j, k, 0)*(u3_grad(0)) + u_old(i, j, k, 1)*(u3_grad(1)) + u_old(i, j, k, 2)*(u3_grad(2))) + 1/rho_old(i, j, k)*(p_grad(2)))*dt;
				    }
				
				    e(i,j,k) = e_old(i,j,k) - ((u_old(i, j, k, 0)*(e_grad(0)) + u_old(i, j, k, 1)*(e_grad(1)) + u_old(i, j, k, 2)*(e_grad(2))) + p_old(i, j, k)/rho_old(i, j, k)*(u1_grad(0) + u2_grad(1) + u3_grad(2)))*dt;
				
				    p(i, j, k) = (gamma - 1)*rho(i, j, k)*e(i, j, k);
				  }
				

				if (std::isnan(rho(amrex::IntVect(AMREX_D_DECL(i,j,k)))))
					Util::Abort(INFO, "NaN encountered");
		  });
    }
}

void Hydro::TagCellsForRefinement (int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
		const amrex::Real *DX = geom[lev].CellSize();
		Set::Scalar dr  = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

		for (amrex::MFIter mfi(*rho_mf[lev], true); mfi.isValid(); ++mfi)
		{
			const amrex::Box &bx = mfi.tilebox();
			amrex::Array4<char>              const &tags = a_tags.array(mfi);
			amrex::Array4<const Set::Scalar> const &rho  = (*rho_mf[lev]).array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				Set::Vector gradrho = Numeric::Gradient(rho,i,j,k,0,DX);
				if (gradrho.lpNorm<2>() * dr > 0.001)
					tags(i,j,k) = amrex::TagBox::SET;
			});
		}
	}
}
