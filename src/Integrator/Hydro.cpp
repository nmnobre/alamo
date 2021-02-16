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

  rhoIC = new IC::Constant(geom);
  pIC   = new IC::Constant(geom);
  eIC   = new IC::Constant(geom);

  RegisterNewFab(rho_mf,     rhoBC, 1, 0, "rho", true);
  RegisterNewFab(rho_old_mf, rhoBC, 1, 0, "rho_old", false);
  RegisterNodalFab(u1_mf, 1, 0, "u1", true);
  RegisterNodalFab(u1_old_mf, 1, 0, "u1_old", false);
  RegisterNodalFab(u2_mf, 1, 0, "u2", true);
  RegisterNodalFab(u2_old_mf, 1, 0, "u2_old", false);
  RegisterNodalFab(u3_mf, 1, 0, "u3", true);
  RegisterNodalFab(u3_old_mf, 1, 0, "u3_old", false);
  RegisterNewFab(p_mf,     pBC, 1, 0, "p", true);
  RegisterNewFab(p_old_mf, pBC, 1, 0, "p_old", false);
  RegisterNewFab(e_mf,     eBC, 1, 0, "e", true);
  RegisterNewFab(e_old_mf, eBC, 1, 0, "e_old", false);
}

void Hydro::Initialize(int lev)
{        
        rho_mf[lev]->setVal(1.0);
	rho_old_mf[lev]->setVal(1.0);
	
	u1_mf[lev].get()->setVal(1.0);
	u2_mf[lev].get()->setVal(0.0);
	u3_mf[lev].get()->setVal(0.0);

	pIC->Initialize(lev, p_mf);
	pIC->Initialize(lev, p_old_mf);

	eIC->Initialize(lev, e_mf);
	eIC->Initialize(lev, e_old_mf);

	// u1BC = BC::Operator::Elastic::Constant::Type::Displacement;
	// u2BC = BC::Operator::Elastic::Constant::Type::Displacement;
	// u3BC = BC::Operator::Elastic::Constant::Type::Displacement;

}


void Hydro::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
  std::swap(rho_old_mf[lev], rho_mf[lev]);
  std::swap(u1_old_mf[lev], u1_mf[lev]);
  std::swap(u2_old_mf[lev], u2_mf[lev]);
  std::swap(u3_old_mf[lev], u3_mf[lev]);
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
      amrex::Array4<Set::Scalar> const &u1		        = (*u1_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &u1_old    	= (*u1_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &u2		        = (*u2_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &u2_old	        = (*u2_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &u3		        = (*u3_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &u3_old    	= (*u3_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &p		        = (*p_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &p_old		= (*p_old_mf[lev]).array(mfi);
      amrex::Array4<Set::Scalar> const &e		        = (*e_mf[lev]).array(mfi);
      amrex::Array4<const Set::Scalar> const &e_old	        = (*e_old_mf[lev]).array(mfi);


		AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
						 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
						 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));

				// e(i,j,k) value at i,j,k
				// e(i,j,k,n)  value at i,j,k for component n

				AMREX_D_TERM(amrex::Real rho_gradx = (rho_old(m+dx) - rho_old(m-dx))/(2*DX[0]);,
								 amrex::Real rho_grady = (rho_old(m+dy) - rho_old(m-dy))/(2*DX[1]);,
								 amrex::Real rho_gradz = (rho_old(m+dz) - rho_old(m-dz))/(2*DX[2]););
				
				AMREX_D_TERM(amrex::Real u1_gradx = (u1_old(m+dx) - u1_old(m-dx))/(2*DX[0]);,
								 amrex::Real u1_grady = (u1_old(m+dy) - u1_old(m-dy))/(2*DX[1]);,
								 amrex::Real u1_gradz = (u1_old(m+dz) - u1_old(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real u2_gradx = (u2_old(m+dx) - u2_old(m-dx))/(2*DX[0]);,
								 amrex::Real u2_grady = (u2_old(m+dy) - u2_old(m-dy))/(2*DX[1]);,
								 amrex::Real u2_gradz = (u2_old(m+dz) - u2_old(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real u3_gradx = (u3_old(m+dx) - u3_old(m-dx))/(2*DX[0]);,
								 amrex::Real u3_grady = (u3_old(m+dy) - u3_old(m-dy))/(2*DX[1]);,
								 amrex::Real u3_gradz = (u3_old(m+dz) - u3_old(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real p_gradx = (p_old(m+dx) - p_old(m-dx))/(2*DX[0]);,
								 amrex::Real p_grady = (p_old(m+dy) - p_old(m-dy))/(2*DX[1]);,
								 amrex::Real p_gradz = (p_old(m+dz) - p_old(m-dz))/(2*DX[2]););
				
				AMREX_D_TERM(amrex::Real e_gradx = (e_old(m+dx) - e_old(m-dx))/(2*DX[0]);,
								 amrex::Real e_grady = (e_old(m+dy) - e_old(m-dy))/(2*DX[1]);,
								 amrex::Real e_gradz = (e_old(m+dz) - e_old(m-dz))/(2*DX[2]););

				rho(m) = rho_old(m) - (u1_old(m)*(rho_gradx) + u2_old(m)*(rho_grady) + u3_old(m)*(rho_gradz) + rho_old(m)*(u1_gradx + u2_grady + u3_gradz))*dt;
				
				u1(m) = u1_old(m) - (u1_old(m)*(u1_old(m)*(u1_gradx) + u2_old(m)*(u1_grady) + u3_old(m)*(u1_gradz)) + 1/rho_old(m)*(p_gradx))*dt;
				
				u2(m) = u2_old(m) - (u2_old(m)*(u1_old(m)*(u2_gradx) + u2_old(m)*(u2_grady) + u3_old(m)*(u2_gradz)) + 1/rho_old(m)*(p_grady))*dt;
				
				u3(m) = u3_old(m) - (u3_old(m)*(u1_old(m)*(u3_gradx) + u2_old(m)*(u3_grady) + u3_old(m)*(u3_gradz)) + 1/rho_old(m)*(p_gradz))*dt;
				
				e(m) = e_old(m) - ((u1_old(m)*(e_gradx) + u2_old(m)*(e_grady) + u3_old(m)*(e_gradz)) + p_old(m)/rho_old(m)*(u1_gradx + u2_grady + u3_gradz))*dt;
				
				p(m) = (gamma - 1)*rho(m)*e(m);

				if (std::isnan(rho(amrex::IntVect(AMREX_D_DECL(i,j,k)))))
					Util::Abort(INFO, "NaN encountered");
			}
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
