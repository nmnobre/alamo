#include "Hydro.H"
#include "IC/Constant.H"

namespace Integrator
{
Hydro::Hydro () : Integrator()
{
  amrex::ParmParse pp("physics"); 
  pp.query("gamma",gamma);

  {
    amrex::ParmParse pp("TempBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    rhoBC = new BC::Constant(1, bc_hi_str, bc_lo_str,
			      AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3),
			      AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    u1BC = new BC::Constant(1,bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    u2BC = new BC::Constant(1,bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    u3BC = new BC::Constant(1,bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    pBC = new BC::Constant(1,bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }
  {
    amrex::ParmParse pp("EtaBC");
    amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
    amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
    pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
    pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
    amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
    if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
    if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
    amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
    if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
    if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
    amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;
    if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
    if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

    eBC = new BC::Constant(1,bc_hi_str, bc_lo_str,
			     AMREX_D_DECL(bc_lo_1,bc_lo_2,bc_lo_3),
			     AMREX_D_DECL(bc_hi_1,bc_hi_2,bc_hi_3));
  }

  rhoIC = new IC::Constant(geom);
  u1IC = new IC::Constant(geom);
  u2IC = new IC::Constant(geom);
  u3IC = new IC::Constant(geom);
  pIC = new IC::Constant(geom);
  eIC = new IC::Constant(geom);

  RegisterNewFab(rho,     rhoBC, 1, 0, "rho", true);
  RegisterNewFab(rho_old, rhoBC, 1, 0, "rho_old", false);
  RegisterNodalFab(u1,     u1BC, 1, 0, "u1", true);
  RegisterNodalFab(u1_old, u1BC, 1, 0, "u1_old", false);
  RegisterNodalFab(u2,     u2BC, 1, 0, "u2", true);
  RegisterNodalFab(u2_old, u2BC, 1, 0, "u2_old", false);
  RegisterNodalFab(u3,     u3BC, 1, 0, "u3", true);
  RegisterNodalFab(u3_old, u3BC, 1, 0, "u3_old", false);
  RegisterNewFab(p,     pBC, 1, 0, "p", true);
  RegisterNewFab(p_old, pBC, 1, 0, "p_old", false);
  RegisterNewFab(e,     eBC, 1, 0, "e", true);
  RegisterNewFab(e_old, eBC, 1, 0, "e_old", false);
}

void Hydro::Initialize (int lev)
{
	for (amrex::MFIter mfi(*rho[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &rho_box		= (*rho[lev])[mfi];
			amrex::BaseFab<amrex::Real> &rho_old_box	= (*rho_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u1_box		= (*u1[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u1_old_box	        = (*u1_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u2_box		= (*u2[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u2_old_box	        = (*u2_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u3_box		= (*u3[lev])[mfi];
			amrex::BaseFab<amrex::Real> &u3_old_box	        = (*u3_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &p_box		= (*p[lev])[mfi];
			amrex::BaseFab<amrex::Real> &p_old_box	        = (*p_old[lev])[mfi];
			amrex::BaseFab<amrex::Real> &e_box		= (*e[lev])[mfi];
			amrex::BaseFab<amrex::Real> &e_old_box	        = (*e_old[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
							 for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
							 for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
				{
					rho_box     (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					rho_old_box (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u1_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u1_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u2_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u2_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u3_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					u3_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					p_box     (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					p_old_box (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					e_box    (amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
					e_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) =  0.0;
				}
		}
	rhoIC->Initialize(lev,rho);
	rhoIC->Initialize(lev,rho_old);
	
	u1IC->Initialize(lev,u1);
	u1IC->Initialize(lev,u1_old);

	u2IC->Initialize(lev,u2);
	u2IC->Initialize(lev,u2_old);

	u3IC->Initialize(lev,u3);
	u3IC->Initialize(lev,u3_old);

	pIC->Initialize(lev,p);
	pIC->Initialize(lev,p_old);

	eIC->Initialize(lev,e);
	eIC->Initialize(lev,e_old);
}



void Hydro::Advance (int lev, amrex::Real time, amrex::Real dt)
{
  std::swap(rho_old[lev], rho[lev]);
  std::swap(u1_old[lev], u1[lev]);
  std::swap(u2_old[lev], u2[lev]);
  std::swap(u3_old[lev], u3[lev]);
  std::swap(p_old[lev], p[lev]);
  std::swap(e_old[lev], e[lev]);

  static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
												 dy(AMREX_D_DECL(0,1,0)),
												 dz(AMREX_D_DECL(0,0,1)));

  const amrex::Real* DX = geom[lev].CellSize();

  for ( amrex::MFIter mfi(*rho[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::FArrayBox &rho_box		        = (*rho[lev])[mfi];
      amrex::FArrayBox &rho_old_box		= (*rho_old[lev])[mfi];
      amrex::FArrayBox &u1_box		        = (*u1[lev])[mfi];
      amrex::FArrayBox &u1_old_box	        = (*u1_old[lev])[mfi];
      amrex::FArrayBox &u2_box		        = (*u2[lev])[mfi];
      amrex::FArrayBox &u2_old_box	        = (*u2_old[lev])[mfi];
      amrex::FArrayBox &u3_box		        = (*u3[lev])[mfi];
      amrex::FArrayBox &u3_old_box	        = (*u3_old[lev])[mfi];
      amrex::FArrayBox &p_box		        = (*p[lev])[mfi];
      amrex::FArrayBox &p_old_box		= (*p_old[lev])[mfi];
      amrex::FArrayBox &e_box		        = (*e[lev])[mfi];
      amrex::FArrayBox &e_old_box	        = (*e_old[lev])[mfi];

		AMREX_D_TERM(for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
						 for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
						 for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++))
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));

				AMREX_D_TERM(amrex::Real rho_gradx = (rho_old_box(m+dx) - rho_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real rho_grady = (rho_old_box(m+dy) - rho_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real rho_gradz = (rho_old_box(m+dz) - rho_old_box(m-dz))/(2*DX[2]););
				
				AMREX_D_TERM(amrex::Real u1_gradx = (u1_old_box(m+dx) - u1_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real u1_grady = (u1_old_box(m+dy) - u1_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real u1_gradz = (u1_old_box(m+dz) - u1_old_box(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real u2_gradx = (u2_old_box(m+dx) - u2_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real u2_grady = (u2_old_box(m+dy) - u2_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real u2_gradz = (u2_old_box(m+dz) - u2_old_box(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real u3_gradx = (u3_old_box(m+dx) - u3_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real u3_grady = (u3_old_box(m+dy) - u3_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real u3_gradz = (u3_old_box(m+dz) - u3_old_box(m-dz))/(2*DX[2]););

				AMREX_D_TERM(amrex::Real p_gradx = (p_old_box(m+dx) - p_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real p_grady = (p_old_box(m+dy) - p_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real p_gradz = (p_old_box(m+dz) - p_old_box(m-dz))/(2*DX[2]););
				
				AMREX_D_TERM(amrex::Real e_gradx = (e_old_box(m+dx) - e_old_box(m-dx))/(2*DX[0]);,
								 amrex::Real e_grady = (e_old_box(m+dy) - e_old_box(m-dy))/(2*DX[1]);,
								 amrex::Real e_gradz = (e_old_box(m+dz) - e_old_box(m-dz))/(2*DX[2]););

				rho_box(m) = rho_old_box(m) - (u1_old_box(m)*(rho_gradx) + u2_old_box(m)*(rho_grady) + u3_old_box(m)*(rho_gradz) + rho_old_box(m)*(u1_gradx + u2_grady + u3_gradz))*dt;
				
				u1_box(m) = u1_old_box(m) - (u1_old_box(m)*(u1_old_box(m)*(u1_gradx) + u2_old_box(m)*(u1_grady) + u3_old_box(m)*(u1_gradz)) + 1/rho_old_box(m)*(p_gradx))*dt;
				
				u2_box(m) = u2_old_box(m) - (u2_old_box(m)*(u1_old_box(m)*(u2_gradx) + u2_old_box(m)*(u2_grady) + u3_old_box(m)*(u2_gradz)) + 1/rho_old_box(m)*(p_grady))*dt;
				
				u3_box(m) = u3_old_box(m) - (u3_old_box(m)*(u1_old_box(m)*(u3_gradx) + u2_old_box(m)*(u3_grady) + u3_old_box(m)*(u3_gradz)) + 1/rho_old_box(m)*(p_gradz))*dt;
				
				e_box(m) = e_old_box(m) - ((u1_old_box(m)*(e_gradx) + u2_old_box(m)*(e_grady) + u3_old_box(m)*(e_gradz)) + p_old_box(m)/rho_old_box(m)*(u1_gradx + u2_grady + u3_gradz))*dt;
				
				p_box(m) = (gamma - 1)*rho_box(m)*e_box(m);

				if (std::isnan(rho_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))))
					Util::Abort(INFO, "NaN encountered");
			}
    }
}



void Hydro::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{

	const amrex::Real* dx      = geom[lev].CellSize();

	const amrex::MultiFab& state = *rho[lev];
	
	for (amrex::MFIter mfi(state,true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();

			amrex::TagBox&     tag  = tags[mfi];
	    
			amrex::BaseFab<amrex::Real> &rho_box = (*rho[lev])[mfi];

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if BL_SPACEDIM>2
					for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
						{
							amrex::Real gradx = (rho_box(amrex::IntVect(AMREX_D_DECL(i+1,j,k))) - rho_box(amrex::IntVect(AMREX_D_DECL(i-1,j,k))))/(2.*dx[0]);
							amrex::Real grady = (rho_box(amrex::IntVect(AMREX_D_DECL(i,j+1,k))) - rho_box(amrex::IntVect(AMREX_D_DECL(i,j-1,k))))/(2.*dx[1]);
							amrex::Real gradz = (rho_box(amrex::IntVect(AMREX_D_DECL(i,j,k+1))) - rho_box(amrex::IntVect(AMREX_D_DECL(i,j,k-1))))/(2.*dx[2]);
							if (dx[0]*dx[1]*dx[2]*(gradx*gradx + grady*grady + gradz*gradz)>0.001) tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
						}
		}

}
}
