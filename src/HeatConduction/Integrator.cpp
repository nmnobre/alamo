#include "Integrator.H"

/// \fn    HeatConduction::Integrator::Integrator
///
/// Read in the following simulation parameters
///
///     heat.alpha                (default 1.0)
///     heat.refinement_threshold (default 0.01)
///     ic.type
///
/// Initialize initial condition pointer #ic, and register
/// the #Temp, #Temp_old Multifab arrays.

Integrator::HeatConduction::HeatConduction() :
  Integrator::Integrator()
{
  amrex::ParmParse pp("heat");
  pp.query("alpha", alpha);
  pp.query("refinement_threshold", refinement_threshold);
  pp.query("ic_type", ic_type);

  // // Determine initial condition
  if (ic_type == "cylinder")
    ic = new IC::Cylinder(geom);
  else if (ic_type == "constant")
    ic = new IC::Constant(geom);
  else
    ic = new IC::Constant(geom);
    
  {
    amrex::ParmParse pp("bc");
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

    mybc = new BC::BC(geom,
		      bc_hi_str, bc_lo_str,
		      bc_lo_1, bc_hi_1,
		      bc_lo_2, bc_hi_2,
		      bc_lo_3, bc_hi_3);
  }


  RegisterNewFab(Temp,     *mybc, number_of_components, number_of_ghost_cells, "Temp");
  RegisterNewFab(Temp_old, *mybc, number_of_components, number_of_ghost_cells, "Temp old");
}

Integrator::HeatConduction::~HeatConduction()
{
}

/// \fn HeatConduction::Integrator::Initialize
///
/// Use the #ic object to initialize #Temp
void
Integrator::HeatConduction::Initialize (int lev)
{
  ic->Initialize(lev,Temp);
}


/// \fn    Integrator::HeatConduction::Advance
///
/// Integrate the heat diffusion equation
/// \f[\nabla^2T = \alpha \frac{\partial T}{\partial t}\f]
/// using an explicit forward Euler method.
/// \f$\alpha\f$ is stored in #alpha
void
Integrator::HeatConduction::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
  std::swap(*Temp[lev], *Temp_old[lev]);

  const amrex::Real* dx = geom[lev].CellSize();

  for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
    {
      const amrex::Box& bx = mfi.tilebox();

      amrex::FArrayBox &Temp_old_box = (*Temp_old[lev])[mfi];
      amrex::FArrayBox &Temp_box = (*Temp[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
	for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
	  {
	    TEMP(i,j,k)
	      = TEMP_OLD(i,j,k)
	      + dt * alpha * (AMREX_D_TERM((TEMP_OLD(i+1,j,k) + TEMP_OLD(i-1,j,k) - 2*TEMP_OLD(i,j,k)) / dx[0] / dx[0],
					   + (TEMP_OLD(i,j+1,k) + TEMP_OLD(i,j-1,k) - 2*TEMP_OLD(i,j,k)) / dx[1] / dx[1],
					   + (TEMP_OLD(i,j,k+1) + TEMP_OLD(i,j,k-1) - 2*TEMP_OLD(i,j,k)) / dx[2] / dx[2]));
		      }
    }
}


/// \fn    Integrator::HeatConduction::TagCellsForRefinement
///
/// The following criterion is used to determine if a cell should be refined:
/// \f[|\nabla T|\,|\mathbf{r}| > h\f]
/// where
/// \f[\mathbf{r} = \sqrt{\Delta x_1^2 + \Delta x_2^2 + \Delta x_3^2}\f]
/// and \f$h\f$ is stored in #refinement_threshold
void
Integrator::HeatConduction::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
  const amrex::Real* dx      = geom[lev].CellSize();

  //amrex::Array<int>  itags;
  for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
    {
      const amrex::Box&  bx  = mfi.tilebox();
      amrex::TagBox&     tag  = tags[mfi];
 	    
      amrex::FArrayBox &Temp_box = (*Temp[lev])[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
	for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
	  {
	    amrex::Real grad1 = (TEMP(i+1,j,k) - TEMP(i-1,j,k));
	    amrex::Real grad2 = (TEMP(i,j+1,k) - TEMP(i,j-1,k));
#if AMREX_SPACEDIM>2
 	    amrex::Real grad3 = (TEMP(i,j,k+1) - TEMP(i,j,k-1));
#endif
	    amrex::Real grad = sqrt(AMREX_D_TERM(grad1*grad1,
						 + grad2*grad2,
						 + grad3*grad3));

	    amrex::Real dr = sqrt(AMREX_D_TERM(dx[0]*dx[0],
					       + dx[1]*dx[1],
					       + dx[2]*dx[2]));

	    if (grad*dr > refinement_threshold)
	      tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
	  }
    }
}

