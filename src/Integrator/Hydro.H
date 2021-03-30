#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include <string>
#include <limits>
#include <memory>

#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/Constant.H"
#include "BC/BC.H"
#include "IC/Constant.H"
#include "IC/IC.H"

namespace Integrator
{
class Hydro : public Integrator::Integrator
{
public:
	Hydro ();

protected:
	void Initialize (int lev) ;
	void Advance (int lev, amrex::Real time, amrex::Real dt);
	void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);
	//void Regrid(int lev, Set::Scalar time) override;
private:

	Set::Field<Set::Scalar> rho_mf;
	Set::Field<Set::Scalar> rho_old_mf;
        Set::Field<Set::Scalar> u_mf;
	Set::Field<Set::Scalar> u_old_mf;
	// Set::Field<Set::Scalar> u1_mf;
	// Set::Field<Set::Scalar> u1_old_mf;
        // Set::Field<Set::Scalar> u2_mf;
	// Set::Field<Set::Scalar> u2_old_mf;
        // Set::Field<Set::Scalar> u3_mf;
	// Set::Field<Set::Scalar> u3_old_mf;
	Set::Field<Set::Scalar> p_mf;
        Set::Field<Set::Scalar> p_old_mf;
	Set::Field<Set::Scalar> e_mf;
	Set::Field<Set::Scalar> e_old_mf;

        BC::BC<Set::Scalar> *rhoBC;
	BC::BC<Set::Scalar> *u1BC;
        BC::BC<Set::Scalar> *u2BC;
	BC::BC<Set::Scalar> *u3BC;
        BC::BC<Set::Scalar> *pBC;
	BC::BC<Set::Scalar> *eBC;
 
	// IC::IC *rhoIC;
	// IC::IC *pIC;
        // IC::IC *eIC;
  // initial velocities are hard coded

	amrex::Real gamma;                  // adiabatic gamma

};
}

#endif
