#ifndef INTEGRATOR_HYDRO_H
#define INTEGRATOR_HYDRO_H

#include <string>
#include <limits>
#include <memory>

#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "BC/Constant.H"
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
	void Regrid(int lev, Set::Scalar time) override;
private:

	Set::Field<Set::Scalar> rho;
	Set::Field<Set::Scalar> rho_old;
	Set::Field<Set::Scalar> u1;
	Set::Field<Set::Scalar> u1_old;
        Set::Field<Set::Scalar> u2;
	Set::Field<Set::Scalar> u2_old;
        Set::Field<Set::Scalar> u3;
	Set::Field<Set::Scalar> u3_old;
	Set::Field<Set::Scalar> p;
        Set::Field<Set::Scalar> p_old;
	Set::Field<Set::Scalar> e;
	Set::Field<Set::Scalar> e_old;

        BC::BC<Set::Scalar> *rhoBC;
	BC::BC<Set::Scalar> *u1BC;
        BC::BC<Set::Scalar> *u2BC;
	BC::BC<Set::Scalar> *u3BC;
        BC::BC<Set::Scalar> *pBC;
	BC::BC<Set::Scalar> *eBC;
 
	IC::Constant *rhoIC;
	IC::Constant *u1IC;
        IC::Constant *u2IC;
        IC::Constant *u3IC;
	IC::Constant *pIC;
        IC::Constant *eIC;

	amrex::Real gamma;                  // adiabatic gamma

};
}

#endif
