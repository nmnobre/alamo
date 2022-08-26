#ifndef INTEGRATOR_MYFLAME_H
#define INTEGRATOR_MYFLAME_H

#include <string>
#include <limits>
#include <memory>


#include "Util/Util.H"
#include "Integrator/Integrator.H"
#include "Integrator/Base/Mechanics.H"
#include "BC/BC.H"
#include "IC/Voronoi.H"
#include "IC/IC.H"
#include "IC/Wedge.H"
#include "IC/Constant.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/Isotropic.H"
#include "IO/ParmParse.H"
#include "BC/Operator/Elastic/Constant.H"
#include "Solver/Nonlocal/Newton.H"



namespace Integrator
{

class MyFlame : virtual public Integrator::Integrator, 
                virtual public Base::Mechanics<Model::Solid::Affine::Isotropic>
{

    using model_type = Model::Solid::Affine::Isotropic;

public:
    MyFlame () : Base::Mechanics<Model::Solid::Affine::Isotropic>() {};
    MyFlame (IO::ParmParse &pp) : Base::Mechanics<Model::Solid::Affine::Isotropic>()
    {pp.queryclass(*this);}

    static void Parse(MyFlame &value, IO::ParmParse &pp)
    {
        value.bc_eta = new BC::Constant(1,pp,"eta.bc.constant"); 
        value.ic_eta = new IC::Constant(value.geom,pp,"eta.ic.constant");
        value.RegisterNewFab(value.eta_mf,value.bc_eta,1,2,"eta",true);
        value.RegisterNewFab(value.eta_old_mf,value.bc_eta,1,2,"eta_old",false);
        value.m_type = Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable;
    };

protected:
    void Initialize (int lev) override
    {
        Base::Mechanics<Model::Solid::Affine::Isotropic>::Initialize(lev);
        // TODO
    }
    void TimeStepBegin(Set::Scalar a_time, int a_iter) override
    {
        Base::Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time,a_iter);
        // TODO
    }
    void Advance (int lev, Set::Scalar time, Set::Scalar dt) override
    {
        Base::Mechanics<Model::Solid::Affine::Isotropic>::Advance(lev,time,dt);
        // TODO
    }
    void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override
    {
        Base::Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev,tags,time,ngrow);
        // TODO
    }
    void UpdateModel(int /*a_step*/) override
    {
        //TODO
    }
private:

    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> eta_old_mf;
    BC::BC<Set::Scalar> *bc_eta;
    IC::IC *ic_eta;
};
}

#endif
