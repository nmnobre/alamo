#ifndef INTEGRATOR_DYNAMICS_H
#define INTEGRATOR_DYNAMICS_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Affine/Isotropic.H"

#include "Operator/Operator.H"

#include "IO/ParmParse.H"

namespace Integrator
{
class Dynamics : public Integrator
{
    using model_type = Model::Solid::Affine::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    Dynamics()
    {
        RegisterNodalFab(u_mf,     AMREX_SPACEDIM, 2, "u",true);
        RegisterNodalFab(u_old_mf, AMREX_SPACEDIM, 2, "uold",false);
        RegisterNodalFab(v_mf,     AMREX_SPACEDIM, 2, "v",true);
        RegisterNodalFab(v_old_mf, AMREX_SPACEDIM, 2, "vold",false);
        RegisterNodalFab(trace_eps_mf, AMREX_SPACEDIM, 1, "treps",true);
        RegisterGeneralFab(sig_mf, 1, 2);

        RegisterNodalFab(eta_mf, 1, 2, "eta",true);
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                ic = new IC::Ellipse(geom);
                pp.queryclass("ellipse",static_cast<IC::Ellipse*>(ic));
            }
            if (pp.contains("Finclusion")) pp.queryarr("Finclusion",Finclusion);
            if (pp.contains("Fmatrix"))    pp.queryarr("Fmatrix",Fmatrix);
        }
        {
            IO::ParmParse pp("elastic");
            {
                // Initialize elastic model
                pp.queryclass("model",elastic.model);
                elastic.model.F0 = Set::Matrix::Zero();

                // TODO: figure out how to initialize solvers in
                // here instead of TimeStepBegin...

                pp.queryclass("bc",elastic.bc);

                pp.query("ref_criterion",elastic.ref_criterion);
            }
        }
        
        RegisterGeneralFab(model_mf, 1, 2);
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        eta_mf[lev]->setVal(0.0);
        ic->Initialize(lev, eta_mf);
        eta_mf[lev]->FillBoundary();

        disp_mf[lev]->setVal(0.);
        rhs_mf[lev]->setVal(0.);

        u_mf[lev]->setVal(0.0);
        u_old_mf[lev]->setVal(0.0);
        v_mf[lev]->setVal(0.0);
        v_old_mf[lev]->setVal(0.0);

        model_mf[lev]->setVal(elastic.model);

        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.growntilebox(2);

            amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Matrix Fgb = (1.0 - eta(i, j, k)) * Fmatrix + eta(i, j, k) * Finclusion;
                model(i, j, k) = elastic.model;
                model(i, j, k).F0 = Fgb;
            });
        }

        Util::RealFillBoundary(*model_mf[lev],geom[lev]);

    }

    void TimeStepBegin(Set::Scalar, int) override
    {
        if (finest_level != rhs_mf.size() - 1)
        {
            Util::Abort(INFO, "amr.max_level is larger than necessary. Set to ", finest_level, " or less");
        }
        
        // Set linear elastic model
        
        for (int lev = 0; lev < rhs_mf.size(); ++lev)
        {
        }

//        elastic.bc.Init(rhs_mf,geom);
//
//        amrex::LPInfo info;
//        Operator::Elastic<Model::Solid::Affine::Isotropic> elastic_op(geom, grids, dmap, info);
//        elastic_op.SetUniform(false);
//        elastic_op.SetBC(&elastic.bc);
//
//        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
//
//        IO::ParmParse pp("elastic");
//        elastic.solver = new Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic>(elastic_op);
//        pp.queryclass("solver",*elastic.solver);
//
//        elastic.solver->solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);
//
//        for (int lev = 0; lev < disp_mf.size(); lev++)
//        {
//            elastic_op.Stress(lev, *stress_mf[lev], *disp_mf[lev]);
//            elastic_op.Strain(lev, *strain_mf[lev], *disp_mf[lev]);
//        }
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
    	BL_PROFILE("PhaseFieldMicrostructure::Advance");
    	/// TODO Make this optional
    	//if (lev != max_level) return;
    	std::swap(u_old_mf[lev], u_mf[lev]);
    	std::swap(v_old_mf[lev], v_mf[lev]);
    	const amrex::Real *DX = geom[lev].CellSize();
        const Set::Scalar vol = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);
        amrex::Box domain(geom[lev].Domain());
    
    	for (amrex::MFIter mfi(*u_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    	{
            Box bx = mfi.tilebox();
		    bx.grow(1);
            Box dom = domain;
		    bx = bx & dom;

    		amrex::Array4<amrex::Real> const &unew = (*u_mf[lev]).array(mfi);
    		amrex::Array4<amrex::Real> const &vnew = (*v_mf[lev]).array(mfi);
    		amrex::Array4<const amrex::Real> const &u = (*u_mf[lev]).array(mfi);
    		amrex::Array4<const amrex::Real> const &v = (*v_mf[lev]).array(mfi);
            amrex::Array4<Set::Matrix>       const &sig = (*sig_mf[lev]).array(mfi);
            amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);         
    		amrex::Array4<amrex::Real> const &treps = (*trace_eps_mf[lev]).array(mfi);
    
    		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                std::array<Numeric::StencilType,AMREX_SPACEDIM>
					sten = Numeric::GetStencil(i,j,k,domain);

                Set::Matrix gradu = Numeric::Gradient(u,i,j,k,DX,sten);
                sig(i,j,k) = model(i,j,k).DW(gradu);
                treps(i,j,k) = gradu.trace();
            });

            dom.grow(-1);
            bx = bx & dom;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

                std::array<Numeric::StencilType,AMREX_SPACEDIM>
					sten = Numeric::GetStencil(i,j,k,domain);

                Set::Matrix sig_x = ((sig(i+1,j,k) - sig(i-1,j,k))/DX[0]);
                Set::Matrix sig_y = ((sig(i,j+1,k) - sig(i,j-1,k))/DX[1]);

                Set::Vector f;
                f(0) = sig_x(0,0) + sig_y(0,1);
                f(1) = sig_x(1,0) + sig_y(1,1);

                vnew(i,j,k,0) = v(i,j,k,0) + dt*f(0);
                vnew(i,j,k,1) = v(i,j,k,1) + dt*f(1);

                unew(i,j,k,0) = u(i,j,k,0) + dt*v(i,j,k,0);
                unew(i,j,k,1) = u(i,j,k,1) + dt*v(i,j,k,1);

            });

            // Nothing to do here.
        }
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/) override
    {
        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        a_tags.setVal(amrex::TagBox::CLEAR);
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &disp = u_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
                {
                    Set::Vector grad = Numeric::Gradient(disp, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > elastic.ref_criterion)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
        }
    }

private:
    Set::Field<Set::Scalar> u_mf;
    Set::Field<Set::Scalar> u_old_mf;
    Set::Field<Set::Scalar> trace_eps_mf;
    Set::Field<Set::Scalar> v_mf;
    Set::Field<Set::Scalar> v_old_mf;
    Set::Field<Set::Matrix> sig_mf;

    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> disp_mf;
    Set::Field<Set::Scalar> rhs_mf;
    Set::Field<Set::Scalar> res_mf;
    Set::Field<Set::Scalar> stress_mf;
    Set::Field<Set::Scalar> strain_mf;

    Set::Field<model_type> model_mf;

    IC::IC *ic;          ///< Pointer to abstract IC object
    BC::BC<Set::Scalar> *mybc;        ///< Stock generic boundary condition object
    
    Set::Matrix Fmatrix = Set::Matrix::Zero();
    Set::Matrix Finclusion = Set::Matrix::Zero();

    struct {
        Model::Solid::Affine::Isotropic model;
        Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic> *solver;//(elastic.op);
        BC::Operator::Elastic<model_type> bc;
        Set::Scalar ref_criterion = 1E-8;
    } elastic;

};
} // namespace Integrator
#endif
