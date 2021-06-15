#ifndef INTEGRATOR_TENSIONTESTSGCP_H
#define INTEGRATOR_TENSIONTESTSGCP_H
#include <iostream>
#include <fstream>
#include <iomanip>
#include <numeric>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"
#include "BC/Operator/Elastic/Constant.H"
#include "BC/Operator/Elastic/TensionTest.H"

#include "IC/Ellipse.H"
#include "Numeric/Stencil.H"

#include "Model/Solid/Solid.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "Model/Solid/Affine/CrystalPlastic.H"
#include "Model/Solid/Affine/StrainGradientCrystalPlastic.H"

#include "Operator/Operator.H"

#include "IO/ParmParse.H"

namespace Integrator
{
template<class MODEL>
class TensionTestSGCP : public Integrator
{
    //using model_type = Model::Solid::Affine::Isotropic;
public:
    /// \brief Read in parameters and register field variables
    TensionTestSGCP()
    {
        
        RegisterNodalFab(disp_mf, AMREX_SPACEDIM, 2, "disp",true);
        RegisterNodalFab(rhs_mf, AMREX_SPACEDIM, 2, "rhs",true);
        RegisterNodalFab(stress_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "stress",true);
        RegisterNodalFab(strain_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strain",true);
        

        {
            IO::ParmParse pp("elastic");
            int nmodels = 1;
            pp.query("nmodels",nmodels);
            for (int i = 0; i<nmodels; i++)
            {
                std::string name = "model" + std::to_string(i+1);
                MODEL tmp_model;
                pp.queryclass(name.data(), tmp_model);
                elastic.models.push_back(tmp_model);
            }
            Util::Assert(INFO,TEST(elastic.models.size() > 0));
            RegisterNodalFab(eta_mf, elastic.models.size(), 2, "eta",true);    

            std::string bc_type = "constant";
            pp.query("bc.type",bc_type);
            if (bc_type == "constant")
            {
                elastic.bc = new BC::Operator::Elastic::Constant();
                pp.queryclass("bc",static_cast<BC::Operator::Elastic::Constant*>(elastic.bc));
            }
            else if (bc_type == "tension_test")
            {
                elastic.bc = new BC::Operator::Elastic::TensionTest();
                pp.queryclass("bc.tension_test",static_cast<BC::Operator::Elastic::TensionTest*>(elastic.bc));
            }
            pp.query("solver.int",elastic.interval);
        }
        RegisterGeneralFab(model_mf, 1, 2);

        {
            IO::ParmParse pp("ic");
            std::string type;
            pp.query("type",type);
            if (type=="ellipse")
            {
                ic = new IC::Ellipse(geom);
                pp.queryclass("ellipse",static_cast<IC::Ellipse*>(ic));
            }
        }
        #if AMREX_SPACEDIM==3
        RegisterIntegratedVariable(&(trac_lo[0].data()[0]),"trac_lo_xx");
        #endif

        {
            IO::ParmParse pp("plastic");
            pp.query("on", plastic.on);
            if (plastic.on)
            {
                pp.query("hardening",plastic.hardening);
                // Hardcoded FCC in here..
                RegisterNodalFab(plastic.gamalpha_mf, 12, 2, "gamalpha", true);
                RegisterNodalFab(plastic.gamdotalpha_mf, 12, 2, "gamdotalpha", true);
                RegisterIntegratedVariable(&(plastic.gamdot_int),"gamdot_int");

                RegisterNodalFab(plastic.strainp_mf, AMREX_SPACEDIM * AMREX_SPACEDIM, 2, "strainp",true);
                RegisterNodalFab(plastic.strainp_eq_mf, 1, 2, "strainp_eq",true);
            }
        }
    }

protected:
    /// \brief Use the #ic object to initialize#Temp
    void Initialize(int lev) override
    {
        eta_mf[lev]->setVal(0.0);
        if (elastic.models.size() > 1 && ic)
        {
            ic->Initialize(lev, eta_mf);
        }
        else 
        {
            eta_mf[lev]->setVal(1.0);
        }

        disp_mf[lev]->setVal(0.);
        rhs_mf[lev]->setVal(0.);
        stress_mf[lev]->setVal(0.);
        strain_mf[lev]->setVal(0.);

        if (plastic.on)
        {
            plastic.strainp_mf[lev]->setVal(0.);
            plastic.strainp_eq_mf[lev]->setVal(0.0);
            plastic.gamalpha_mf[lev]->setVal(0.0);
            plastic.gamdotalpha_mf[lev]->setVal(0.0);
        }
        
        eta_mf[lev]->FillBoundary();
        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox();
            amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                model(i,j,k) = MODEL::Zero();
                if(plastic.on) model(i,j,k).SetHardening(plastic.hardening);
                for (unsigned int n = 0; n < elastic.models.size(); n++) 
                    model(i,j,k) += eta(i,j,k,n) * elastic.models[n];
            });
        }
        Util::RealFillBoundary(*model_mf[lev],geom[lev]);
    }

    void TimeStepBegin(Set::Scalar a_time, int a_iter) override
    {
        Util::Message(INFO,plastic.gamdot_int," ",plastic.gamdot_int_ref," ",plastic.gamdot_rel);

        // Has gamma dot calculation converged?
        if (plastic.on && plastic.gamdot_int/plastic.gamdot_int_ref < plastic.gamdot_rel)
            elastic.do_elastic_solve = true;

        if (!elastic.do_elastic_solve) return;
        // Has gamma 
        //if(a_iter%elastic.interval) return;
        // Set linear elastic model
        
        for (int lev = 0; lev <= finest_level; ++lev)
        {
//            rhs_mf[lev]->setVal(0.0);
//
//            eta_mf[lev]->FillBoundary();
//
//            Set::Vector DX(geom[lev].CellSize());
//
//            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
//            {
//                amrex::Box bx = mfi.grownnodaltilebox();
//
//                amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
//                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
//
//                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
//                    model(i,j,k) = MODEL::Zero();
//                    if(plastic.on) model(i,j,k).SetHardening(plastic.hardening);
//                    for (unsigned int n = 0; n < elastic.models.size(); n++) 
//                        model(i,j,k) += eta(i,j,k,n) * elastic.models[n];
//                });
//            }

            disp_mf[lev]->setVal(0.0);

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }

        elastic.bc->SetTime(a_time);
        elastic.bc->Init(rhs_mf,geom);

        amrex::LPInfo info;
        Operator::Elastic<MODEL::sym> elastic_op(Geom(0,finest_level), grids, DistributionMap(0,finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetBC(elastic.bc);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;

        IO::ParmParse pp("elastic");
        elastic.solver = new Solver::Nonlocal::Newton<MODEL>(elastic_op);
        pp.queryclass("solver",*elastic.solver);

        elastic.solver->solve(disp_mf,rhs_mf,model_mf,tol_rel,tol_abs);

        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            const amrex::Real* DX = geom[lev].CellSize();
            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                amrex::Array4<MODEL>        const &model = model_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>       const &stress  = stress_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>       const &strain  = strain_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &disp  = disp_mf[lev]->array(mfi);


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    std::array<Numeric::StencilType,AMREX_SPACEDIM>
					    sten = Numeric::GetStencil(i,j,k,bx);

                    if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp,i,j,k,DX,sten);
                        Set::Matrix P = model(i,j,k).DW(F);
                        Numeric::MatrixToField(stress,i,j,k,P);
                        Numeric::MatrixToField(strain,i,j,k,F);
                    }
                    else
                    {                    
                        Set::Matrix gradu = Numeric::Gradient(disp,i,j,k,DX,sten);
                        Set::Matrix sigma = model(i,j,k).DW(gradu);
                        Set::Matrix eps = 0.5*(gradu + gradu.transpose());
                        Numeric::MatrixToField(stress,i,j,k,sigma);
                        Numeric::MatrixToField(strain,i,j,k,eps);
                    }
                });
            }
            Util::RealFillBoundary(*stress_mf[lev],geom[lev]);
            Util::RealFillBoundary(*strain_mf[lev],geom[lev]);
        }
    }

    void Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt) override
    {
        if (plastic.on)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            Util::Warning(INFO,"Need to determine why tiling isn't working...");

            if (!elastic.do_elastic_solve)
            {
                for (amrex::MFIter mfi(*eta_mf[lev],false); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.nodaltilebox();
                    amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& stress = stress_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const& gamalpha = plastic.gamalpha_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& gamdotalpha = plastic.gamdotalpha_mf[lev]->array(mfi);

                    amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        std::array<Numeric::StencilType,AMREX_SPACEDIM> sten = Numeric::GetStencil(i,j,k,bx);
                        Set::Matrix sigma = Numeric::FieldToMatrix(stress,i,j,k);
                        std::array<Set::Vector, 12> a_gamdot_grad, a_gam_grad;
                        std::array<Set::Matrix, 12> a_gamdot_gradgrad, a_gam_gradgrad;
                        for (int m = 0; m<12; m++)
                        {
                            a_gamdot_grad[m] = Numeric::Gradient(gamdotalpha,i,j,k,m,DX,sten);
                            a_gam_grad[m] = Numeric::Gradient(gamalpha,i,j,k,m,DX,sten);
                            a_gamdot_gradgrad[m] = Numeric::Hessian(gamdotalpha,i,j,k,m,DX);
                            a_gam_gradgrad[m] = Numeric::Hessian(gamalpha,i,j,k,m,DX);
                        }
                        model(i,j,k).ComputePlasticStrainRate(sigma, a_gamdot_grad, a_gam_grad, a_gamdot_gradgrad, a_gam_gradgrad, dt);
                        for (int m=0; m<12; m++)
                            gamdotalpha(i,j,k,m) = model(i,j,k).gamdotalpha[m];
                    });
                }
            }
            else //if (plastic.do_gam)
            {
                for ( amrex::MFIter mfi(*eta_mf[lev],false); mfi.isValid(); ++mfi )
                {
                    const amrex::Box& bx = mfi.nodaltilebox();
                    amrex::Array4<MODEL> const &model = model_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& gamalpha = plastic.gamalpha_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& strainp = plastic.strainp_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& strainp_eq = plastic.strainp_eq_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &disp  = disp_mf[lev]->array(mfi);
                    amrex::Array4<Set::Scalar> const& stress = stress_mf[lev]->array(mfi);

                    amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        model(i,j,k).Advance(Set::Matrix::Zero(),Set::Matrix::Zero(),Set::Matrix3::Zero(),dt);
                        Set::Matrix epsp = model(i,j,k).GetPlasticState();
                        for (int m=0; m<12; m++)
                            gamalpha(i,j,k,m) = model(i,j,k).gamalpha[m];

                        std::array<Numeric::StencilType,AMREX_SPACEDIM> sten = Numeric::GetStencil(i,j,k,bx);
                        if (model(i,j,k).kinvar == Model::Solid::KinematicVariable::F)
                        {
                            Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp,i,j,k,DX,sten);
                            Set::Matrix P = model(i,j,k).DW(F);
                            Numeric::MatrixToField(stress,i,j,k,P);
                        }
                        else
                        {
                            Set::Matrix gradu = Numeric::Gradient(disp,i,j,k,DX,sten);
                            Set::Matrix sigma = model(i,j,k).DW(gradu);
                            Numeric::MatrixToField(stress,i,j,k,sigma);
                        }

                        Numeric::MatrixToField(strainp,i,j,k,epsp);
                        strainp_eq(i,j,k) = epsp.norm();
                    });
                }
            }
            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
            Util::RealFillBoundary(*plastic.strainp_mf[lev],geom[lev]);
            Util::RealFillBoundary(*plastic.gamdotalpha_mf[lev],geom[lev]);
            Util::RealFillBoundary(*plastic.gamalpha_mf[lev],geom[lev]);
            Util::RealFillBoundary(*plastic.strainp_eq_mf[lev],geom[lev]);
            Util::RealFillBoundary(*stress_mf[lev],geom[lev]);
        }
    }
    
    void Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
    				 const amrex::MFIter &mfi, const amrex::Box &box) override
    {
    	const amrex::Real *DX = geom[amrlev].CellSize();
        amrex::Box domain = geom[amrlev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());


    	Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        #if AMREX_SPACEDIM == 3
        Set::Vector da(AMREX_D_DECL(DX[1]*DX[2], DX[2]*DX[0], DX[0]*DX[1]));
        #endif

		const Dim3 lo= amrex::lbound(domain)/*, hi = amrex::ubound(domain)*/;
		const Dim3 /*boxlo= amrex::lbound(box),*/ boxhi = amrex::ubound(box);

    	amrex::Array4<const amrex::Real> const &stress = (*stress_mf[amrlev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& gamdotalpha = plastic.gamdotalpha_mf[amrlev]->array(mfi);
    	//amrex::Array4<const amrex::Real> const &disp   = (*disp_mf[amrlev]).array(mfi);
    	amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

            #if AMREX_SPACEDIM == 3

            for (int m = 0; m<12; m++)
                plastic.gamdot_int += gamdotalpha(i,j,k,m)*gamdotalpha(i,j,k,m)*dv;

            if (i == lo.x && (j < boxhi.y && k < boxhi.z))
            {
                trac_lo[0](0) += 0.25 * (stress(i,j,k,0) + stress(i,j+1,k,0) + stress(i,j,k+1,0) + stress(i,j+1,k+1,0)) * da(0);
            } 

            #endif

        });
            
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
            amrex::Array4<Set::Scalar> const &disp = eta_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& strainp_eq = plastic.strainp_eq_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
                {
                    Set::Vector grad = Numeric::Gradient(disp, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
                {
                    Set::Vector grad = Numeric::Gradient(strainp_eq, i, j, k, 0, DX.data());
                    if (grad.lpNorm<2>() * DXnorm > 0.01)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
        }
    }

    void TimeStepComplete(amrex::Real /*time*/,int /*iter*/)
    {
        if (elastic.do_elastic_solve && plastic.gamdot_int > 0)
            plastic.gamdot_int_ref = plastic.gamdot_int;
        
        elastic.do_elastic_solve = false;
    }

private:
    Set::Field<Set::Scalar> eta_mf;
    Set::Field<Set::Scalar> disp_mf;
    Set::Field<Set::Scalar> rhs_mf;
    Set::Field<Set::Scalar> res_mf;
    Set::Field<Set::Scalar> stress_mf;
    Set::Field<Set::Scalar> strain_mf;

    Set::Vector trac_lo[AMREX_SPACEDIM];
    Set::Vector trac_hi[AMREX_SPACEDIM];


    Set::Field<MODEL> model_mf;

    IC::IC *ic = nullptr;          ///< Pointer to abstract IC object
    BC::BC<Set::Scalar> *mybc;        ///< Stock generic boundary condition object
    
    struct {
        bool do_elastic_solve = true;
        std::vector<MODEL> models;
        Solver::Nonlocal::Newton<MODEL> *solver;//(elastic.op);
        BC::Operator::Elastic::Elastic *bc;
        int interval=1;
    } elastic;

    struct{
        bool on = false;
        bool hardening = false;
        Set::Scalar gamdot_int = 0.0, gamdot_int_ref = 1.0;
        Set::Scalar gamdot_rel = 1.E-4;
        Set::Field<Set::Scalar> gamalpha_mf;
        Set::Field<Set::Scalar> gamdotalpha_mf;
        Set::Field<Set::Scalar> strainp_mf;
        Set::Field<Set::Scalar> strainp_eq_mf;
    } plastic;

};
} // namespace Integrator
#endif
