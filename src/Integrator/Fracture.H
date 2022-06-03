#ifndef INTEGRATOR_FRACTURE_H
#define INTEGRATOR_FRACTURE_H

#include "Integrator/Mechanics.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Model/Interface/Crack/Constant.H"

#include "IC/IC.H"
#include "IC/Ellipsoid.H"
#include "IC/Notch.H"
#include "IC/Laminate.H"
#include "IC/PerturbedInterface.H"
#include "IC/BMP.H"

#include <eigen3/Eigen/Dense>

namespace Integrator
{
using brittle_model = Model::Solid::Linear::Isotropic;
using crack_type = Model::Interface::Crack::Constant;

class Fracture: 
    virtual public MechanicsBase<brittle_model>
{

public:
    Fracture(): MechanicsBase<brittle_model>(){}

    static void Parse(Fracture &value, IO::ParmParse &pp)
    {
        MechanicsBase<brittle_model>::Parse(value,pp);

        // Material field related parsing
        pp.query("material.refinement_threshold", value.material.m_eta_ref_threshold);
        std::string mat_ic_type;
        pp.query("material.ic.type", mat_ic_type);
        if(mat_ic_type == "laminate")
        {
            IC::Laminate *tmpic = new IC::Laminate(value.geom);
            pp.queryclass("material.ic.laminate", *tmpic);
            pp.query("material.ic.laminate.number_of_inclusions", value.material.num_mat);
            value.material.num_mat++;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "ellipse")
        {
            IC::Ellipse *tmpic = new IC::Ellipse(value.geom);
            pp.queryclass("material.ic.ellipse", *tmpic);
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "perturbed_interface")
        {
            IC::PerturbedInterface *tmpic = new IC::PerturbedInterface(value.geom);
            pp.queryclass("material.ic.perturbed_interface", *tmpic);
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else if(mat_ic_type == "bmp")
        {
            IC::BMP *tmpic = new IC::BMP(value.geom);
            pp.queryclass("material.ic.bmp", *tmpic);
            value.material.num_mat = 2;
            value.material.ic = tmpic;
            value.material.is_ic = true;
        }
        else
        {
            value.material.num_mat = 1;
            value.material.is_ic = false;
        }

        for (int i = 0; i<value.material.num_mat; i++)
        {
            std::string name = "material.model" + std::to_string(i+1);
            std::string name_type = "material.model" + std::to_string(i+1) + ".type";
            std::string mat_model_type;
            pp.query(name_type.c_str(),mat_model_type);
            if(mat_model_type == "isotropic")
            {
                brittle_model tmp_model;
                pp.queryclass(std::string((name+".isotropic").data()), tmp_model);
                value.material.models.push_back(tmp_model);
            }
            else
                Util::Abort(INFO, "This type of material model is not supported for fracture.");
        }
        Util::Assert(INFO,TEST(value.material.models.size() > 0));

        value.RegisterNodalFab(value.material.eta_mf, value.material.num_mat, 2, "eta", true);

        // Crack related parsing
        pp.query("crack.modulus_scaling_max", value.crack.scaleModulusMax);
        pp.query("crack.refinement_threshold",value.crack.refinement_threshold);
        pp.query("crack.df.tol_rel",value.crack.tol_rel);
        pp.query("crack.df.tol_abs",value.crack.tol_abs);
        pp.query("crack.df.mult_Gc", value.crack.mult_Gc);
        pp.query("crack.df.el_mult", value.crack.el_mult);

        pp.queryarr("crack.ic.type",value.crack.ic_type);
        if(value.crack.ic_type.size() == 0) value.crack.is_ic = false;
        else
        {
            value.crack.is_ic= true;
            for (unsigned int i = 0; i<value.crack.ic_type.size(); i++)
            {
                if(value.crack.ic_type[i] == "notch")
                {
                    IC::Notch *tmpic = new IC::Notch(value.geom);
                    pp.queryclass("crack.ic.notch",*tmpic);
                    value.crack.ic.push_back(tmpic);
                }
                else if(value.crack.ic_type[i] == "ellipsoid")
                {
                    IC::Ellipsoid *tmpic = new IC::Ellipsoid(value.geom);
                    pp.queryclass("crack.ic.ellipsoid",*tmpic);
                    value.crack.ic.push_back(tmpic);
                }
                else
                    Util::Abort(INFO, "This IC hasn't been implemented yet");
            }
        }
        value.RegisterNodalFab(value.crack.c_mf, 1, 2, "crack", true);
        value.RegisterNodalFab(value.crack.c_old_mf, 1, 2, "crack_old", true);
        value.RegisterNodalFab(value.crack.driving_force_mf, 5, 2, "driving_force", true);
        value.RegisterNodalFab(value.crack.energy_pristine_mf, 1, 2, "energy_pristine", true);
        value.RegisterNodalFab(value.crack.energy_pristine_old_mf, 1, 2, "energy_pristine_old", false);
        value.RegisterIntegratedVariable(&(value.crack.driving_force_norm),"driving_force_norm");

        for (int i = 0; i<value.material.num_mat; i++)
        {
            std::string name = "crack.model" + std::to_string(i+1);
            std::string name_type = "crack.model" + std::to_string(i+1) + ".type";
            std::string mat_model_type;
            pp.query(name_type.c_str(),mat_model_type);
            if(mat_model_type == "constant")
            {
                crack_type tmp_model;
                pp.queryclass(std::string((name+".constant").data()), tmp_model);
                value.crack.cracktype.push_back(tmp_model);
            }
            else
                Util::Abort(INFO, "This type of material model is not supported for fracture.");
        }
        Util::Assert(INFO,TEST(value.material.models.size() > 0));
    }

    void Initialize(int lev) override
    {
        MechanicsBase<brittle_model>::Initialize(lev);

        crack.c_mf[lev]->setVal(1.0);
        crack.c_old_mf[lev]->setVal(1.0);
        if(crack.is_ic)
        {
            for (int i=0; i < crack.ic.size(); i++)
            {
                crack.ic[i]->Add(lev,crack.c_mf);
                crack.ic[i]->Add(lev,crack.c_old_mf);
            }
        }
        if (material.is_ic) material.ic->Initialize(lev,material.eta_mf);
        else material.eta_mf[lev]->setVal(1.0);
    }

    virtual void UpdateModel(int a_step) override
    {
        if (m_type == MechanicsBase<brittle_model>::Type::Disable) return;

        // if (a_step == 0)
        {
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                material.eta_mf[lev]->FillBoundary();

                for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
                {
                    amrex::Box bx = mfi.grownnodaltilebox();
                    amrex::Array4<brittle_model> const &model = model_mf[lev]->array(mfi);
                    amrex::Array4<const Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                                model(i,j,k) = brittle_model::Zero();
                                                for (int n = 0; n < material.num_mat; n++) 
                                                    model(i,j,k) += eta(i,j,k,n) * material.models[n];
                                            });
                }
                Util::RealFillBoundary(*model_mf[lev],geom[lev]);
            }
        }

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            // model_mf[lev]->setVal(material.models[0]);
            crack.c_mf[lev]->FillBoundary();
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                amrex::Array4<brittle_model> const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                            Set::Scalar _tmp = 0.0;
                                            _tmp = crack.cracktype[0].g_phi(c(i,j,k,0),0) + crack.scaleModulusMax;
                                            if (std::isnan(_tmp)) Util::Abort(INFO);
                                            if(_tmp < 0.0) _tmp = crack.scaleModulusMax;
                                            if(_tmp > 1.0) _tmp = 1.0;
					
                                            model(i,j,k).DegradeModulus(1-_tmp);
                                        });
            }
            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }
    }

    void TimeStepBegin(Set::Scalar a_time, int a_step) override
    {
        // Deciding whether to do an elastic solve or not
        Util::Message(INFO,crack.driving_force_norm," ",crack.driving_force_reference," ",crack.tol_rel);
        if (crack.driving_force_norm / crack.driving_force_reference < crack.tol_rel)
            elastic_do_solve_now = true;
        if (crack.driving_force_norm < crack.tol_abs)
           elastic_do_solve_now = true;
        if (!elastic_do_solve_now) return;

        // Doing an elastic solve
        MechanicsBase<brittle_model>::TimeStepBegin(a_time,a_step);

        // Computing pristine energy
        for (int lev = 0; lev <= disp_mf.finest_level; lev++)
        {
            std::swap(crack.energy_pristine_mf[lev],crack.energy_pristine_old_mf[lev]);
            amrex::Box domain = geom[lev].Domain();
            domain.convert(amrex::IntVect::TheNodeVector());

            const amrex::Real* DX = geom[lev].CellSize();
            for (MFIter mfi(*disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::Array4<Set::Scalar>     const &eta   = material.eta_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>     const &strain  = strain_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>     const &energy  = crack.energy_pristine_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>     const &energy_old  = crack.energy_pristine_old_mf[lev]->array(mfi);
                
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                {
                    Set::Matrix eps = strain(i,j,k);
                    Eigen::SelfAdjointEigenSolver<Set::Matrix> eigensolver(eps);
                    Set::Vector eValues = eigensolver.eigenvalues();
                    Set::Matrix eVectors = eigensolver.eigenvectors();

                    Set::Matrix epsp = Set::Matrix::Zero();
                    Set::Matrix epsn = Set::Matrix::Zero();

                    for (int n = 0; n<AMREX_SPACEDIM; n++)
                    {
                        if(eValues(n) > 0.0) epsp += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                        else epsn += epsp += eValues(n)*(eVectors.col(n)*eVectors.col(n).transpose());
                    }

                    for (int n = 0; n < material.num_mat; n++)
                        energy(i,j,k) += eta(i,j,k,n)*material.models[n].W(epsp);
                    
                    if (energy(i,j,k) < energy_old(i,j,k)) energy(i,j,k) = energy_old(i,j,k);
                });
            }
            Util::RealFillBoundary(*crack.energy_pristine_mf[lev],geom[lev]);
        }

        integrate_variables_before_advance = false;
        integrate_variables_after_advance = true;
    }

    void Advance(int a_lev, amrex::Real a_time, amrex::Real a_dt)
    {
        // advance for crack field
        std::swap (crack.c_old_mf[a_lev], crack.c_mf[a_lev]);

        const Set::Scalar* DX = geom[a_lev].CellSize();
        amrex::Box domain(geom[a_lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        const amrex::Dim3 lo= amrex::lbound(domain), hi = amrex::ubound(domain);

        for (amrex::MFIter mfi(*crack.c_mf[a_lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.nodaltilebox();
            amrex::Array4<const Set::Scalar> const& c_old = crack.c_old_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& energy = crack.energy_pristine_mf[a_lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = material.eta_mf[a_lev]->array(mfi);

            amrex::Array4<Set::Scalar> const& c = crack.c_mf[a_lev]->array(mfi);
            amrex::Array4<Set::Scalar> const& df = crack.driving_force_mf[a_lev]->array(mfi);

            amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
                #if AMREX_SPACEDIM !=2
                Util::Abort(INFO, "This does not work for 1D or 3D yet.");
                #endif

                // First we are taking care of zero gradient BC
                if      (i == lo.x && j == lo.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i+1,j+1,k) - 0.5*c(i+2,j+2,k));
                else if (i == lo.x && j == hi.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i+1,j-1,k) - 0.5*c(i+2,j-2,k));
                else if (i == hi.x && j == lo.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i-1,j+1,k) - 0.5*c(i-2,j+2,k));
                else if (i == hi.x && j == hi.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i-1,j-1,k) - 0.5*c(i-2,j-2,k));
                else if (i == lo.x) c(i,j,k) = (2.0/3.0)*(2.0*c(i+1,j,k) - 0.5*c(i+2,j,k));
                else if (j == lo.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i,j+1,k) - 0.5*c(i,j+2,k));
                else if (i == hi.x) c(i,j,k) = (2.0/3.0)*(2.0*c(i-1,j,k) - 0.5*c(i-2,j,k));
                else if (j == hi.y) c(i,j,k) = (2.0/3.0)*(2.0*c(i,j-1,k) - 0.5*c(i,j-2,k));

                // Next, we are doing crack evolution
                else
                {
                    Set::Scalar rhs = 0.0;
                    Set::Scalar bilap = 
                        Numeric::Stencil<Set::Scalar, 4, 0, 0>::D(c_old,i,j,k,0,DX) +
                        Numeric::Stencil<Set::Scalar, 2, 2, 0>::D(c_old,i,j,k,0,DX) * 2.0 +
                        Numeric::Stencil<Set::Scalar, 0, 4, 0>::D(c_old,i,j,k,0,DX);
                    
                    df(i,j,k,0) = crack.cracktype[0].Dg_phi(c_old(i,j,k),0.0)*energy(i,j,k)*crack.el_mult;
                    rhs += df(i,j,k,0);

                    // std::array<Numeric::StencilType,AMREX_SPACEDIM> stencil = 
                    //     {AMREX_D_DECL(Numeric::StencilType::Central,Numeric::StencilType::Central,Numeric::StencilType::Central)};
                    Set::Scalar laplacian = Numeric::Laplacian(c_old,i,j,k,0,DX);

                    Set::Scalar Gc = 0.0;
                    Set::Scalar Zeta = 0.0;
                    Set::Scalar Threshold = 0.0;
                    Set::Scalar Mobility = 0.0;
                    Set::Scalar _temp_product = 1.0, _temp_product2 = 1.0;
                    for (int m = 0; m < material.num_mat; m++)
                    {
                        Gc += eta(i,j,k,m)*crack.cracktype[m].Gc(c_old(i,j,k,0));
                        Zeta += eta(i,j,k,m)*crack.cracktype[m].Zeta(c_old(i,j,k,0));
                        Threshold += eta(i,j,k,m)*crack.cracktype[m].DrivingForceThreshold(c_old(i,j,k,0));
                        Mobility += eta(i,j,k,m)*crack.cracktype[m].Mobility(c_old(i,j,k,0));
                        _temp_product *= eta(i,j,k,m);
                        _temp_product2 *= 0.5;
                    }
                    if (material.num_mat > 1) Gc *= (1.0 - _temp_product*(1.-crack.mult_Gc)/_temp_product2);

                    df(i,j,k,1) = Gc*crack.cracktype[0].Dw_phi(c_old(i,j,k,0),0.0)/(4.0*Zeta);
                    rhs += df(i,j,k,1);

                    df(i,j,k,2) = 2.0*Zeta*Gc*laplacian;
                    rhs -= df(i,j,k,2);

                    df(i,j,k,3) = crack.beta*bilap;
                    rhs += df(i,j,k,3);
                    
                    df(i,j,k,4) =                    std::max(0.,rhs - Threshold);
                    c(i,j,k,0) = c_old(i,j,k,0) - a_dt*df(i,j,k,4)*Mobility;

                    if (c (i,j,k,0) < 0.0) c(i,j,k,0) = 0.0;
                    if (c (i,j,k,0) > 1.0) c(i,j,k,0) = 1.0;
                }
            });
        }
        Util::RealFillBoundary(*crack.c_mf[a_lev],geom[a_lev]);

        // advance for mechanics field
        MechanicsBase<brittle_model>::Advance(a_lev,a_time,a_dt);
    }

    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar a_time, int a_ngrow) override
    {
        MechanicsBase<brittle_model>::TagCellsForRefinement(lev,a_tags,a_time,a_ngrow);

        Set::Vector DX(geom[lev].CellSize());
        Set::Scalar DXnorm = DX.lpNorm<2>();
        for (amrex::MFIter mfi(*material.eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &eta = material.eta_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                auto sten = Numeric::GetStencil(i,j,k,bx);
                {
                    Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX.data(),sten);
                    if (grad.lpNorm<2>() * DXnorm > material.m_eta_ref_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
        }
        for (amrex::MFIter mfi(*material.eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.nodaltilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar> const &c = crack.c_mf[lev]->array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                auto sten = Numeric::GetStencil(i,j,k,bx);
                {
                    Set::Vector grad = Numeric::Gradient(c, i, j, k, 0, DX.data(),sten);
                    if (grad.lpNorm<2>() * DXnorm > crack.refinement_threshold)
                        tags(i, j, k) = amrex::TagBox::SET;
                }
            });
        }
    }

    void Integrate(int amrlev, Set::Scalar time, int step,
                    const amrex::MFIter &mfi, const amrex::Box &a_box) override
    {
        Set::Vector DX(geom[amrlev].CellSize());
        const Set::Scalar DV = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        amrex::Array4<const Set::Scalar> const &df = crack.driving_force_mf[amrlev]->array(mfi);
        amrex::ParallelFor(a_box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            crack.driving_force_norm += df(i,j,k,4)*DV;
        });

        MechanicsBase<brittle_model>::Integrate(amrlev,time,step,mfi,a_box);
    }

    void TimeStepComplete(Set::Scalar /*time*/, int /* iter*/)
    {
        if (elastic_do_solve_now)
            crack.driving_force_reference = crack.driving_force_norm;
        elastic_do_solve_now = false;
    }

// Member variables
protected:
    struct{
        Set::Field<Set::Scalar> c_mf;
        Set::Field<Set::Scalar> c_old_mf;
        Set::Field<Set::Scalar> energy_pristine_mf;
        Set::Field<Set::Scalar> energy_pristine_old_mf;
        Set::Field<Set::Scalar> driving_force_mf;
        std::vector<crack_type> cracktype;

        std::vector<std::string> ic_type;
        std::vector<IC::IC*> ic;
        bool is_ic = true;

        Set::Scalar scaleModulusMax = 0.02;
        Set::Scalar refinement_threshold = 0.0001;
        Set::Scalar mult_Gc = 1.0;
        Set::Scalar el_mult = 1.0;

        Set::Scalar driving_force_reference = 1.0;
        Set::Scalar driving_force_norm = 0.0;
        Set::Scalar tol_rel = 1.e-3;
        Set::Scalar tol_abs = 1.e-3;

        Set::Scalar beta = 0.0;
    } crack;
    struct{
        Set::Field<Set::Scalar> eta_mf;
        Set::Scalar m_eta_ref_threshold = 0.01;
        std::vector<brittle_model> models;
        IC::IC *ic;
        bool is_ic = false;
        int num_mat = 1;
    } material;

    bool elastic_do_solve_now = true;

    using MechanicsBase<brittle_model>::m_type;
    using MechanicsBase<brittle_model>::finest_level;
    using MechanicsBase<brittle_model>::geom;
    using MechanicsBase<brittle_model>::model_mf;
};
}


#endif