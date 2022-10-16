
#include <eigen3/Eigen/Eigenvalues>

#include <omp.h>
#include <cmath>
#include <chrono>
#include <random>

#include <AMReX_SPACE.H>
#include <AMReX_ParallelReduce.H>

#include "PhaseFieldMicrostructure.H"
#include "Integrator/Base/Mechanics.H"
#include "BC/Constant.H"
#include "BC/Step.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Random.H"
#include "IC/Trig.H"
#include "IC/Sphere.H"
#include "IC/Expression.H"
#include "Model/Interface/GB/SH.H"
#include "Numeric/Stencil.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "IC/Trig.H"
//#include "Numeric/Function/Mix.H"

#include "Util/MPI.H"

namespace Integrator
{
    void PhaseFieldMicrostructure::Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {
        BL_PROFILE("PhaseFieldMicrostructure::Advance");
        Base::Mechanics<model_type>::Advance(lev, time, dt);
        /// TODO Make this optional
        // if (lev != max_level) return;
        std::swap(eta_old_mf[lev], eta_new_mf[lev]);
        const Set::Scalar *DX = geom[lev].CellSize();

        Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

        for (amrex::MFIter mfi(*eta_new_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            // if (m_type == MechanicsBase<model_type>::Type::Static)
            // bx.grow(number_of_ghost_cells-1);
            // bx = bx & domain;
            amrex::Array4<const amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);
            amrex::Array4<amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
            // amrex::Array4<amrex::Real> const &fluct = (*fluct_mf[lev]).array(mfi);
            amrex::Array4<const amrex::Real> const &elasticdf = (*elasticdf_mf[lev]).array(mfi);
            //amrex::Array4<Set::Scalar> const &totaldf = (*totaldf_mf[lev]).array(mfi);
            //amrex::Array4<const Set::Matrix> const &sigma = (*stress_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (anisotropy.on && time >= anisotropy.tstart) 
                    boundary->SetMix(eta,i,j,k);
                for (int m = 0; m < number_of_grains; m++)
                {
                    Set::Scalar driving_force = 0.0;
                    Set::Scalar kappa = NAN, mu = NAN;

                    //
                    // BOUNDARY TERM and SECOND ORDER REGULARIZATION
                    //

                    Set::Vector Deta = Numeric::Gradient(eta, i, j, k, m, DX);
                    Set::Scalar normgrad = Deta.lpNorm<2>();
                    if (normgrad < 1E-4)  continue; // This ought to speed things up.
                    if (std::isnan(normgrad)) Util::Abort(INFO,"normgrad nan");

                    Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, m, DX);
                    Set::Scalar laplacian = DDeta.trace();

                    if (!anisotropy.on || time < anisotropy.tstart)
                    {
                        kappa = pf.l_gb * 0.75 * pf.sigma0;
                        mu = 0.75 * (1.0 / 0.23) * pf.sigma0 / pf.l_gb;
                        driving_force += -kappa * laplacian;
                    }
                    else
                    {
                        Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDeta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);
                        Util::Message(INFO,Deta);
                        Util::Message(INFO,DDeta);
                        //Util::Message(INFO,DDDDeta);
                        auto anisotropic_df = boundary->DrivingForce(Deta, DDeta, DDDDeta);
                        driving_force += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                        if (std::isnan(std::get<0>(anisotropic_df))) Util::Abort(INFO); // this is returning nan
                        driving_force += anisotropy.beta * std::get<1>(anisotropic_df);
                        if (std::isnan(std::get<1>(anisotropic_df))) Util::Abort(INFO);
                        mu = 0.75 * (1.0/0.23) * boundary->W(Deta) / pf.l_gb;
                    }

                    //
                    // CHEMICAL POTENTIAL
                    //

                    Set::Scalar sum_of_squares = 0.;
                    for (int n = 0; n < number_of_grains; n++)
                    {
                        if (m == n)
                            continue;
                        sum_of_squares += eta(i, j, k, n) * eta(i, j, k, n);
                    }
                    driving_force += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);

                    //
                    // LAGRANGE MULTIPLIER TERM FOR CONSERVING VOLUME
                    //
                    if (lagrange.on && m == 0 && time > lagrange.tstart)
                    {
                        driving_force += lagrange.lambda * (volume - lagrange.vol0);
                    }

                    //
                    // SYNTHETIC DRIVING FORCE
                    //
                    if (sdf.on && time > sdf.tstart)
                    {
                        driving_force += sdf.val[m](time);
                    }

                    //
                    // ELASTIC
                    //
                    if (m_type != Base::Mechanics<model_type>::Disable && time > mechanics.tstart)
                    {
                        //Util::Abort(INFO,"Need to implement");
                        driving_force += pf.elastic_mult * Numeric::Interpolate::NodeToCellAverage(elasticdf,i,j,k,m);
                        //Numeric::Interpolate::NodeToCellAverage(elasticdf,i,j,k,m);
                    }

                    //
                    // EVOLVE ETA
                    //
                    if (time < mechanics.tstart)
                        etanew(i, j, k, m) = eta(i, j, k, m) - pf.M * dt * driving_force;
                    else
                    {

                        if (driving_force > pf.elastic_threshold)
                        {
                            //totaldf(i,j,k,m) = (driving_force-pf.elastic_threshold);
                            etanew(i, j, k, m) = eta(i, j, k, m) - pf.M * dt * (driving_force-pf.elastic_threshold);
                        }
                        else if (driving_force < -pf.elastic_threshold)
                        {
                            //totaldf(i,j,k,m) = (driving_force + pf.elastic_threshold);
                            etanew(i, j, k, m) = eta(i, j, k, m) - pf.M * dt * (driving_force + pf.elastic_threshold);
                        }
                        else
                        {
                            //totaldf(i,j,k,m) = 0.0;
                        }
                    }
                    //totaldf(i,j,k,m) = driving_force;


                    //
                    // FLUCTUATION TERM
                    //

                    if (fluctuation.on && time > fluctuation.tstart)
                    {
                        etanew(i,j,k,m) += (fluctuation.amp * fluctuation.norm_dist(fluctuation.rand_num_gen) / DX[0]) * dt;
                    }

                    if (std::isnan(driving_force)) 
                        Util::Abort(INFO, i, " ", j, " ", k, " ", m);
                } 
            });
        }

        if (time < disconnection.tstart) return;

        UpdateEigenstrain(lev);
    }

void PhaseFieldMicrostructure::UpdateEigenstrain(int lev)
{
    eta_new_mf[lev]->FillBoundary();
    for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box bx = mfi.grownnodaltilebox();// & domain;
        amrex::Array4<const Set::Scalar> const &etaold = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const &etanew = (*eta_new_mf[lev]).array(mfi);
        //amrex::Array4<Set::Scalar> const &disc = (*disc_mf[lev]).array(mfi);
        amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            for (int m = 0; m < number_of_grains; m++)
                for (int n = m+1; n < number_of_grains; n++)
                {
                    Set::Scalar emnew = Numeric::Interpolate::CellToNodeAverage(etanew, i, j, k, m, sten);
                    Set::Scalar emold = Numeric::Interpolate::CellToNodeAverage(etaold, i, j, k, m, sten);
                    Set::Scalar ennew = Numeric::Interpolate::CellToNodeAverage(etanew, i, j, k, n, sten);
                    Set::Scalar enold = Numeric::Interpolate::CellToNodeAverage(etaold, i, j, k, n, sten);
                    Set::Scalar gmnew = (emnew * emnew) / (emnew * emnew + ennew * ennew);
                    Set::Scalar gmold = (emold * emold) / (emold * emold + enold * enold);
                    model(i, j, k, m).F0 -= (gmnew - gmold) * shearcouple.Fgb[m*number_of_grains + n]; 
                }
        });
    }
    Util::RealFillBoundary(*model_mf[lev],geom[lev]);

}

void PhaseFieldMicrostructure::Initialize(int lev)
{
    BL_PROFILE("PhaseFieldMicrostructure::Initialize");
    Base::Mechanics<model_type>::Initialize(lev);
    ic->Initialize(lev, eta_new_mf);
    ic->Initialize(lev, eta_old_mf);
}

void PhaseFieldMicrostructure::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("PhaseFieldMicrostructure::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev,a_tags,time,ngrow);
    const amrex::Real *DX = geom[lev].CellSize();
    const Set::Vector dx(DX);
    const Set::Scalar dxnorm = dx.lpNorm<2>();

    for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
        amrex::Array4<char> const &tags = a_tags.array(mfi);

        for (int n = 0; n < number_of_grains; n++)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        Set::Vector grad = Numeric::Gradient(etanew, i, j, k, n, DX);

                                        if (dxnorm * grad.lpNorm<2>() > ref_threshold)
                                            tags(i, j, k) = amrex::TagBox::SET;
                                    });
    }
}

void PhaseFieldMicrostructure::TimeStepComplete(Set::Scalar /*time*/, int /*iter*/)
{
    // TODO: remove this function, it is no longer needed.
}

void PhaseFieldMicrostructure::UpdateModel(int a_step)
{
    BL_PROFILE("PhaseFieldMicrostructure::UpdateModel");
    if (a_step % m_interval) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

        eta_new_mf[lev]->FillBoundary();
        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            
            amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &eta = eta_new_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                Set::Scalar sumsq = 0;
                for (int n = 0; n < number_of_grains; n++) sumsq += eta(i,j,k,n)*eta(i,j,k,n);

                model_type mix = model_type::Zero();
                for (int n = 0; n < number_of_grains; n++)
                {
                    mix += eta(i,j,k,n)*eta(i,j,k,n)*mechanics.model[n];
                }

                if (a_step == 0) // If we are just starting, do a complete model initialization
                    model(i,j,k) = mix;
                else // Otherwise, we will set the modulus only and leave the eigenstrain alone
                    model(i,j,k).ddw = mix.ddw;
            });
        }

        Util::RealFillBoundary(*model_mf[lev],geom[lev]);
    }

}

void PhaseFieldMicrostructure::TimeStepBegin(Set::Scalar time, int iter)
{
    BL_PROFILE("PhaseFieldMicrostructure::TimeStepBegin");
    for (int lev = 0; lev < totaldf_mf.size(); lev++) totaldf_mf[lev]->setVal(0.0);

    //
    // Manual Disconnection Nucleation
    //

    if (disconnection.on && time > disconnection.tstart && !(iter % disconnection.interval))
    {
        disconnection.sitex.clear();
        disconnection.sitey.clear();
        disconnection.phases.clear();

        int lev = max_level;
        const Set::Scalar *DX = geom[lev].CellSize();
        Set::Scalar exponent = DX[0] * DX[0] * (timestep * disconnection.interval) / disconnection.tau_vol;

        if (!disconnection.fixed.on)
        {
            // Determine the nucleation sites in my portion of the mesh
            for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                amrex::Array4<amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar E0 = 2.0*disconnection.nucleation_energy;
                    E0 /= disconnection.epsilon + 256.0*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1);
                    Set::Scalar p = std::exp(-E0/(disconnection.K_b*disconnection.temp));
                    Set::Scalar P = 1.0 - std::pow(1.0 - p,exponent);
                    if (eta(i,j,k,0) < 0 || eta(i,j,k,0) > 1.0 || eta(i,j,k,1) < 0 || eta(i,j,k,1) > 1.0) P = 0.0;
                    Set::Scalar q = 0.0;
                    q = disconnection.unif_dist(disconnection.rand_num_gen);
                    if (q < P)
                    {
                        disconnection.sitex.push_back(geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0]);
                        disconnection.sitey.push_back(geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1]);
                        int phase = disconnection.int_dist(disconnection.rand_num_gen);
                        disconnection.phases.push_back(phase);
                    } });
            }
            // Sync up all the nucleation sites among processors
            Util::MPI::Allgather(disconnection.sitex);
            Util::MPI::Allgather(disconnection.sitey);
            Util::MPI::Allgather(disconnection.phases);
            Util::Message(INFO, "Nucleating ", disconnection.phases.size(), " disconnections");
        }
        else
        {
            for (int n = 0; n < disconnection.fixed.sitex.size(); n++)
            {
                if (time > disconnection.fixed.time[n] && !disconnection.fixed.done[n])
                {
                    disconnection.sitex.push_back(disconnection.fixed.sitex[n]);
                    disconnection.sitey.push_back(disconnection.fixed.sitey[n]);
                    disconnection.phases.push_back(disconnection.fixed.phases[n]);
                    disconnection.fixed.done[n] = true;
                }
            }
        }


        if (disconnection.sitex.size() > 0)
        {
            // Now that we all know the nucleation locations, perform the nucleation
            for (int lev = 0; lev <= max_level; lev++)
            {
                amrex::Box domain = geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());

                const amrex::Real *DX = geom[lev].CellSize();
                // for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                for (amrex::MFIter mfi(*model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box bx = mfi.grownnodaltilebox() & domain;
                    //amrex::Array4<const Set::Scalar> const &etaold = (*eta_old_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const &etanew = (*eta_new_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const &disc = (*disc_mf[lev]).array(mfi);
                    //amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                       {
                                           //auto sten = Numeric::GetStencil(i, j, k, bx);
                                           Set::Vector x;
                                           AMREX_D_TERM(
                                               x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0];,
                                                                                                         x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1];
                                               ,
                                               x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * DX[2];);
                                           disc(i, j, k, 0) = 0.0;
                                           for (unsigned int m = 0; m < disconnection.phases.size(); m++)
                                           {
                                               // Util::ParallelMessage(INFO,disconnection.phases[m]);
                                               amrex::Real r_squared = 0;
                                               Set::Vector nucleation_site(disconnection.sitex[m], disconnection.sitey[m]);
                                               for (int n = 0; n < AMREX_SPACEDIM; n++)
                                               {
                                                   amrex::Real dist = nucleation_site(n) - x(n);
                                                   r_squared += dist * dist;
                                               }
                                               // amrex::Real bump = exp(1 - 1 / (1 - 2/disconnection.box_size * r_squared));
                                               amrex::Real bump = exp(-r_squared / disconnection.box_size);
                                               disc(i, j, k, 0) += bump;
                                               etanew(i, j, k, disconnection.phases[m]) = bump * (1 - etanew(i, j, k, disconnection.phases[m])) + etanew(i, j, k, disconnection.phases[m]);
                                               etanew(i, j, k, 1 - disconnection.phases[m]) = (1. - bump) * etanew(i, j, k, 1 - disconnection.phases[m]);
                                           }
                                       });
                }

                UpdateEigenstrain(lev);
                Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            }
        }
    }

    // 
    // Elastic solve
    //
    if (anisotropy.on && time >= anisotropy.tstart)
    {
        SetTimestep(anisotropy.timestep);
        if (anisotropy.elastic_int > 0) 
            if (iter % anisotropy.elastic_int) return;
    }
    if (time < mechanics.tstart) return;
    Base::Mechanics<model_type>::TimeStepBegin(time, iter);

    // Calculate strain energy
    for (int lev = 0; lev < rhs_mf.size(); ++lev)
    {
        amrex::Box domain(geom[lev].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());
        domain.grow(-1);
        elasticdf_mf[lev]->setVal(0.0);
        const Set::Scalar *DX = geom[lev].CellSize();
        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            amrex::Array4<const model_type> const &model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Vector> const &disp = disp_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &elasticdf = elasticdf_mf[lev]->array(mfi);
            //amrex::Array4<Set::Matrix> const &strain = strain_mf[lev]->array(mfi);
            amrex::Array4<Set::Scalar> const &eta = (*eta_new_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                auto sten = Numeric::GetStencil(i,j,k,bx);
                Set::Matrix gradu = Numeric::Gradient(disp, i, j, k, DX,sten);
                Set::Matrix P = model(i,j,k).DW(gradu);

                for (int m = 0; m < number_of_grains; m++)
                    for (int n = m+1; n < number_of_grains; n++)
                    {
                        Set::Scalar etam = Numeric::Interpolate::CellToNodeAverage(eta,i,j,k,m);
                        Set::Scalar etan = Numeric::Interpolate::CellToNodeAverage(eta,i,j,k,n);

                        Set::Scalar sumsq = (etam*etam + etan*etan);
                        sumsq = sumsq * sumsq;

                        Set::Scalar dgm = 2.0*etan*etan*etam / sumsq;
                        Set::Scalar dgn = 2.0*etan*etam*etam / sumsq;

                        Set::Matrix Fgbm = shearcouple.Fgb[m*number_of_grains + n];
                        Set::Matrix Fgbn = shearcouple.Fgb[n*number_of_grains + m];

                        elasticdf(i,j,k,m) = (P.transpose() * Fgbm).trace() * dgm;
                        elasticdf(i,j,k,n) = (P.transpose() * Fgbn).trace() * dgn;
                    }
                

                
            });
        }
        Util::RealFillBoundary(*elasticdf_mf[lev], geom[lev]);
    }
}

void PhaseFieldMicrostructure::Integrate(int amrlev, Set::Scalar time, int step,
                                        const amrex::MFIter &mfi, const amrex::Box &box)
{
    BL_PROFILE("PhaseFieldMicrostructure::Integrate");
    Base::Mechanics<model_type>::Integrate(amrlev,time,step,mfi,box);

    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);
    const amrex::Real *DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

    amrex::Array4<amrex::Real> const &eta = (*eta_new_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        #if AMREX_SPACEDIM == 2
        auto sten = Numeric::GetStencil(i,j,k,box);
        #endif

        volume += eta(i, j, k, 0) * dv;

        Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
        Set::Scalar normgrad = grad.lpNorm<2>();

        if (normgrad > 1E-8)
        {
            Set::Vector normal = grad / normgrad;

            Set::Scalar da = normgrad * dv;
            area += da;

            if (!anisotropy.on || time < anisotropy.tstart)
            {
                gbenergy += pf.sigma0 * da;
                Set::Scalar k = 0.75 * pf.sigma0 * pf.l_gb;
                realgbenergy += 0.5 * k * normgrad * normgrad * dv;
                regenergy = 0.0;
            }
            else
            {
            #if AMREX_SPACEDIM == 2
                Set::Scalar theta = atan2(grad(1), grad(0));
                Set::Scalar sigma = boundary->W(theta);
                gbenergy += sigma * da;

                Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, 0, DX, sten);
                Set::Vector tangent(normal[1], -normal[0]);
                Set::Scalar k2 = (DDeta * tangent).dot(tangent);
                regenergy += 0.5 * anisotropy.beta * k2 * k2;
            #elif AMREX_SPACEDIM == 3
                gbenergy += gbmodel.W(normal) * da;
            #endif
            }
        }
    });
}

} // namespace Integrator
