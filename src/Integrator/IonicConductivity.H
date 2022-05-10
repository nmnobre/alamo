#ifndef INTEGRATOR_IONICCONDUCTIVITY_H
#define INTEGRATOR_IONICCONDUCTIVITY_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/Constant.H"

#include "IC/Cylinder.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"

#include "Numeric/Stencil.H"

namespace Integrator
{
class IonicConductivity : virtual public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    IonicConductivity() : Integrator()
    {}

    static void Parse(IonicConductivity &value, IO::ParmParse &pp)
    {
        //pp.query("heat.alpha", value.alpha);
        pp.query("heat.refinement_threshold", value.refinement_threshold);


        std::string type = "unset";
        pp.query("chi.ic.type",type);
        if (type == "sphere")          value.ic_chi = new IC::Sphere(value.geom,pp,"chi.ic.sphere");
        else if (type == "constant")   value.ic_chi = new IC::Constant(value.geom, pp, "chi.ic.constant");
        else if (type == "expression") value.ic_chi = new IC::Expression(value.geom,pp,"chi.ic.expression");
        else  Util::Abort(INFO,"Invalid ic.type ", type);

        type = "unset";
        pp.query("phi.ic.type",type);
        if (type == "sphere")          value.ic_phi = new IC::Sphere(value.geom,pp,"phi.ic.sphere");
        else if (type == "constant")   value.ic_phi = new IC::Constant(value.geom, pp, "phi.ic.constant");
        else if (type == "expression") value.ic_phi = new IC::Expression(value.geom,pp,"phi.ic.expression");
        else  Util::Abort(INFO,"Invalid ic.type ", type);
            
        type = "unset";
        pp.query("alpha.ic.type",type);
        if (type == "sphere")          value.ic_alpha = new IC::Sphere(value.geom,pp,"alpha.ic.sphere");
        else if (type == "constant")   value.ic_alpha = new IC::Constant(value.geom, pp, "alpha.ic.constant");
        else if (type == "expression") value.ic_alpha = new IC::Expression(value.geom,pp,"alpha.ic.expression");
        else  Util::Abort(INFO,"Invalid ic.type ", type);

            
        value.bc = new BC::Constant(1);
        pp.queryclass("chi.bc",*static_cast<BC::Constant *>(value.bc));
                
        value.RegisterNewFab(value.chi_mf,     value.bc, 1, 1, "chi",true);
        value.RegisterNewFab(value.chi_old_mf, value.bc, 1, 1, "chi_old",false);
        value.RegisterNewFab(value.phi_mf,     value.bc, 1, 1, "phi", true);
        value.RegisterNewFab(value.alpha_mf,   value.bc, 1, 1, "alpha", true);
    }

protected:

    /// \brief Use the #ic object to initialize #Temp
    void Initialize(int lev)
    {
        ic_chi->Initialize(lev,chi_old_mf);
        ic_chi->Initialize(lev,chi_mf);
        ic_phi->Initialize(lev,phi_mf);
        ic_alpha->Initialize(lev,alpha_mf);
    }

    /// \brief Integrate the heat equation
    void Advance(int lev, amrex::Real /*time*/, amrex::Real dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*chi_mf[lev], *chi_old_mf[lev]);
        
        // Get the cell size corresponding to this level
        const Set::Scalar *DX = geom[lev].CellSize();

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*chi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box &bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            amrex::Array4<const Set::Scalar> const &chi     = (*chi_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const &chi_new = (*chi_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const &phi     = (*phi_mf[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const &alpha   = (*alpha_mf[lev]).array(mfi);
        
            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                // Do the physics!
                // Note that Numeric::Laplacian is an inlined function so there is no overhead.
                // You can calculate the derivatives yourself if you want.
                //temp(i,j,k) = temp_old(i,j,k) + dt * alpha * Numeric::Laplacian(temp_old,i,j,k,0,DX);

                chi_new(i,j,k) = chi(i,j,k) + dt*(
                    Numeric::Gradient(alpha,i,j,k,0,DX).dot( Numeric::Gradient(chi,i,j,k,0,DX)) + 
                    alpha(i,j,k) * Numeric::Laplacian(chi,i,j,k,0,DX));

                //chi_new(i,j,k) = 

            });
        }
    }

    /// \brief Tag cells for mesh refinement based on temperature gradient
    void TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {
        // Get cell dimensions as done above.
        const Set::Scalar *DX = geom[lev].CellSize();
        // Calculate the diagonal.
        Set::Scalar dr  = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Iterate over the patches on this level
        for (amrex::MFIter mfi(*chi_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char>         const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Scalar>  const &phi = (*phi_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar>  const &alpha = (*alpha_mf[lev]).array(mfi);
            
            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                    {
                                        // Calculate the temperature gradient.
                                        Set::Vector gradphi = Numeric::Gradient(phi,i,j,k,0,DX);
                                        Set::Vector gradalpha = Numeric::Gradient(alpha,i,j,k,0,DX);

                                        if (gradphi.lpNorm<2>() * dr > refinement_threshold)
                                            tags(i,j,k) = amrex::TagBox::SET;
                                        if (gradalpha.lpNorm<2>() * dr > refinement_threshold)
                                            tags(i,j,k) = amrex::TagBox::SET;
                                    });
        }
    }

protected:
    Set::Field<Set::Scalar> chi_mf;         // charge density
    Set::Field<Set::Scalar> chi_old_mf;     // old charge density
    Set::Field<Set::Scalar> phi_mf;         // species indicator
    Set::Field<Set::Scalar> alpha_mf;       // polymer grading indicator


private:

    Set::Scalar refinement_threshold = 0.01; ///< Criterion for cell refinement

    IC::IC *ic_chi;
    IC::IC *ic_phi;
    IC::IC *ic_alpha;
    BC::BC<Set::Scalar> *bc;
};
} // namespace Integrator
#endif
