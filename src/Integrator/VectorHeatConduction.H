#ifndef INTEGRATOR_VECTORHEATCONDUCTION_H
#define INTEGRATOR_VECTORHEATCONDUCTION_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "IO/ParmParse.H"
#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/ConstantVector.H"

#include "IC/Cylinder.H"
#include "IC/Sphere.H"
#include "IC/Constant.H"
#include "IC/Expression.H"

#include "Numeric/Stencil.H"

namespace Integrator
{
class VectorHeatConduction : virtual public Integrator
{
public:
    /// \brief Read in parameters and register field variables
    VectorHeatConduction() : Integrator()
    {}
    VectorHeatConduction(IO::ParmParse &pp) : VectorHeatConduction()
    {Parse(*this,pp);}

    static void Parse(VectorHeatConduction &value, IO::ParmParse &pp)
    {
        pp.query("heat.alpha", value.alpha);
        pp.query("heat.refinement_threshold", value.refinement_threshold);

        std::string type = "";
        pp.query("ic.type",type);
        if (type == "sphere")          value.ic = new IC::Sphere(value.geom,pp,"ic.sphere");
        else if (type == "constant")   value.ic = new IC::Constant(value.geom, pp, "ic.constant");
        else if (type == "expression") value.ic = new IC::Expression(value.geom,pp,"ic.expression");
            
        value.bc = new BC::ConstantVector(pp,"bc.temp");
                
        value.AddField<Set::Vector,Set::Hypercube::Cell>(value.temp_mf,     value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp",true,true);
        value.AddField<Set::Vector,Set::Hypercube::Cell>(value.temp_old_mf, value.bc, value.number_of_components, value.number_of_ghost_cells, "Temp_old",false,true);
    }

protected:

    /// \brief Use the #ic object to initialize #Temp
    void Initialize(int lev)
    {
        if (ic) 
        {
            ic->Initialize(lev,temp_mf);
            ic->Initialize(lev,temp_old_mf);
        }
        else 
        {
            temp_mf[lev]->setVal(Set::Vector::Zero());
            temp_old_mf[lev]->setVal(Set::Vector::Zero());
        }
    }

    /// \brief Integrate the heat equation
    void Advance(int lev, amrex::Real /*time*/, amrex::Real dt)
    {
        // Swap the old temp fab and the new temp fab so we use
        // the new one.
        std::swap(*temp_mf[lev], *temp_old_mf[lev]);
        
        // Get the cell size corresponding to this level
        const amrex::Real *DX = geom[lev].CellSize();

        // Iterate over all of the patches on this level
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box (index dimensions) for this patch
            const amrex::Box &bx = mfi.tilebox();

            // Get an array-accessible handle to the data on this patch.
            amrex::Array4<const Set::Vector> const &temp_old = (*temp_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Vector>       const &temp     = (*temp_mf[lev]).array(mfi);
        
            // Iterate over the grid on this patch
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                    {
                                        // Do the physics!
                                        // Note that Numeric::Laplacian is an inlined function so there is no overhead.
                                        // You can calculate the derivatives yourself if you want.
                                        temp(i,j,k) = temp_old(i,j,k) + dt * alpha * Numeric::Laplacian(temp_old,i,j,k,DX);
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
        for (amrex::MFIter mfi(*temp_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            // Get the box and handles as done above.
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char>         const &tags = a_tags.array(mfi);
            amrex::Array4<Set::Vector>  const &temp = (*temp_mf[lev]).array(mfi);
            
            // Iterate over the grid as done above.
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                                    {
                                        // Calculate the temperature gradient.
                                        Set::Matrix grad = Numeric::Gradient(temp,i,j,k,DX);

                                        Set::Scalar norm = sqrt((grad.transpose()*grad).trace());
    
                                        // Is the gradient * cell_size too big? If so, then
                                        // mark this cell as needing refinement.
                                        if (norm * dr > refinement_threshold)
                                            tags(i,j,k) = amrex::TagBox::SET;
                                    });
        }
    }

protected:
    Set::Field<Set::Vector> temp_mf;         ///< Temperature field variable (current timestep)
    Set::Field<Set::Vector> temp_old_mf;     ///< Temperature field variable (previous timestep)

private:
    int number_of_components = 1;            ///< Number of components
    int number_of_ghost_cells = 2;           ///< Number of ghost cells


    Set::Scalar alpha = 1.0;                 ///< Thermal diffusivity
    Set::Scalar refinement_threshold = 0.01; ///< Criterion for cell refinement

    IC::IC *ic = nullptr;                    ///< Object used to initialize temperature field
    BC::BC<Set::Vector> *bc;                 ///< Object used to update temp field boundary ghost cells
};
} // namespace Integrator
#endif
