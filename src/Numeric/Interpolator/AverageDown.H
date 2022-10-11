#ifndef NUMERIC_INTERPOLATOR_AVERAGEDOWN_H
#define NUMERIC_INTERPOLATOR_AVERAGEDOWN_H

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_MultiFabUtil.H>

namespace Numeric
{
namespace Interpolator
{
    using namespace amrex;

    template<class FAB>
    void average_down (const FAB& S_fine, FAB& S_crse,
                       int scomp, int ncomp, const IntVect& ratio)
    {
        BL_PROFILE("amrex::average_down");
        AMREX_ASSERT(S_crse.nComp() == S_fine.nComp());
        AMREX_ASSERT((S_crse.is_cell_centered() && S_fine.is_cell_centered()) ||
                     (S_crse.is_nodal()         && S_fine.is_nodal()));

        bool is_cell_centered = S_crse.is_cell_centered();

        using value_type = typename FAB::value_type;

        //
        // Coarsen() the fine stuff on processors owning the fine data.
        //
        BoxArray crse_S_fine_BA = S_fine.boxArray(); crse_S_fine_BA.coarsen(ratio);

        if (crse_S_fine_BA == S_crse.boxArray() && S_fine.DistributionMap() == S_crse.DistributionMap())
        {
#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion() && S_crse.isFusingCandidate()) {
                auto const& crsema = S_crse.arrays();
                auto const& finema = S_fine.const_arrays();
                if (is_cell_centered) {
                    ParallelFor(S_crse, IntVect(0), ncomp,
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                    {
                        amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],scomp,scomp,ratio);
                    });
                }
                Gpu::streamSynchronize();
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(S_crse,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    //  NOTE: The tilebox is defined at the coarse level.
                    const Box& bx = mfi.tilebox();
                    Array4<value_type> const& crsearr = S_crse.array(mfi);
                    Array4<value_type const> const& finearr = S_fine.const_array(mfi);

                    if (is_cell_centered) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                        {
                            amrex_avgdown(i,j,k,n,crsearr,finearr,scomp,scomp,ratio);
                        });
                    }
                }
            }
        }
        else
        {
            //const DefaultFabFactory<FAB> defaultFabFactory;
            FAB crse_S_fine(crse_S_fine_BA, S_fine.DistributionMap(), ncomp, 0, MFInfo());

#ifdef AMREX_USE_GPU
            if (Gpu::inLaunchRegion() && crse_S_fine.isFusingCandidate()) {
                auto const& crsema = crse_S_fine.arrays();
                auto const& finema = S_fine.const_arrays();
                if (is_cell_centered) {
                    ParallelFor(crse_S_fine, IntVect(0), ncomp,
                    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
                    {
                        amrex_avgdown(i,j,k,n,crsema[box_no],finema[box_no],0,scomp,ratio);
                    });
                }
                Gpu::streamSynchronize();
            } else
#endif
            {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
                for (MFIter mfi(crse_S_fine,TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    //  NOTE: The tilebox is defined at the coarse level.
                    const Box& bx = mfi.tilebox();
                    Array4<value_type> const& crsearr = crse_S_fine.array(mfi);
                    Array4<value_type const> const& finearr = S_fine.const_array(mfi);

                    //  NOTE: We copy from component scomp of the fine fab into component 0 of the crse fab
                    //        because the crse fab is a temporary which was made starting at comp 0, it is
                    //        not part of the actual crse multifab which came in.

                    if (is_cell_centered) {
                        AMREX_HOST_DEVICE_PARALLEL_FOR_4D(bx, ncomp, i, j, k, n,
                        {
                            amrex_avgdown(i,j,k,n,crsearr,finearr,0,scomp,ratio);
                        });
                    } 
                }
            }

            S_crse.ParallelCopy(crse_S_fine,0,scomp,ncomp);
        }
   }



}
}

#endif
