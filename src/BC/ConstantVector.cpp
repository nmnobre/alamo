#include "ConstantVector.H"

namespace BC
{


void
ConstantVector::FillBoundary (amrex::BaseFab<Set::Vector> &a_in,
            const amrex::Box &a_box,
            int ngrow, int /*dcomp*/, int /*ncomp*/, amrex::Real /*time*/,
            Orientation face, const amrex::Mask * /*mask*/)
{
    const amrex::Real* DX = m_geom.CellSize();

    //Util::Assert(INFO,TEST(a_in.nComp() == (int)m_ncomp));

    amrex::Box box = a_box;
    box.grow(ngrow);
    const amrex::Dim3 lo= amrex::lbound(m_geom.Domain()), hi = amrex::ubound(m_geom.Domain());

    amrex::Array4<Set::Vector> const& in = a_in.array();

    for (unsigned int n = 0; n < m_ncomp; n++)
    amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k)
    {
        amrex::IntVect glevel;
        AMREX_D_TERM(glevel[0] = std::max(std::min(0,i-lo.x),i-hi.x); ,
                    glevel[1] = std::max(std::min(0,j-lo.y),j-hi.y); ,
                    glevel[2] = std::max(std::min(0,k-lo.z),k-hi.z); );
        
        if (glevel[0]<0 && (face == Orientation::xlo || face == Orientation::All)) // Left boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::XLO][n]))
                in(i,j,k)(n) = m_bc_val[Face::XLO](n);
            else if(BCUtil::IsNeumann(m_bc_type[Face::XLO][n]))
                in(i,j,k)(n) = in(i+1,j,k)(n) - (m_bc_val[Face::XLO].size() > 0 ? m_bc_val[Face::XLO](n)*DX[0] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::XLO][n]))
                in(i,j,k)(n) = in(1-glevel[0],j,k)(n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::XLO][n]))
                in(i,j,k)(n) = -in(1-glevel[0],j,k)(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::XLO][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions");
        }
        else if (glevel[0]>0 && (face == Orientation::xhi || face == Orientation::All)) // Right boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::XHI][n]))
                in(i,j,k)(n) = m_bc_val[Face::XHI](n);
            else if(BCUtil::IsNeumann(m_bc_type[Face::XHI][n]))
                in(i,j,k)(n) = in(i-1,j,k)(n) - (m_bc_val[Face::XHI].size() > 0 ? m_bc_val[Face::XHI](n)*DX[0] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::XHI][n]))
                in(i,j,k)(n) = in(hi.x-glevel[0],j,k)(n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::XHI][n]))
                in(i,j,k)(n) = -in(hi.x-glevel[0],j,k)(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::XHI][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions");
        }
        
        else if (glevel[1]<0 && (face == Orientation::ylo || face == Orientation::All)) // Bottom boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::YLO][n]))
                in(i,j,k)(n) = m_bc_val[Face::YLO](n);
            else if (BCUtil::IsNeumann(m_bc_type[Face::YLO][n]))
                in(i,j,k)(n) = in(i,j+1,k)(n) - (m_bc_val[Face::YLO].size() > 0 ? m_bc_val[Face::YLO](n)*DX[1] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::YLO][n]))
                in(i,j,k)(n) = in(i,j-glevel[1],k)(n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::YLO][n]))
                in(i,j,k)(n) = -in(i,j-glevel[1],k)(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::YLO][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions");
        }
        else if (glevel[1]>0 && (face == Orientation::yhi || face == Orientation::All)) // Top boundary
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::YHI][n]))
                in(i,j,k)(n) = m_bc_val[Face::YHI](n);
            else if (BCUtil::IsNeumann(m_bc_type[Face::YHI][n]))
                in(i,j,k)(n) = in(i,j-1,k)(n) - (m_bc_val[Face::YHI].size() > 0 ? m_bc_val[Face::YHI](n)*DX[1] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::YHI][n]))
                in(i,j,k)(n) = in(i,hi.y-glevel[1],k)(n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::YHI][n]))
                in(i,j,k)(n) = -in(i,hi.y-glevel[1],k)(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::YHI][n])) {}
            else
                Util::Abort(INFO, "Incorrect boundary conditions");
        }

#if AMREX_SPACEDIM>2
        else if (glevel[2]<0 && (face == Orientation::zlo || face == Orientation::All))
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::ZLO][n]))
                in(i,j,k)(n) = m_bc_val[Face::ZLO](n);
            else if (BCUtil::IsNeumann(m_bc_type[Face::ZLO][n]))
                in(i,j,k)(n) = in(i,j,k+1)(n) - (m_bc_val[Face::ZLO].size() > 0 ? m_bc_val[Face::ZLO](n)*DX[2] : 0);
            else if (BCUtil::IsReflectEven(m_bc_type[Face::ZLO][n]))
                in(i,j,k)(n) = in(i,j,1-glevel[2])(n);
            else if (BCUtil::IsReflectOdd(m_bc_type[Face::ZLO][n]))
                in(i,j,k)(n) = -in(i,j,1-glevel[2])(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::ZLO][n])) {}
            else Util::Abort(INFO, "Incorrect boundary conditions");
        }
        else if (glevel[2]>0 && (face == Orientation::zhi || face == Orientation::All))
        {
            if (BCUtil::IsDirichlet(m_bc_type[Face::ZHI][n]))
                in(i,j,k)(n) = m_bc_val[Face::ZHI](n);
            else if(BCUtil::IsNeumann(m_bc_type[Face::ZHI][n]))
                in(i,j,k)(n) = in(i,j,k-1)(n) - (m_bc_val[Face::ZHI].size() > 0 ? m_bc_val[Face::ZHI](n)*DX[2] : 0);
            else if(BCUtil::IsReflectEven(m_bc_type[Face::ZHI][n]))
                in(i,j,k)(n) = in(i,j,hi.z-glevel[2])(n);
            else if(BCUtil::IsReflectOdd(m_bc_type[Face::ZHI][n]))
                in(i,j,k)(n) = -in(i,j,hi.z-glevel[2])(n);
            else if(BCUtil::IsPeriodic(m_bc_type[Face::ZHI][n])) {}
            else Util::Abort(INFO, "Incorrect boundary conditions");
        }
#endif


    });
}

amrex::BCRec
ConstantVector::GetBCRec() 
{
    int bc_lo[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XLO][0],m_bc_type[Face::YLO][0],m_bc_type[Face::XLO][0])};
    int bc_hi[BL_SPACEDIM] = {AMREX_D_DECL(m_bc_type[Face::XHI][0],m_bc_type[Face::YHI][0],m_bc_type[Face::XHI][0])};

    return amrex::BCRec(bc_lo,bc_hi);
}

amrex::Array<int,AMREX_SPACEDIM>
ConstantVector::IsPeriodic()
{
    return {AMREX_D_DECL(BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))};
}
amrex::Periodicity ConstantVector::Periodicity () const
{
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(m_geom.Domain().length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                            m_geom.Domain().length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                            m_geom.Domain().length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));
}
amrex::Periodicity ConstantVector::Periodicity (const amrex::Box& b) {
    return amrex::Periodicity(amrex::IntVect(AMREX_D_DECL(b.length(0) * BCUtil::IsPeriodic(m_bc_type[Face::XLO][0]),
                                                        b.length(1) * BCUtil::IsPeriodic(m_bc_type[Face::YLO][0]),
                                                        b.length(2) * BCUtil::IsPeriodic(m_bc_type[Face::ZLO][0]))));

}


}
