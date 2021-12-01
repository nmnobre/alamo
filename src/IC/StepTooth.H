#ifndef IC_STEPTOOTH_H_
#define IC_STEPTOOTH_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

namespace IC
{
class StepTooth : public IC
{
public:
	StepTooth (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
	}
  
	void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
	{
		// Set::Scalar AMREX_D_DECL(L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0],
		// 						 L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1],
		// 						 L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2]);

		for (amrex::MFIter mfi(*a_field[lev],true); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::IndexType type = a_field[lev]->ixType();

			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				
				Set::Vector x;
                // NODE
                if (type == amrex::IndexType::TheNodeType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
                }
                else if (type == amrex::IndexType::TheCellType())
                {
				    AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
				    			 x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
				    			 x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
                }

#if AMREX_SPACEDIM == 2				

                // if(x(0) < 5.e-3 && x(0) > 3.e-3 && std::abs(x(1)) < 0.5*step_height)
                //     Util::Message(INFO, "x = ", x(0), ", y = ", x(1), ". Value = ", std::floor((std::abs(x(0) - center[0]) + 0.5*step_width)/step_width) % 2);

				if (x(1) - center[1] > offset + 0.5*step_height)
                {
                    field(i,j,k,0) = 1.;
                    field(i,j,k,1) = 0.;
                }
                else if (x(1) - center[1] < offset - 0.5*step_height)
                {
                    field(i,j,k,0) = 0.;
                    field(i,j,k,1) = 1.;
                }
                else if ((int)((std::abs(x(0) - center[0]) + 0.5*step_width)/step_width) % 2 ==  0)
                {
                    field(i,j,k,0) = 1.;
                    field(i,j,k,1) = 0.;
                }
                else
                {
                    field(i,j,k,0) = 0.;
                    field(i,j,k,1) = 1.;
                }
#endif
			});
		}

        const int W = 5;
        Eigen::Matrix<amrex::Real,W,W> kernel = Eigen::Matrix<amrex::Real,W,W>::Zero();
        Set::Scalar mean = W/2.0;
        
        for (int m = 0; m < W; ++m) 
            for (int n = 0; n < W; ++n) {
                kernel(m,n) = std::exp( -0.5 * (((m-mean)/sigma)*((m-mean)/sigma) + ((n-mean)/sigma)*((n-mean)/sigma))) / (2.0 * PI * sigma * sigma);
            }

        kernel /= kernel.sum();

        for (amrex::MFIter mfi(*a_field[lev],true); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
                for (int num_comp = 0; num_comp < 2; num_comp++)
                {
                    Eigen::Matrix<amrex::Real,W,W> field_matrix = Eigen::Matrix<amrex::Real,W,W>::Zero();
                    for (int m = -2; m < 3; m++)
                        for (int n=-2; n< 3; n++)
                            field_matrix(m+2,n+1) = field(i+m,j+n,k,num_comp);
                
                    field(i,j,k,num_comp) = (field_matrix.transpose()*kernel).trace();
                }
            });
        }

	};
  
private:
	enum Direction {X,Y,Z};
    Set::Vector center = Set::Vector::Zero();
	Set::Scalar offset = 0.0;
	Set::Scalar step_height = 0.0;
    Set::Scalar step_width = 0.0;
    Set::Scalar sigma = 0.0;

public:
	static void Parse(StepTooth & value, IO::ParmParse & pp)
    {
        amrex::Vector<Set::Scalar> nCenter;
        pp.queryarr("center",nCenter);
        if(nCenter.size() == AMREX_SPACEDIM)
        {
        	value.center = Set::Vector(AMREX_D_DECL(nCenter[0], nCenter[1], nCenter[2]));
        }
        else
        	Util::Abort(INFO, "Incompatible values in center");

        pp.query("step_height", value.step_height);
        pp.query("step_offset", value.offset);
        pp.query("step_width", value.step_width);
        pp.query("sigma",value.sigma);
    }

};
}
#endif
