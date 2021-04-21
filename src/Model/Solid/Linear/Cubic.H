#ifndef MODEL_SOLID_LINEAR_CUBIC_H_
#define MODEL_SOLID_LINEAR_CUBIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Cubic : public Solid<Set::Sym::MajorMinor>
{
public:

	Cubic() {};
	virtual ~Cubic() {};
    
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (DDW(gradu)*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return DDW(gradu)*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & /*gradu*/) const override
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> ddw;


        Eigen::Matrix3d R = q.normalized().toRotationMatrix();

    	for(int p = 0; p < AMREX_SPACEDIM; p++) 
    		for(int q = 0; q < AMREX_SPACEDIM; q++) 
    			for(int s = 0; s < AMREX_SPACEDIM; s++) 
    				for(int t = 0; t < AMREX_SPACEDIM; t++)
    				{
    					ddw(p,q,s,t) = 0.0;
    					for(int i = 0; i < 3; i++) 
    						for(int j = 0; j < 3; j++) 
    							for(int k = 0; k < 3; k++) 
    								for(int l = 0; l < 3; l++) 
                                    {
                    					if(i == j && j == k && k == l) ddw(p,q,s,t) += R(p,i)*R(s,k)*C11*R(q,j)*R(t,l);
                    					else if (i==k && j==l)         ddw(p,q,s,t) += R(p,i)*R(s,k)*C44*R(q,j)*R(t,l);
                    					else if (i==j && k==l)         ddw(p,q,s,t) += R(p,i)*R(s,k)*C12*R(q,j)*R(t,l);
                    					else continue;
                                    }
    				}
        return ddw;
    }

    virtual void Print(std::ostream &out) const //override 
    {
        out << "C11=" << C11 << ", C12=" << C12 << ", C44=" << C44 << " R=" << q.toRotationMatrix();;
    }

protected:

    Set::Scalar C11=NAN, C12=NAN, C44=NAN;
    Eigen::Quaterniond q;    

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;

    AMREX_FORCE_INLINE
    static Cubic Combine(const std::vector<Cubic> &models, const std::vector<Set::Scalar> &eta)
    {
        Cubic ret;
        Util::Abort(INFO,"not implemented yet");
        //Eigen::Quaterniond::Identity();
        //ret.q = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        //Set::Scalar etasum = 0.;
        //for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        //for (unsigned int n = 0 ; n < models.size(); n++)
        //{
        //    ret.ddw += models[n].ddw * (eta[n] / etasum);
        //}
        return ret;
    }

    static Cubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static Cubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        Cubic ret;
        ret.q = Eigen::Quaterniond::UnitRandom();
	    return ret;
    }

    static void Parse(Cubic & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = 1.68, C12 = 1.21, C44 = 0.75;
        pp.query("C11",C11);
        pp.query("C12",C12);
        pp.query("C44",C44);

        if (pp.contains("random"))
        {
            value = Cubic::Random(C11,C12,C44);
            return;
        }

        if (pp.contains("phi1") || pp.contains("Phi") || pp.contains("phi2"))
        {
            Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
            pp.query("phi1",phi1);
            pp.query("Phi",Phi);
            pp.query("phi2",phi2);

            Eigen::Matrix3d R;
            R = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
    		    Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
    	 	    Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
            value.q = Eigen::Quaterniond(R);
        }        
    }
//
//    AMREX_FORCE_INLINE
//	void operator += (const Cubic &rhs)
//    {
//        C11 += rhs.C11; C12 += rhs.C12; C44 += rhs.C44;
//        q.w += rhs.q.w; q.x += rhs.q.x; q.y += rhs.q.y; q.z += rhs.q.z;
//	}    
//    AMREX_FORCE_INLINE
//   	Cubic operator * (const Set::Scalar alpha) const
//	{
//        Cubic ret;
//        ret.C11 = C11*alpha; ret.C12 = C12*alpha; ret.C44 = C44*alpha;
//        ret.q.w = q.w*alpha; ret.q.x = q.x*alpha; ret.q.y = q.y*alpha; ret.q.z = q.z*alpha;        
//        return ret;
//	}    
//    friend Cubic operator * (const Set::Scalar alpha, const Cubic b);
//    friend Cubic operator + (const Cubic a, const Cubic b);
//    friend Cubic operator - (const Cubic a, const Cubic b);

    #define OP_CLASS Cubic
    #define OP_VARS X(C11) X(C12) X(C44) X(q.coeffs())
    #include "Model/Solid/InClassOperators.H"
    
    
};
#include "Model/Solid/ExtClassOperators.H"
//AMREX_FORCE_INLINE
//Cubic operator * (const Set::Scalar alpha, const Cubic b)
//{
//    Cubic ret;
//    ret.C11 = b.C11*alpha; ret.C12 = b.C12*alpha; ret.C44 = b.C44*alpha;
//    ret.q.w = b.q.w*alpha; ret.q.x = b.q.x*alpha; ret.q.y = b.q.y*alpha; ret.q.x = b.q.z*alpha;
//    return ret;
//}    
//AMREX_FORCE_INLINE
//Cubic operator + (const Cubic a, const Cubic b)
//{
//    Cubic ret;
//    ret.C11 = a.C11+b.C11; ret.C12 = a.C12+b.C12; ret.C44 = a.C44+b.C44;
//    ret.q.w = a.q.w+b.q.w; ret.q.x = a.q.x+b.q.x; ret.q.y = a.q.y+b.q.y; ret.q.z = a.q.z+b.q.z;
//    return ret;
//}    
//AMREX_FORCE_INLINE
//Cubic operator - (const Cubic a, const Cubic b)
//{
//    Cubic ret;
//    ret.C11 = a.C11-b.C11; ret.C12 = a.C12-b.C12; ret.C44 = a.C44-b.C44;
//    ret.q.w = a.q.w-b.q.w; ret.q.x = a.q.x-b.q.x; ret.q.y = a.q.y-b.q.y; ret.q.z = a.q.z-b.q.z;
//    return ret;
//}

}
}
}

#endif

