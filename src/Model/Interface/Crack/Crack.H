#ifndef MODEL_INTERFACE_CRACK_H
#define MODEL_INTERFACE_CRACK_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace Crack
{
class Crack
{
    public:
    enum GType {GSQUARE, GMULTIWELL, GPHI4C3, GCUBICM, GSQUAREP, GSQUAREPM};
    enum WType {WSQUARE, WMULTIWELL, WMULTIWELL2, WPHI4C3};
    
    Crack() {};
    
    AMREX_FORCE_INLINE
    virtual Set::Scalar w_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:        return (1.-c)*(1.-c);
            case WMULTIWELL:        return (1.-c)*(1.-c)*c*c;
            case WMULTIWELL2:    return (1.+c)*(1.+c)*(1.-c)*(1.-c);
            case WPHI4C3:        return 1. - 4.*c*c*c + 3.*c*c*c*c;
            
            default:            return (1.-c)*(1.-c);
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar g_phi(Set::Scalar c, Set::Scalar p=0.)
    {
        switch(g_type)
        {
            case GSQUARE:        return c*c;
            case GMULTIWELL:     return (2.-c)*(2.-c)*c*c;
            case GPHI4C3:         return 4.*c*c*c - 3.*c*c*c*c;
            case GSQUAREP:        return std::pow(c,2.*p);
            case GSQUAREPM:        return std::pow(c,2.*(std::pow(p,m_d_exponent)));
            case GCUBICM:        return m_d_exponent*(c*c*c - c*c) + 3.*c*c - 2.*c*c*c;

            default:             return c*c;
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar Dw_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:        return -2.*(1.-c);
            case WMULTIWELL:        return 4.*c*c*c - 6.*c*c + 2.*c;
            case WMULTIWELL2:    return 4.*c*c*c - 4.*c;
            case WPHI4C3:        return 12.*(c-1.)*c*c;
            
            default:            return -2.*(1.-c);
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar Dg_phi(Set::Scalar c, Set::Scalar p=0.)
    {
        switch(g_type)
        {
            case GSQUARE:        return 2.*c;
            case GMULTIWELL:     return 4.*c*c*c - 12.*c*c + 8.*c;
            case GPHI4C3:         return 12.*(1.-c)*c*c;
            case GSQUAREP:        return 2.*p*std::pow(c,2*p -1.);
            case GSQUAREPM:        return 2.*std::pow(p,m_d_exponent)*(std::pow(c, 2*std::pow(p,m_d_exponent)-1));
            case GCUBICM:        return m_d_exponent*(3.*c*c - 2.*c) + 6.*c - 6.*c*c;

            default:             return 2*c;
        }
    }

    virtual Set::Scalar Gc(Set::Scalar theta) = 0;
    virtual Set::Scalar DGc(Set::Scalar theta) = 0;
    virtual Set::Scalar DDGc(Set::Scalar theta) = 0;
    virtual Set::Scalar Zeta(Set::Scalar theta) = 0;
    virtual Set::Scalar Mobility(Set::Scalar theta) = 0;
    virtual Set::Scalar DrivingForceThreshold(Set::Scalar theta) = 0;

    void ExportToFile(std::string filename, amrex::Real dTheta)
    {
        std::ofstream outFile;
        outFile.open(filename);
    
        for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
            {
                outFile << theta << " " << Gc(theta) << std::endl;   
            }
        outFile.close();
    
    }
    void SetGType(const GType a_type)
    {
        g_type = a_type;
    }

    void SetWType(const WType a_type)
    {
        w_type = a_type;
    }
    
    void SetDuctileExponent(const Set::Scalar m)
    {
        m_d_exponent = m;
    }

	static void Parse(Crack & value, IO::ParmParse & pp)
    {
		std::string gtype = "";
        std::string wtype = "";
        pp.query("gtype", gtype);
        pp.query("wtype", wtype);

		std::map<std::string,Model::Interface::Crack::Crack::GType>  g_map;
	    g_map["square"] = Model::Interface::Crack::Crack::GType::GSQUARE;
	    g_map["multiwell"] = Model::Interface::Crack::Crack::GType::GMULTIWELL;
        g_map["4c3"] = Model::Interface::Crack::Crack::GType::GPHI4C3;
        g_map["squarep"] = Model::Interface::Crack::Crack::GType::GSQUAREP;
	    g_map["squarepexp"] = Model::Interface::Crack::Crack::GType::GSQUAREPM;
        g_map["cubicm"] = Model::Interface::Crack::Crack::GType::GCUBICM;

        std::map<std::string,Model::Interface::Crack::Crack::WType>  w_map;
        w_map["square"] = Model::Interface::Crack::Crack::WType::WSQUARE;
	    w_map["multiwell"] = Model::Interface::Crack::Crack::WType::WMULTIWELL;
        w_map["multiwell2"] = Model::Interface::Crack::Crack::WType::WMULTIWELL2;
        w_map["4c3"] = Model::Interface::Crack::Crack::WType::WPHI4C3;

        value.g_type = g_map[gtype];
        value.w_type = w_map[wtype];

		Set::Scalar ductile_exponent = 1.;
        if(gtype == "squarepexp"|| gtype == "cubicm") pp.query("exponent",ductile_exponent);
		value.m_d_exponent = ductile_exponent;
	}

protected:
    static constexpr amrex::Real pi = 3.14159265359;
    GType g_type = GType::GSQUARE;
    WType w_type = WType::WSQUARE;
    Set::Scalar m_d_exponent = 1.;
};
}
}
}

#endif
