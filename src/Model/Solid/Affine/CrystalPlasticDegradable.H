#ifndef MODEL_SOLID_PLASTIC_CRYSTALPLASTICDEGRADABLE_H_
#define MODEL_SOLID_PLASTIC_CRYSTALPLASTICDEGRADABLE_H_

#include "AMReX.H"
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>
#include "Affine.H"
#include "Set/Set.H"
#include "Model/Solid/Affine/CubicDegradable.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class CrystalPlasticDegradable : public CubicDegradable
{
public:
    CrystalPlasticDegradable(){};
    CrystalPlasticDegradable(Solid<Set::Sym::MajorMinor> base) : CubicDegradable(base) {};
    CrystalPlasticDegradable(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Set::Scalar a_phi1, Set::Scalar a_Phi, Set::Scalar a_phi2, LatticeType lattice=LatticeType::FCC)
    {
        Eigen::Matrix3d m;
    	m =     Eigen::AngleAxisd(a_phi2, Eigen::Vector3d::UnitX()) *
    		Eigen::AngleAxisd(a_Phi,  Eigen::Vector3d::UnitZ()) *
    	 	Eigen::AngleAxisd(a_phi1, Eigen::Vector3d::UnitX());
        Define(a_C11, a_C12, a_C44, m, lattice);
    }
    CrystalPlasticDegradable(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m, LatticeType lattice=LatticeType::FCC)
    {
        Define(a_C11, a_C12, a_C44, m, lattice);
    }

    void Define(Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m, LatticeType lattice=LatticeType::FCC)
    {
        Affine::CubicDegradable::Define(a_C11, a_C12, a_C44, m);
        amrex::Vector<Set::Scalar> g_alpha, exp0, gam0dot;
        if(lattice == LatticeType::FCC)
        {
            g_alpha.resize(12); exp0.resize(12); gam0dot.resize(12);
            for (int i=0; i<12; i++)
            {
                g_alpha[i] = 20.0; // default values.
                exp0[i] = 10;
                gam0dot[i] = 0.001;
            }
        }
        else 
            Util::Abort(INFO, "This type of lattice has not been implemented yet");
        DefineSlipSystem(g_alpha, exp0, gam0dot, lattice);
    }
    
    void Define (Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m, amrex::Vector<Set::Scalar> &a_g_alpha, amrex::Vector<Set::Scalar> &a_exp, amrex::Vector<Set::Scalar> &a_gam0dot, LatticeType lattice=LatticeType::FCC)
    {
        Affine::CubicDegradable::Define(a_C11, a_C12, a_C44, m);
        DefineSlipSystem(a_g_alpha, a_exp, a_gam0dot, lattice);
    }

    void Define (Set::Scalar a_C11, Set::Scalar a_C12, Set::Scalar a_C44, Eigen::Matrix3d m, amrex::Vector<SlipSystem> &a_slip)
    {
        Affine::CubicDegradable::Define(a_C11, a_C12, a_C44, m);
        DefineSlipSystem(a_slip);
    }

    void DefineSlipSystem(amrex::Vector<Set::Scalar> &a_g_alpha, amrex::Vector<Set::Scalar> &a_exp, amrex::Vector<Set::Scalar> &a_gam0dot, LatticeType lattice=LatticeType::FCC)
    {
        if(lattice == LatticeType::FCC)
        {
            if(a_g_alpha.size()!=12 || a_exp.size() !=12 || a_gam0dot.size() != 12) 
                Util::Abort(INFO, "Incorrect size of slip system input information");
            
            int size = 12; 
            Set::Matrix R = Set::Matrix::Zero();
            AMREX_D_TERM(
                R(0,0) = m(0,0);
                ,
                R(0,1) = m(0,1); R(1,0) = m(1,0); R(1,1) = m(1,1);
                ,
                R(0,2) = m(0,2); R(1,2) = m(1,2); R(2,0) = m(2,0); R(2,1) = m(2,1); R(2,2) = m(2,2);
            );
            Set::Scalar dir = AMREX_D_PICK(1.0,1.0,INVSQ2);
            Set::Scalar nor = AMREX_D_PICK(1.0, INVSQ2, INVSQ3);
            
            slip.resize(size);
            slip[0].s = R*Set::Vector(AMREX_D_DECL(nor, nor, nor));
            slip[1].s = R*Set::Vector(AMREX_D_DECL(nor, nor, nor));
            slip[2].s = R*Set::Vector(AMREX_D_DECL(nor, nor, nor));
            slip[0].m = R*Set::Vector(AMREX_D_DECL(0.0, -dir, dir));
            slip[1].m = R*Set::Vector(AMREX_D_DECL(dir, 0.0, -dir));
            slip[2].m = R*Set::Vector(AMREX_D_DECL(-dir, dir, 0.0));

            slip[3].s = R*Set::Vector(AMREX_D_DECL(-nor, -nor, nor));
            slip[4].s = R*Set::Vector(AMREX_D_DECL(-nor, -nor, nor));
            slip[5].s = R*Set::Vector(AMREX_D_DECL(-nor, -nor, nor));
            slip[3].m = R*Set::Vector(AMREX_D_DECL(0.0, dir, dir));
            slip[4].m = R*Set::Vector(AMREX_D_DECL(-dir, 0.0, -dir));
            slip[5].m = R*Set::Vector(AMREX_D_DECL(dir, -dir, 0.0));

            slip[6].s = R*Set::Vector(AMREX_D_DECL(-nor, nor, nor));
            slip[7].s = R*Set::Vector(AMREX_D_DECL(-nor, nor, nor));
            slip[8].s = R*Set::Vector(AMREX_D_DECL(-nor, nor, nor));
            slip[6].m = R*Set::Vector(AMREX_D_DECL(0.0, -dir, dir));
            slip[7].m = R*Set::Vector(AMREX_D_DECL(-dir, 0.0, -dir));
            slip[8].m = R*Set::Vector(AMREX_D_DECL(dir, dir, 0.0));

            slip[9].s = R*Set::Vector(AMREX_D_DECL(nor, -nor, nor));
            slip[10].s = R*Set::Vector(AMREX_D_DECL(nor, -nor, nor));
            slip[11].s = R*Set::Vector(AMREX_D_DECL(nor, -nor, nor));
            slip[9].m = R*Set::Vector(AMREX_D_DECL(0.0, dir, dir));
            slip[10].m = R*Set::Vector(AMREX_D_DECL(dir, 0.0, -dir));
            slip[11].m = R*Set::Vector(AMREX_D_DECL(-dir, -dir, 0.0));

            for (int i = 0; i < size; i++)
            {
                slip[i].on = true;
                slip[i].g_alpha = a_g_alpha[i];
                slip[i].g_alpha0 = a_g_alpha[i];
                slip[i].m_exp = a_exp[i];
                slip[i].gam0dot = a_gam0dot[i];
            }
        }
        else
            Util::Abort(INFO, "Slip system not implemented yet");

        prev = PlasticState::Zero();
        curr = PlasticState::Zero();
    }

    void DefineSlipSystem(amrex::Vector<SlipSystem> &a_slip)
    {
        slip = a_slip;

        prev = PlasticState::Zero();
        curr = PlasticState::Zero();
    }

    void DegradeYieldSurface(Set::Scalar eta)
    {
        if (std::isnan(eta) || std::isinf(eta)) Util::Abort(INFO, "Nans or Infs in eta. ", eta);
        if (eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta. ", eta);
        for (int i=0; i<slip.size(); i++)
            slip[i].g_alpha = slip[i].g_alpha0/eta;
    }

    void DegradeYieldSurface(amrex::Vector<Set::Scalar> &eta)
    {
        if(eta.size() ==0 || eta.size() != slip.size()) Util::Abort(INFO, "Inconsistent size of eta. ", eta.size());
        for (int i=0; i<slip.size(); i++)
            slip[i].g_alpha = slip[i].g_alpha0/eta[i];
    }

    Set::Scalar OriginalYieldSurface()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);
            yield += slip[i].g_alpha0*R;
        }
        return yield.norm();
    }

    Set::Scalar YieldSurface()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);
            yield += slip[i].g_alpha*R;
        }
        return yield.norm();
    }

    Set::Scalar OriginalPlasticEnergy()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);
            yield += slip[i].g_alpha0*R;
        }
        return yield.norm()*curr.alpha;
    }

    Set::Scalar OriginalPlasticEnergy(Set::Scalar alpha)
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);
            yield += slip[i].g_alpha0*R;
        }
        return yield.norm()*alpha;
    }

    Set::Scalar PlasticEnergy()
    {
        Set::Matrix yield = Set::Matrix::Zero(), R = Set::Matrix::Zero();
        for(int i=0; i<slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);
            yield += slip[i].g_alpha*R;
        }
        return yield.norm()*curr.alpha;
    }

    void EvolvePlasticStrain(Set::Matrix sigma, Set::Matrix /*strain*/, Set::Scalar dt)
    {
        curr.epsp = prev.epsp;
        Set::Matrix R = Set::Matrix::Zero();

        amrex::Vector<Set::Scalar> gamdot;
        gamdot.resize(slip.size());

        for (int i = 0; i < slip.size(); i++)
        {
            gamdot[i] = 0.0;
            Set::Scalar tau = slip[i].m.dot(sigma*slip[i].s);
            Set::Scalar strength = slip[i].g_alpha;
            if(tau/strength < 1.0) continue;
            gamdot[i] = slip[i].gam0dot*std::pow((tau/strength), slip[i].m_exp);
        }

        for (int i = 0; i < slip.size(); i++)
        {
            for (int j = 0; j < AMREX_SPACEDIM; j++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    R(j,k) = slip[i].m(j)*slip[i].s(k);

            curr.epsp += dt*gamdot[i]*R;
            gam_bar += dt*std::abs(gamdot[i]);
        }

        if (hardening)
        {
            Set::Scalar sech_gambar = 1.0/(std::cosh(((h0-hs)/(gs-g0))*gam_bar));
            Set::Scalar h_gambar = hs + (hs-h0)*sech_gambar*sech_gambar;
            for (int i = 0; i < slip.size(); i++)
            {
                Set::Scalar gdotalpha = 0;
                for (int j = 0; j < slip.size(); j++)
                {
                    Set::Scalar h_ab = ((slip[i].s == slip[j].s) ? 1.0 : q) * h_gambar;
                    gdotalpha += h_ab*std::abs(gamdot[j]);
                }
                slip[i].g_alpha += dt*gdotalpha;
            }
        }
        curr.alpha = sqrt(2.0/3.0)*curr.epsp.norm();
        Cubic::SetF0(curr.epsp);
    }

    void SetPlasticStrains(PlasticState &a_state)
    {
        prev = a_state;
        curr = a_state;
    }

    PlasticState GetPlasticState()
    {
        return curr;
    }
private: // this is to fix an issue with the intel compiler and constexpr
  constexpr static Set::Scalar INVSQ2 = 0.70710678118654752440084436210484903928483593768847403658833986899536623923105351942519376716382078636750692311545614851246241802792536860632206074854996791570661133296375279637789997525057639103028574;
  constexpr static Set::Scalar INVSQ3 = 0.57735026918962576450914878050195745564760175127012687601860232648397767230293334569371539558574952522520871380513556767665664836499965082627055183736479121617603107730076852735599160670036155830775501;
public:
    amrex::Vector<SlipSystem> slip; 
    bool hardening = false;
    Set::Scalar g0, gs, h0, hs, q, gam_bar;
    PlasticState curr, prev;

public:
    static void Parse(CrystalPlasticDegradable & value, IO::ParmParse & pp)
    {
        Affine::Cubic::Parse(value,pp);
        
        std::string lattice;
        LatticeType input_lattice;
        int number_of_slip_systems;
        bool latticeinput = false;

        std::map<std::string, LatticeType> latticemap;
        std::map<LatticeType, int> latticeslipsystemmap;

        amrex::Vector<Set::Scalar> slipplanes, slipdirections;
        amrex::Vector<Set::Vector> input_slipplanes, input_slipdirections;

        latticemap["fcc"] = LatticeType::FCC;
        latticemap["FCC"] = LatticeType::FCC;
        latticemap["bcc"] = LatticeType::BCC;
        latticemap["BCC"] = LatticeType::BCC;
        latticemap["hcp"] = LatticeType::HCP;
        latticemap["HCP"] = LatticeType::HCP;
        
        latticeslipsystemmap[LatticeType::FCC] = 12;
        latticeslipsystemmap[LatticeType::BCC] = 48;
        latticeslipsystemmap[LatticeType::HCP] = 2;
        if(pp.contains("lattice"))
        {
            pp.query("lattice",lattice);
            input_lattice = latticemap[lattice];
            number_of_slip_systems = latticeslipsystemmap[input_lattice];
            latticeinput = true;
        }
        else if (pp.contains("slipplane") && pp.contains("slipdirections"))
        {
            pp.queryarr("slipplane",slipplanes);
            pp.queryarr("slipdirections",slipdirections);
            if(slipplanes.size() <=0 || slipplanes.size()%AMREX_SPACEDIM != 0 || slipdirections.size()%AMREX_SPACEDIM !=0 || slipplanes.size() != slipdirections.size()) Util::Abort(INFO, "Insufficient information on slip system");
            for(int i = 0; i < slipdirections.size()/AMREX_SPACEDIM; i++)
            {
                input_slipdirections.push_back(Set::Vector(AMREX_D_DECL(slipdirections[i],slipdirections[i+1],slipdirections[i+2])));
                input_slipplanes.push_back(Set::Vector(AMREX_D_DECL(slipplanes[i],slipplanes[i+1],slipplanes[i+2])));
            }
            for (int i = 0; i < input_slipdirections.size(); i++)
            {
                if(input_slipdirections[i].norm() ==0 || input_slipplanes[i].norm() == 0) Util::Abort(INFO, "invalid values of slip planes  or directions");
                input_slipdirections[i] = input_slipdirections[i]/input_slipdirections[i].norm();
                input_slipplanes[i] = input_slipplanes[i]/input_slipplanes[i].norm();
            }
            number_of_slip_systems = input_slipdirections.size();
        }
        else
        {
            input_lattice = LatticeType::FCC;
            number_of_slip_systems = latticeslipsystemmap[input_lattice];
            latticeinput = true;
        }

        amrex::Vector<Set::Scalar> input_g_alpha, g_alpha;
        g_alpha.resize(number_of_slip_systems); 

        pp.query("hardening.on", value.hardening);
        if (value.hardening)
        {
            pp.query("hardening.g0",value.g0);
            pp.query("hardening.gs",value.gs);
            pp.query("hardening.h0",value.h0);
            pp.query("hardening.hs",value.hs);
            pp.query("hardening.q",value.q);
            value.gam_bar = 0.0;
            for (int i = 0; i< number_of_slip_systems; i++) g_alpha[i] = value.g0;
        }
        else
        {
            if(pp.contains("g_alpha"))
            {
                pp.queryarr("g_alpha", input_g_alpha);
                if(input_g_alpha.size() == 1)
                {
                    for (int i = 0; i< number_of_slip_systems; i++) 
                        g_alpha[i] = input_g_alpha[0];
                }
                else if(input_g_alpha.size() <= 0 || input_g_alpha.size() != number_of_slip_systems) Util::Abort(INFO, "Invalid number of values in g_alpha");
                else
                {
                    for (int i = 0; i< number_of_slip_systems; i++) 
                        g_alpha[i] = input_g_alpha[i];
                }
            }
            else
                for (int i = 0; i< number_of_slip_systems; i++) g_alpha[i] = 16.0;
        }

        amrex::Vector<Set::Scalar> g_exp, input_exp;
        g_exp.resize(number_of_slip_systems); 
        if(pp.contains("exp"))
        {
            pp.queryarr("exp",input_exp);
            if(input_exp.size() == 1)
            {
                for (int i = 0; i< number_of_slip_systems; i++) 
                    g_exp[i] = input_exp[0];
            }
            else if(input_exp.size() <= 0 || input_exp.size() != number_of_slip_systems) Util::Abort(INFO, "Invalid number of values in exp");
            else
            {
                for (int i = 0; i< number_of_slip_systems; i++) 
                    g_exp[i] = input_exp[i];
            }
        }
        else
            for (int i = 0; i< number_of_slip_systems; i++) g_exp[i] = 10.0;

        amrex::Vector<Set::Scalar> input_gam0dot, gam0dot;
        gam0dot.resize(number_of_slip_systems);
        if(pp.contains("gam0dot"))
        {
            pp.queryarr("gam0dot",input_gam0dot);
            if(input_gam0dot.size() == 1)
            {
                for (int i = 0; i< number_of_slip_systems; i++) 
                    gam0dot[i] = input_gam0dot[0];
            }
            else if(input_gam0dot.size() <= 0 || input_gam0dot.size() != number_of_slip_systems) Util::Abort(INFO, "Invalid number of values in gam0dot");
            else
            {
                for (int i = 0; i< number_of_slip_systems; i++) 
                    gam0dot[i] = input_gam0dot[i];
            }
        }
        else
            for (int i = 0; i< number_of_slip_systems; i++) gam0dot[i] = 0.001;

        if(latticeinput)
        {
            value.DefineSlipSystem(g_alpha,g_exp,gam0dot, input_lattice);
        }
        else
        {
            amrex::Vector<SlipSystem> a_slipsystem;
            a_slipsystem.resize(number_of_slip_systems);
            for (int i = 0; i<number_of_slip_systems; i++)
            {
                a_slipsystem[i].on = true;
                a_slipsystem[i].s = input_slipplanes[i];
                a_slipsystem[i].m = input_slipdirections[i];
                a_slipsystem[i].g_alpha = g_alpha[i];
                a_slipsystem[i].g_alpha0 = g_alpha[i];
                a_slipsystem[i].m_exp = g_exp[i];
                a_slipsystem[i].gam0dot = gam0dot[i];
            }
            value.DefineSlipSystem(a_slipsystem);
        }
    }
    AMREX_FORCE_INLINE
	void operator += (const CrystalPlasticDegradable &rhs)
    {
        curr += rhs.curr;
        prev += rhs.prev;
        for (int i = 0; i< std::min(rhs.slip.size(),slip.size()); i++) slip[i] += rhs.slip[i];
        g0 += rhs.g0; 
        gs += rhs.gs;
        h0 += rhs.h0;
        hs += rhs.hs;
        q += rhs.q;
        gam_bar += rhs.gam_bar;
        ddw += rhs.ddw;
        F0  += rhs.F0;
    }
    AMREX_FORCE_INLINE
    CrystalPlasticDegradable operator * (Set::Scalar alpha)
    {
        CrystalPlasticDegradable ret;
        ret.curr = curr*alpha;
        ret.prev = prev*alpha;
        ret.slip.resize(slip.size());
        for (int i=0; i<slip.size(); i++) ret.slip[i] = slip[i]*alpha;
        ret.g0 = g0*alpha; 
        ret.gs = gs*alpha; 
        ret.h0 = h0*alpha; 
        ret.hs = hs*alpha;
        ret.q = q*alpha;
        ret.gam_bar = gam_bar*alpha;
        ret.hardening = hardening;
        ret.ddw = ddw*alpha;
        ret.F0  = alpha*F0;
    
        return ret;
    }
    friend CrystalPlasticDegradable operator * (const Set::Scalar alpha, const CrystalPlasticDegradable b);
    friend CrystalPlasticDegradable operator + (const CrystalPlasticDegradable a, const CrystalPlasticDegradable b);
    friend CrystalPlasticDegradable operator - (const CrystalPlasticDegradable a, const CrystalPlasticDegradable b);
};
AMREX_FORCE_INLINE
CrystalPlasticDegradable operator * (const Set::Scalar alpha, const CrystalPlasticDegradable b)
{
    CrystalPlasticDegradable ret;
    ret.curr = b.curr*alpha;
    ret.prev = b.prev*alpha;
    ret.slip.resize(b.slip.size());
    for (int i=0; i<b.slip.size(); i++) ret.slip[i] = b.slip[i]*alpha;
    ret.g0 = b.g0*alpha; 
    ret.gs = b.gs*alpha; 
    ret.h0 = b.h0*alpha; 
    ret.hs = b.hs*alpha;
    ret.q = b.q*alpha;
    ret.gam_bar = b.gam_bar*alpha;
    ret.hardening = b.hardening;

    ret.ddw = b.ddw*alpha;
    ret.F0  = alpha*b.F0;
    
    return ret;
}

AMREX_FORCE_INLINE
CrystalPlasticDegradable operator + (const CrystalPlasticDegradable a, const CrystalPlasticDegradable b)
{
    CrystalPlasticDegradable ret;
    ret.curr = a.curr + b.curr;
    ret.prev = a.prev + b.prev;
    ret.slip.resize(std::min(a.slip.size(),b.slip.size()));
    for (int i=0; i< std::min(a.slip.size(),b.slip.size()); i++) ret.slip[i] = a.slip[i] + b.slip[i];

    ret.g0 = a.g0 + b.g0; 
    ret.gs = a.gs + b.gs; 
    ret.h0 = a.h0 + b.h0; 
    ret.hs = a.hs + b.hs;
    ret.q = a.q + b.q;
    ret.gam_bar = a.gam_bar + b.gam_bar;
    ret.hardening = a.hardening && b.hardening;

    ret.ddw = a.ddw + b.ddw;
    ret.F0  = a.F0 + b.F0;
    
    return ret;
}

AMREX_FORCE_INLINE
CrystalPlasticDegradable operator - (const CrystalPlasticDegradable a, const CrystalPlasticDegradable b)
{
    CrystalPlasticDegradable ret;
    ret.curr = a.curr - b.curr;
    ret.prev = a.prev - b.prev;
    ret.slip.resize(std::min(a.slip.size(),b.slip.size()));
    for (int i=0; i< std::min(a.slip.size(),b.slip.size()); i++) ret.slip[i] = a.slip[i] - b.slip[i];

    ret.g0 = a.g0 - b.g0; 
    ret.gs = a.gs - b.gs; 
    ret.h0 = a.h0 - b.h0; 
    ret.hs = a.hs - b.hs;
    ret.q = a.q - b.q;
    ret.gam_bar = a.gam_bar - b.gam_bar;
    ret.hardening = a.hardening && b.hardening;

    ret.ddw = a.ddw - b.ddw;
    ret.F0  = a.F0 - b.F0;
    
    return ret;
}

}
}
}

#endif
