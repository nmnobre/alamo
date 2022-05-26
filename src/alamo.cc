#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Elastic/NeoHookean.H"
#include "Model/Solid/Linear/Laplacian.H"
#include "Model/Solid/Affine/J2.H"
#include "Model/Solid/Affine/J2Plastic.H"
#include "Model/Solid/Affine/CrystalPlastic.H"
#include "Model/Solid/Affine/StrainGradientCrystalPlastic.H"

#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mechanics.H"
#include "Integrator/Flame.H"
#include "Integrator/PolymerDegradation.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"
#include "Integrator/ThermoElastic.H"
#include "Integrator/Dynamics.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp.query("alamo.program",program);

    if (program == "microstructure")
    {
        srand(2);
        Integrator::Integrator *pfm = new Integrator::PhaseFieldMicrostructure(pp);
        //Integrator::PhaseFieldMicrostructure pfm;
        pfm->InitData();
        pfm->Evolve();
        delete pfm;
    }
    else if (program == "mechanics")
    {
        Integrator::Integrator *integrator;
        std::string model = "linear.isotropic";
        pp.query("alamo.program.mechanics.model",model);
        if (model == "linear.isotropic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Linear::Isotropic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Linear::Isotropic>*>(integrator));
        }
        else if (model == "affine.isotropic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::Isotropic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::Isotropic>*>(integrator));
        }
        else if (model == "linear.laplacian") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Linear::Laplacian>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Linear::Laplacian>*>(integrator));
        }
        else if (model == "elastic.neohookean") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Elastic::NeoHookean>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Elastic::NeoHookean>*>(integrator));
        }
        else if (model == "affine.j2") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::J2>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::J2>*>(integrator));
        }
        else if (model == "affine.j2plastic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::J2Plastic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::J2Plastic>*>(integrator));
        }
        #if AMREX_SPACEDIM==3
        else if (model == "affine.crystalplastic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::CrystalPlastic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::CrystalPlastic>*>(integrator));
        }
        else if (model == "affine.sgcp") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::StrainGradientCrystalPlastic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::StrainGradientCrystalPlastic>*>(integrator));
        }
        #endif
        else
        {
            Util::Abort(INFO,model," is not a valid model");
        }

        integrator->InitData();
        integrator->Evolve();
        delete integrator;
    }
    else if (program == "eshelby")
    {
        IO::ParmParse pp;
        Integrator::Mechanics<Model::Solid::Affine::Isotropic> eshelby;
        pp.queryclass(eshelby);
        eshelby.InitData();
        eshelby.Evolve();        
    }
    // else if (program == "crystalplastictensiontest")
    // {
    //     Integrator::Integrator *cp_tt = new Integrator::TensionTest<Model::Solid::Affine::CrystalPlastic>();
    //     cp_tt->InitData();
    //     cp_tt->Evolve();        
    //     delete cp_tt;
    // }
    // else if (program == "sgcptensiontest")
    // {
    //     Integrator::Integrator *sgcp_tt = new Integrator::TensionTestSGCP<Model::Solid::Affine::StrainGradientCrystalPlastic>();
    //     sgcp_tt->InitData();
    //     sgcp_tt->Evolve();        
    //     delete sgcp_tt;
    // }
    // else if (program == "j2plastictensiontest")
    // {
    //     Integrator::Integrator *j2_tt = new Integrator::TensionTest<Model::Solid::Affine::J2PlasticDegradable>();
    //     j2_tt->InitData();
    //     j2_tt->Evolve();        
    //     delete j2_tt;
    // }
    else if (program == "finitekinematics")
    {
        Integrator::Integrator *fk = new Integrator::Mechanics<Model::Solid::Elastic::NeoHookean>();
        fk->InitData();
        fk->Evolve();        
        delete fk;
    }
    else if (program == "flame")
    {
        Integrator::Integrator *flame = new Integrator::Flame();
        flame->InitData();
        flame->Evolve();
        delete flame;

    }
    else if (program == "heat")
    {
        IO::ParmParse pp;
        Integrator::HeatConduction heatconduction;
        pp.queryclass(heatconduction);
        heatconduction.InitData();
        heatconduction.Evolve();
    }
    else if (program == "degradation")
    {
        srand(1.0*amrex::ParallelDescriptor::MyProc());
        Integrator::PolymerDegradation model;
        model.InitData();
        model.Evolve();
    }
    else if (program == "fracture")
    {
        IO::ParmParse pp;
        Integrator::Fracture model;
        pp.queryclass(model);
        model.InitData();
        model.Evolve();
        //delete model;
    }
    else if (program == "thermoelastic")
    {
        IO::ParmParse pp;
        Integrator::ThermoElastic te;
        pp.queryclass(te);
        te.InitData();
        te.Evolve();
    }
    else if (program == "dynamics")
    {
        IO::ParmParse pp;
        Integrator::Dynamics integrator;
        pp.queryclass(integrator);
        integrator.InitData();
        integrator.Evolve();
    }
    else
    {
        Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");
    }

    Util::Finalize();
} 
