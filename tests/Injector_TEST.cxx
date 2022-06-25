
#include <ios>
#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>
#include <fstream>

#include <gtest/gtest.h>

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Constants.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/LeptonInjector.h"
#include "LeptonInjector/Weighter.h"

#include "earthmodel-service/Geometry.h"
#include "earthmodel-service/EulerQuaternionConversions.h"
#include "earthmodel-service/Placement.h"

//#define AUSTIN

using namespace LeptonInjector;
bool z_samp = true;
bool in_invGeV = true;
bool inelastic = true;
bool miniboone = true;


std::string diff_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
#ifdef AUSTIN
    ss << "/home/austin/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
#else
    if(z_samp) ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/diff_xsec_z_Enu/";
    else ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
#endif
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::string tot_xs(int Z, int A, std::string mHNL) {
    std::stringstream ss;
#ifdef AUSTIN
    ss << "/home/austin/nu-dipole/xsecs/xsec_tables/tot_xsec_Enu/";
#else
    ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/xsec_tables_dipoleFF/tot_xsec_Enu/";
#endif
    ss << "xsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_" << mHNL;
    return ss.str();
}

std::vector<std::array<int, 2>> gen_ZA() {
    return std::vector<std::array<int, 2>>{
        {1, 1},
				{6, 12},
				{8, 16},
				{13, 27},
				{14, 28},
				{20, 40},
				{26, 56},
				{29, 63},
				{29, 65},
				{82, 208},
    };
}

std::vector<LeptonInjector::Particle::ParticleType> gen_TargetPIDs() {
    using ParticleType = LeptonInjector::Particle::ParticleType;
    return std::vector<ParticleType>{
        ParticleType::HNucleus,
            ParticleType::C12Nucleus,
            ParticleType::O16Nucleus,
            ParticleType::Al27Nucleus,
            ParticleType::Si28Nucleus,
            ParticleType::Ca40Nucleus,
            ParticleType::Fe56Nucleus,
            ParticleType::Cu63Nucleus,
            ParticleType::Cu65Nucleus,
            ParticleType::Pb208Nucleus
    };
}

std::vector<std::string> gen_diff_xs_hf(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hf(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL) + "_hf.dat");
    }
    return res;
}

std::vector<std::string> gen_diff_xs_hc(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

std::vector<std::string> gen_tot_xs_hc(std::string mHNL) {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1], mHNL) + "_hc.dat");
    }
    return res;
}

bool inFiducial(std::array<double,3> & int_vtx, earthmodel::ExtrPoly & fidVol) {
    earthmodel::Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    earthmodel::Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

bool inFiducial(std::array<double,3> & int_vtx, earthmodel::Sphere & fidVol) {
    earthmodel::Vector3D pos(int_vtx[0], int_vtx[1], int_vtx[2]);
    earthmodel::Vector3D dir(0,0,1);
    return fidVol.IsInside(pos,dir);
}

double ComputeInteractionLengths(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<LeptonInjector::CrossSectionCollection const> cross_sections, std::pair<earthmodel::Vector3D, earthmodel::Vector3D> const & bounds, InteractionRecord const & record) {
    earthmodel::Vector3D interaction_vertex = record.interaction_vertex;
    earthmodel::Vector3D direction(
            record.primary_momentum[1],
            record.primary_momentum[2],
            record.primary_momentum[3]);
    direction.normalize();

    earthmodel::Geometry::IntersectionList intersections = earth_model->GetIntersections(earth_model->GetEarthCoordPosFromDetCoordPos(interaction_vertex), direction);
	std::map<Particle::ParticleType, std::vector<std::shared_ptr<CrossSection>>> const & cross_sections_by_target = cross_sections->GetCrossSectionsByTarget();
    std::vector<double> total_cross_sections;
    std::vector<LeptonInjector::Particle::ParticleType> targets;
	InteractionRecord fake_record = record;
	for(auto const & target_xs : cross_sections_by_target) {
        targets.push_back(target_xs.first);
		fake_record.target_mass = earth_model->GetTargetMass(target_xs.first);
		fake_record.target_momentum = {fake_record.target_mass,0,0,0};
		std::vector<std::shared_ptr<CrossSection>> const & xs_list = target_xs.second;
		double total_xs = 0.0;
		for(auto const & xs : xs_list) {
			std::vector<InteractionSignature> signatures = xs->GetPossibleSignaturesFromParents(record.signature.primary_type, target_xs.first);
			for(auto const & signature : signatures) {
				fake_record.signature = signature;
				// Add total cross section
				total_xs += xs->TotalCrossSection(fake_record);
			}
		}
		total_cross_sections.push_back(total_xs);
	}
    std::vector<double> particle_depths = earth_model->GetParticleColumnDepth(intersections, bounds.first, bounds.second, targets);

    double interaction_depth = 0.0;
    for(unsigned int i=0; i<targets.size(); ++i) {
        interaction_depth += particle_depths[i] * total_cross_sections[i];
    }
    return interaction_depth;
}

std::vector<double> p_LE_FHC_nue = {1.94e+00, 9.57e-01, 3.86e-01, 1.38e+01, 1.41e-01};
std::vector<double> p_LE_FHC_numu = {2.06e+00, 6.52e-01, 3.36e-01, 7.50e+00, 1.19e-01};
std::vector<double> p_LE_FHC_nuebar = {1.80e+00, 2.95e+00, 3.80e-01, 2.19e+01, 3.12e-01};
std::vector<double> p_LE_FHC_numubar = {2.75e+00, 2.46e+00, 4.90e-01, 4.44e+00, 5.15e-02};

std::vector<double> p_LE_RHC_nue = {2.25e+00, 4.38e+00, 5.82e-01, 2.33e+02, 1.00e+00};
std::vector<double> p_LE_RHC_numu = {3.75e+00, 3.04e+00, 5.53e-01, 1.50e+02, 3.11e-14};
std::vector<double> p_LE_RHC_nuebar = {1.89e+00, 9.06e-01, 3.95e-01, 8.79e+00, 1.02e-01};
std::vector<double> p_LE_RHC_numubar = {1.95e+00, 6.09e-01, 3.49e-01, 5.74e+00, 8.92e-02};

std::vector<double> p_ME_FHC_numu = {4.65e+00, 1.35e+00, 7.24e-02, 3.07e+00, 4.45e-03};


TEST(Injector, Generation)
{
    using ParticleType = LeptonInjector::Particle::ParticleType;

#ifdef AUSTIN
    std::string material_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    std::string earth_file = "/home/austin/programs/LIDUNE/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    std::string flux_file = "/home/austin/nu-dipole/fluxes/LE_FHC_numu.txt";
    z_samp = false;
    in_invGeV = false;
#else
    std::string material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/materials/Minerva.dat";
    std::string earth_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_minerva.dat";
    std::string flux_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/NUMI_Flux_Tables/ME_FHC_numu.txt";
    if(miniboone) {
			material_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/materials/MiniBooNE.dat";
			earth_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/LeptonInjectorDUNE/resources/earthparams/densities/PREM_miniboone.dat";
			flux_file = "/home/nwkamp/Research/Pheno/Neutrissimos2/Sandbox/BNB_Flux_Tables/BNB_numu_flux.txt";
			inelastic = true;
    }
#endif

    double hnl_mass = 0.4; // in GeV; The HNL mass we are injecting
    double dipole_coupling = 1.0e-6; // in GeV^-1; the effective dipole coupling strength
    std::string mHNL = "0.4";

    // Decay parameters used to set the max range when injecting an HNL
    double HNL_decay_width = std::pow(dipole_coupling,2)*std::pow(hnl_mass,3)/(4*Constants::pi); // in GeV; decay_width = d^2 m^3 / (4 * pi)
    double n_decay_lengths = 3.0; // Number of decay lengths to consider
    double max_distance = 240; // Maximum distance, set by distance from MiniBooNE to the decay pipe
    // This should encompass Minerva, should probably be smaller? Depends on how long Minerva is...
    double disk_radius = 1.4; // in meters
    double endcap_length = 5; // in meters
    if(miniboone) {
			max_distance = 541;
			disk_radius = 6.2;
			endcap_length = 6.2;
    }


    // Events to inject
    unsigned int events_to_inject = 5e5;
    Particle::ParticleType primary_type = ParticleType::NuMu;

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = gen_TargetPIDs();
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass, dipole_coupling, DipoleFromTable::HelicityChannel::Flipping, z_samp, in_invGeV, inelastic);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass, dipole_coupling, DipoleFromTable::HelicityChannel::Conserving, z_samp, in_invGeV, inelastic);
    std::vector<std::string> hf_diff_fnames = gen_diff_xs_hf(mHNL);
    std::vector<std::string> hc_diff_fnames = gen_diff_xs_hc(mHNL);
    std::vector<std::string> hf_tot_fnames = gen_tot_xs_hf(mHNL);
    std::vector<std::string> hc_tot_fnames = gen_tot_xs_hc(mHNL);
    for(unsigned int i=0; i < target_types.size(); ++i) {
        //std::cerr << hf_diff_fnames[i] << std::endl;
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        ////std::cerr << hf_tot_fnames[i] << std::endl;
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        //std::cerr << hc_diff_fnames[i] << std::endl;
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        //std::cerr << hc_tot_fnames[i] << std::endl;
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    // Load the earth model
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->LoadMaterialModel(material_file);
    earth_model->LoadEarthModel(earth_file);

    // Setup the primary type and mass
    //std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type, hnl_mass);
    std::shared_ptr<LeptonInjector::PrimaryInjector> primary_injector = std::make_shared<LeptonInjector::PrimaryInjector>(primary_type, 0);

    // Setup power law
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();

    std::vector<double> moyal_exp_params = p_ME_FHC_numu;

    // Setup NUMI flux
    std::shared_ptr<LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution> pdf = std::make_shared<LeptonInjector::ModifiedMoyalPlusExponentialEnergyDistribution>(hnl_mass, 20, moyal_exp_params[0], moyal_exp_params[1], moyal_exp_params[2], moyal_exp_params[3], moyal_exp_params[4]);

    // Setup tabulated flux
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(flux_file, true);
    std::shared_ptr<LeptonInjector::TabulatedFluxDistribution> tab_pdf_gen = std::make_shared<LeptonInjector::TabulatedFluxDistribution>(hnl_mass, 10, flux_file);

    // Change the flux units from cm^-2 to m^-2
    std::shared_ptr<LeptonInjector::WeightableDistribution> flux_units = std::make_shared<LeptonInjector::NormalizationConstant>(1e4);

    // Pick energy distribution
    std::shared_ptr<PrimaryEnergyDistribution> edist = tab_pdf_gen;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::FixedDirection>(earthmodel::Vector3D{0.0, 0.0, 1.0});

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to the decay distribution
    std::shared_ptr<RangeFunction> range_func = std::make_shared<LeptonInjector::DecayRangeFunction>(hnl_mass, HNL_decay_width, n_decay_lengths, max_distance);

    // Helicity distribution
    std::shared_ptr<PrimaryNeutrinoHelicityDistribution> helicity_distribution = std::make_shared<LeptonInjector::PrimaryNeutrinoHelicityDistribution>();

    // Put it all together!
    //RangedLeptonInjector injector(events_to_inject, primary_type, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length);
    std::shared_ptr<InjectorBase> injector = std::make_shared<RangedLeptonInjector>(events_to_inject, primary_injector, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length, helicity_distribution);

    std::vector<std::shared_ptr<WeightableDistribution>> physical_distributions = {
        std::shared_ptr<WeightableDistribution>(tab_pdf),
        std::shared_ptr<WeightableDistribution>(flux_units),
        std::shared_ptr<WeightableDistribution>(ddist),
        std::shared_ptr<WeightableDistribution>(target_momentum_distribution),
        std::shared_ptr<WeightableDistribution>(helicity_distribution)
    };

    LeptonWeighter weighter(std::vector<std::shared_ptr<InjectorBase>>{injector}, earth_model, injector->GetCrossSections(), physical_distributions);

    // MINERvA Fiducial Volume
    std::vector<std::vector<double>> poly;
    poly.push_back({0.0, 1.01758});
    poly.push_back({0.88125, 0.50879});
    poly.push_back({0.88125, -0.50879});
    poly.push_back({0.0, -1.01758});
    poly.push_back({-0.88125, -0.50879});
    poly.push_back({-0.88125, 0.50879});

    double offset[2];
    offset[0] = 0;
    offset[1] = 0;
    std::vector<earthmodel::ExtrPoly::ZSection> zsecs;
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(-2.0672,offset,1));
    zsecs.push_back(earthmodel::ExtrPoly::ZSection(2.0672,offset,1));
    earthmodel::Placement placement(earthmodel::Vector3D(0,0,2.0672), earthmodel::QFromZXZr(0,0,0));
    earthmodel::ExtrPoly MINERvA_fiducial = earthmodel::ExtrPoly(placement, poly, zsecs);
    
    // MiniBooNE Fiducial Volume
    earthmodel::Placement placementMB(earthmodel::Vector3D(0,0,0), earthmodel::QFromZXZr(0,0,0));
    earthmodel::Sphere MiniBooNE_fiducial = earthmodel::Sphere(placement, 5.0, 0.0);

    std::ofstream myFile("injector_test_events.csv");
    // myFile << std::fixed << std::setprecision(6);
    myFile << std::scientific << std::setprecision(16);
    myFile << "intX intY intZ ";
    myFile << "decX decY decZ ";
    myFile << "ppX ppY ppZ ";
    myFile << "p4nu_0 p4nu_1 p4nu_2 p4nu_3 ";
    myFile << "helnu ";
    myFile << "p4itgt_0 p4itgt_1 p4itgt_2 p4itgt_3 ";
    myFile << "helitgt ";
    myFile << "p4hnl_0 p4hnl_1 p4hnl_2 p4hnl_3 ";
    myFile << "helhnl ";
    myFile << "p4ftgt_0 p4ftgt_1 p4ftgt_2 p4ftgt_3 ";
    myFile << "helftgt ";
    myFile << "p4gamma_0 p4gamma_1 p4gamma_2 p4gamma_3 ";
    myFile << "p4gamma_hnlRest_0 p4gamma_hnlRest_1 p4gamma_hnlRest_2 p4gamma_hnlRest_3 ";
    myFile << "helgamma ";
    myFile << "decay_length decay_fid_weight decay_ang_weight prob_nopairprod basic_weight simplified_weight interaction_lengths interaction_prob y target fid\n";
    myFile << std::endl;
    int i = 0;
    while(*injector) {
        LeptonInjector::InteractionRecord event = injector->GenerateEvent();
        LeptonInjector::DecayRecord decay;
        LeptonInjector::InteractionRecord pair_prod;
        double basic_weight, simplified_weight, interaction_lengths, interaction_prob = 0;
        if(event.signature.target_type != LeptonInjector::Particle::ParticleType::unknown) {
            injector->SampleSecondaryDecay(event, decay, HNL_decay_width, 1, 0, &MiniBooNE_fiducial, 0.1);
            injector->SamplePairProduction(decay, pair_prod);
            //basic_weight = weighter.EventWeight(event);
            simplified_weight = weighter.SimplifiedEventWeight(event);
            interaction_lengths = ComputeInteractionLengths(earth_model, injector->GetCrossSections(), injector->InjectionBounds(event), event);
            interaction_prob = weighter.InteractionProbability(injector->InjectionBounds(event), event);
        }
        if(event.secondary_momenta.size() > 0) {
            myFile << event.interaction_vertex[0] << " ";
            myFile << event.interaction_vertex[1] << " ";
            myFile << event.interaction_vertex[2] << " ";

            myFile << decay.decay_vertex[0] << " ";
            myFile << decay.decay_vertex[1] << " ";
            myFile << decay.decay_vertex[2] << " ";

            myFile << pair_prod.interaction_vertex[0] << " ";
            myFile << pair_prod.interaction_vertex[1] << " ";
            myFile << pair_prod.interaction_vertex[2] << " ";

            myFile << event.primary_momentum[0] << " ";
            myFile << event.primary_momentum[1] << " ";
            myFile << event.primary_momentum[2] << " ";
            myFile << event.primary_momentum[3] << " ";

            myFile << event.primary_helicity << " ";

            myFile << event.target_momentum[0] << " ";
            myFile << event.target_momentum[1] << " ";
            myFile << event.target_momentum[2] << " ";
            myFile << event.target_momentum[3] << " ";

            myFile << event.target_helicity << " ";

            myFile << event.secondary_momenta[0][0] << " ";
            myFile << event.secondary_momenta[0][1] << " ";
            myFile << event.secondary_momenta[0][2] << " ";
            myFile << event.secondary_momenta[0][3] << " ";

            myFile << event.secondary_helicity[0] << " ";

            myFile << event.secondary_momenta[1][0] << " ";
            myFile << event.secondary_momenta[1][1] << " ";
            myFile << event.secondary_momenta[1][2] << " ";
            myFile << event.secondary_momenta[1][3] << " ";

            myFile << event.secondary_helicity[1] << " ";

            myFile << decay.secondary_momenta[0][0] << " ";
            myFile << decay.secondary_momenta[0][1] << " ";
            myFile << decay.secondary_momenta[0][2] << " ";
            myFile << decay.secondary_momenta[0][3] << " ";

            myFile << decay.secondary_momenta[1][0] << " ";
            myFile << decay.secondary_momenta[1][1] << " ";
            myFile << decay.secondary_momenta[1][2] << " ";
            myFile << decay.secondary_momenta[1][3] << " ";

            myFile << decay.secondary_helicity[0] << " ";

            myFile << decay.decay_parameters[0] << " "; // decay length
            myFile << decay.decay_parameters[1] << " "; // decay fid weight
            myFile << decay.decay_parameters[2] << " "; // decay anglular weight
            myFile << pair_prod.interaction_parameters[0] << " "; // probability of no pair production
            myFile << basic_weight << " ";
            myFile << simplified_weight << " ";
            myFile << interaction_lengths << " ";
            myFile << interaction_prob << " ";
            myFile << event.interaction_parameters[1] << " "; // sampled y
            myFile << event.signature.target_type << " "; // target type
            if(miniboone) myFile << int(inFiducial(pair_prod.interaction_vertex, MiniBooNE_fiducial)) << "\n"; // fid vol
            else myFile << int(inFiducial(pair_prod.interaction_vertex, MINERvA_fiducial)) << "\n"; // fid vol
            myFile << "\n";
        }
        if((++i) % (events_to_inject/10)==0)
            std::cerr << (int)(i*100 / (events_to_inject)) << "%" << std::endl;
    }
    myFile.close();
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

