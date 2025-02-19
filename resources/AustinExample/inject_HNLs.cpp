#include <LeptonInjector/Particle.h>
#include <LeptonInjector/LeptonInjector.h>
#include <LeptonInjector/Constants.h>
#include <earthmodel-service/EarthModel.h>
#include <earthmodel-service/Geometry.h>
#include <string>
#include <memory>
#include <chrono>
#include <ctime>
#include <argagg.hpp>
#include "date.h"

#include <cmath>
#include <math.h>
#include <memory>
#include <vector>
#include <array>
#include <iostream>


#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"
#include "LeptonInjector/LeptonInjector.h"

template <class Precision>
std::string getISOCurrentTimestamp() {
    auto now = std::chrono::system_clock::now();
    return date::format("%FT%TZ", date::floor<Precision>(now));
}


using namespace LeptonInjector;

std::string diff_xs(int Z, int A) {
  std::stringstream ss;
  ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/diff_xsec_y_Enu/";
    ss << "dxsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_0.001";
    return ss.str();
}

std::string tot_xs(int Z, int A) {
  std::stringstream ss;
  ss << "/home/nwkamp/Research/Pheno/Neutrissimos2/sources/nu-dipole/xsecs/xsec_tables/tot_xsec_y_Enu/";
    ss << "xsec_";
    ss << "Z_" << Z << "_";
    ss << "A_" << A << "_";
    ss << "mHNL_0.001";
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

std::vector<std::string> gen_diff_xs_hf() {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1]) + "_hf.dat");
    }
}

std::vector<std::string> gen_tot_xs_hf() {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1]) + "_hf.dat");
    }
}

std::vector<std::string> gen_diff_xs_hc() {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(diff_xs(za[0], za[1]) + "_hc.dat");
    }
}

std::vector<std::string> gen_tot_xs_hc() {
    std::vector<std::string> res;
    for(auto const & za : gen_ZA()) {
        res.push_back(tot_xs(za[0], za[1]) + "_hc.dat");
    }
}




int main(int argc, char ** argv) {

		using ParticleType = LeptonInjector::Particle::ParticleType;

    std::string material_file = "";
    std::string earth_file = "";
    double powerLawIndex = 2;
    double energyMin = 1; // in GeV
    double energyMax = 20; // in GeV

    double hnl_mass = 0.001; // in GeV; The HNL mass we are injecting

    // Decay parameters used to set the max range when injecting an HNL, decay width is likely wrong, set accordingly...
    double HNL_decay_mass = 0.001; // in GeV
    double HNL_decay_width = 0.001; // in GeV; decay_width == 1.0 / decay_time; fixme
    double n_decay_lengths = 3.0;

    // This should encompass Minerva, should probably be smaller? Depends on how long Minerva is...
    double disk_radius = 10; // in meters
    double endcap_length = 10; // in meters


    // Events to inject
    unsigned int events_to_inject = 1000;
    Particle::ParticleType primary_type = ParticleType::NuE;

    // Load cross sections
    std::vector<std::shared_ptr<CrossSection>> cross_sections;
    std::vector<Particle::ParticleType> primary_types = {Particle::ParticleType::NuE, Particle::ParticleType::NuMu, Particle::ParticleType::NuTau};
    std::vector<Particle::ParticleType> target_types = gen_TargetPIDs();
    std::shared_ptr<DipoleFromTable> hf_xs = std::make_shared<DipoleFromTable>(hnl_mass);
    std::shared_ptr<DipoleFromTable> hc_xs = std::make_shared<DipoleFromTable>(hnl_mass);
    std::vector<std::string> hf_diff_fnames = gen_diff_xs_hf();
    std::vector<std::string> hc_diff_fnames = gen_diff_xs_hc();
    std::vector<std::string> hf_tot_fnames = gen_tot_xs_hf();
    std::vector<std::string> hc_tot_fnames = gen_tot_xs_hc();
    for(unsigned int i=0; i < target_types.size(); ++i) {
        hf_xs->AddDifferentialCrossSectionFile(hf_diff_fnames[i], target_types[i]);
        hf_xs->AddTotalCrossSectionFile(hf_tot_fnames[i], target_types[i]);
        hc_xs->AddDifferentialCrossSectionFile(hc_diff_fnames[i], target_types[i]);
        hc_xs->AddTotalCrossSectionFile(hc_tot_fnames[i], target_types[i]);
    }
    cross_sections.push_back(hf_xs);
    cross_sections.push_back(hc_xs);

    // Load the earth model
    std::shared_ptr<earthmodel::EarthModel> earth_model = std::make_shared<earthmodel::EarthModel>();
    earth_model->LoadMaterialModel(material_file);
    earth_model->LoadEarthModel(earth_file);

    // Setup power law
    std::shared_ptr<LI_random> random = std::make_shared<LI_random>();
    std::shared_ptr<LeptonInjector::PowerLaw> power_law = std::make_shared<LeptonInjector::PowerLaw>();
    power_law->powerLawIndex = powerLawIndex;
    power_law->energyMin = energyMin;
    power_law->energyMax = energyMax;
    std::shared_ptr<PrimaryEnergyDistribution> edist = power_law;

    // Choose injection direction
    std::shared_ptr<PrimaryDirectionDistribution> ddist = std::make_shared<LeptonInjector::FixedDirection>(earthmodel::Vector3D{1.0, 0.0, 0.0});

    // Targets should be stationary
    std::shared_ptr<LeptonInjector::TargetMomentumDistribution> target_momentum_distribution = std::make_shared<LeptonInjector::TargetAtRest>();

    // Let us inject according to the decay distribution
    std::shared_ptr<RangeFunction> range_func = std::make_shared<LeptonInjector::DecayRangeFunction>(HNL_decay_mass, HNL_decay_width, n_decay_lengths);

    // Put it all together!
    RangedLeptonInjector injector(events_to_inject, primary_type, cross_sections, earth_model, random, edist, ddist, target_momentum_distribution, range_func, disk_radius, endcap_length);

    while(injector) {
        LeptonInjector::InteractionRecord event = injector.GenerateEvent();
    }
}
