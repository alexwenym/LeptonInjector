#include <set>
#include <vector>
#include <string>
#include <memory>
#include <algorithm>

#include <rk/rk.hh>
#include <rk/geom3.hh>

#include "LeptonInjector/dataclasses/Particle.h"

#include "LeptonInjector/utilities/Random.h"
#include "LeptonInjector/utilities/Integration.h"

#include "LeptonInjector/crosssections/CrossSection.h"
#include "LeptonInjector/crosssections/ElasticScattering.h"

namespace LI {
namespace crosssections {

bool ElasticScattering::equal(CrossSection const & other) const {
    const ElasticScattering* x = dynamic_cast<const ElasticScattering*>(&other);

    if(!x)
        return false;
    else
        return primary_types == x->primary_types;
}

double ElasticScattering::InteractionThreshold(dataclasses::InteractionRecord const & record) const {
    return 0.0;
}

double ElasticScattering::DifferentialCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    double CLL;
    if(primary_type==LI::dataclasses::Particle::ParticleType::NuE) CLL = 0.7276;
    else if(primary_type==LI::dataclasses::Particle::ParticleType::NuMu) CLL = -0.2730;
    else {
        std::cout << "Faulty primary: " << primary_type << std::endl;
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double s = std::pow(rk::invMass(p1, p2), 2);
    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
        p1_lab = p1;
        p2_lab = p2;
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }
    assert(interaction.signature.secondary_types.size() == 2);
    assert(interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuE or interaction.signature.secondary_types[1] == LI::dataclasses::Particle::ParticleType::NuE or interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuMu or interaction.signature.secondary_types[1] == LI::dataclasses::Particle::ParticleType::NuMu);
    unsigned int nu_index = (interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuE or interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuMu) ? 0 : 1;
    unsigned int electron_index = 1 - nu_index;

    std::array<double, 4> const & mom3 = interaction.secondary_momenta[nu_index];
    std::array<double, 4> const & mom4 = interaction.secondary_momenta[electron_index];
    rk::P4 p3(geom3::Vector3(mom3[1], mom3[2], mom3[3]), interaction.secondary_masses[nu_index]);
    rk::P4 p4(geom3::Vector3(mom4[1], mom4[2], mom4[3]), interaction.secondary_masses[electron_index]);

    double y = 1.0 - p2.dot(p3) / p2.dot(p1);
    double m = interaction.secondary_masses[electron_index];
    double E = primary_energy;

    // use tree level result
    double term1 = CLL*CLL;// * (1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X1);
    double term2 = CLR*CLR * (1-y)*(1-y);// * ( 1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X2);
    double term3 = -CLL*CLR*m*y/E;// * (1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X3);
    double ret =  std::pow(LI::utilities::Constants::FermiConstant,2) * s / LI::utilities::Constants::pi * (term1 + term2 + term3) / LI::utilities::Constants::invGeVsq_per_cmsq;
    return std::max(ret,0.);
}

// Assume initial electron at rest
double ElasticScattering::DifferentialCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, double y) const {
    double CLL;
    if(primary_type==LI::dataclasses::Particle::ParticleType::NuE) CLL = 0.7276;
    else if(primary_type==LI::dataclasses::Particle::ParticleType::NuMu) CLL = -0.2730;
    else {
        std::cout << "Faulty primary: " << primary_type << std::endl;
        throw std::runtime_error("Supplied primary not supported by cross section!");
    }

    double m = LI::utilities::Constants::electronMass;
    double E = primary_energy;
    double s = 2*m*E + m*m;

    double term1 = CLL*CLL;// * (1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X1);
    double term2 = CLR*CLR * (1-y)*(1-y);// * ( 1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X2);
    double term3 = -CLL*CLR*m*y/E;// * (1 + LI::utilities::Constants::fineStructure / LI::utilities::Constants::pi * X3);
    double ret = std::pow(LI::utilities::Constants::FermiConstant,2) * s / LI::utilities::Constants::pi * (term1 + term2 + term3) / LI::utilities::Constants::invGeVsq_per_cmsq;
    return std::max(ret,0.);
}

double ElasticScattering::TotalCrossSection(dataclasses::InteractionRecord const & interaction) const {
    LI::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;
    LI::dataclasses::Particle::ParticleType target_type = interaction.signature.target_type;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        rk::P4 p1_lab = boost_start_to_lab * p1;
        primary_energy = p1_lab.e();
    }
    // if we are below threshold, return 0
    if(primary_energy < InteractionThreshold(interaction))
        return 0;
    return TotalCrossSection(primary_type, primary_energy, target_type);
}

double ElasticScattering::TotalCrossSection(LI::dataclasses::Particle::ParticleType primary_type, double primary_energy, LI::dataclasses::Particle::ParticleType target_type) const {
    double ymax = 2*primary_energy / (2*primary_energy + LI::utilities::Constants::electronMass);
    std::function<double(double)> integrand = [&] (double y) -> double {
        return DifferentialCrossSection(primary_type, primary_energy, y);
    };
    return LI::utilities::rombergIntegrate(integrand, 0, ymax);
}

void ElasticScattering::SampleFinalState(dataclasses::InteractionRecord& interaction, std::shared_ptr<LI::utilities::LI_random> random) const {

    // Uses Metropolis-Hastings Algorithm!
    // useful for cases where we don't know the supremum of our distribution, and the distribution is multi-dimensional

    LI::dataclasses::Particle::ParticleType primary_type = interaction.signature.primary_type;

    unsigned int nu_index = (interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuE or interaction.signature.secondary_types[0] == LI::dataclasses::Particle::ParticleType::NuMu) ? 0 : 1;
    unsigned int electron_index = 1 - nu_index;
    rk::P4 p1(geom3::Vector3(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]), interaction.primary_mass);
    rk::P4 p2(geom3::Vector3(interaction.target_momentum[1], interaction.target_momentum[2], interaction.target_momentum[3]), interaction.target_mass);
    double primary_energy;
    rk::P4 p1_lab;
    rk::P4 p2_lab;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        primary_energy = interaction.primary_momentum[0];
        p1_lab = p1;
        p2_lab = p2;
    } else {
        rk::Boost boost_start_to_lab = p2.restBoost();
        p1_lab = boost_start_to_lab * p1;
        p2_lab = boost_start_to_lab * p2;
        primary_energy = p1_lab.e();
    }

    double yMin = 1e-15;
    double yMax = 2*primary_energy / (2*primary_energy + LI::utilities::Constants::electronMass);
    assert(yMin > 0);
    double log_yMax = log10(yMax);
    double log_yMin = log10(yMin);

    bool accept;

    double y = 0, test_y;

    // values of cross_section from the table
    double cross_section, test_cross_section;

    // sample an intial point
    test_y = std::pow(10.0,random->Uniform(log_yMin, log_yMax));

    // evalutates the differential spline at that point
    test_cross_section = ElasticScattering::DifferentialCrossSection(primary_type, primary_energy, test_y);

    cross_section = test_cross_section;

    // this is the magic part. Metropolis Hastings Algorithm.
    // MCMC method!
    const size_t burnin = 40; // converges to the correct distribution over multiple samplings.
                              // big number means more accurate, but slower
    for(size_t j = 0; j <= burnin; j++) {
        // repeat the sampling from above to get a new valid point
        test_y = std::pow(10.0,random->Uniform(log_yMin, log_yMax));

        // Load the differential cross section depending on sampling variable
        test_cross_section = ElasticScattering::DifferentialCrossSection(primary_type, primary_energy, test_y);
        if(std::isnan(test_cross_section) or test_cross_section <= 0)
            continue;

        //double odds = ((test_y*test_cross_section) / (y*cross_section));
        double odds = (test_cross_section / cross_section); // this gives a better match to the y distribution
        accept = (cross_section == 0 || (odds > 1.) || random->Uniform(0, 1) < odds);

        if(accept) {
            y  = test_y;
            cross_section = test_cross_section;
        }
    }
    double final_y = y + 1e-16; // to account for machine epsilon when adding to O(1) numbers

    interaction.interaction_parameters.resize(1);
    interaction.interaction_parameters[0] = final_y;

    geom3::UnitVector3 x_dir = geom3::UnitVector3::xAxis();
    geom3::Vector3 p1_mom = p1_lab.momentum();
    geom3::UnitVector3 p1_lab_dir = p1_mom.direction();
    geom3::Rotation3 x_to_p1_lab_rot = geom3::rotationBetween(x_dir, p1_lab_dir);

    double phi = random->Uniform(0, 2.0 * M_PI);
    geom3::Rotation3 rand_rot(p1_lab_dir, phi);

    double m3 = LI::utilities::Constants::electronMass;

    double E3_lab = primary_energy * (final_y) + m3;
    double p3_lab_sq = E3_lab * E3_lab- m3 * m3;
    double costh_lab = (m3+primary_energy)/primary_energy * sqrt((E3_lab - m3)/(E3_lab + m3));
    //std::cout << final_y << " " << costh_lab << std::endl;
    double p3x_lab = costh_lab * sqrt(p3_lab_sq);
    double p3y_lab = sqrt(p3_lab_sq - p3x_lab * p3x_lab);

    rk::P4 p3_lab(geom3::Vector3(p3x_lab, p3y_lab, 0), m3);
    p3_lab.rotate(x_to_p1_lab_rot);
    p3_lab.rotate(rand_rot);
    // doing something dumb, ignore outgoing neutrino
    rk::P4 p4_lab = p1_lab;//p2_lab + (p1_lab - p3_lab);

    rk::P4 p3;
    rk::P4 p4;
    if(interaction.target_momentum[1] == 0 and interaction.target_momentum[2] == 0 and interaction.target_momentum[3] == 0) {
        p3 = p3_lab;
        p4 = p4_lab;
    } else {
        rk::Boost boost_lab_to_start = p2.labBoost();
        p3 = boost_lab_to_start * p3_lab;
        p4 = boost_lab_to_start * p4_lab;
    }

    interaction.secondary_momenta.resize(2);
    interaction.secondary_masses.resize(2);
    interaction.secondary_helicity.resize(2);

    interaction.secondary_momenta[electron_index][0] = p3.e(); // p3_energy
    interaction.secondary_momenta[electron_index][1] = p3.px(); // p3_x
    interaction.secondary_momenta[electron_index][2] = p3.py(); // p3_y
    interaction.secondary_momenta[electron_index][3] = p3.pz(); // p3_z
    interaction.secondary_masses[electron_index] = p3.m();


    interaction.secondary_helicity[electron_index] = interaction.target_helicity;

    interaction.secondary_momenta[nu_index][0] = p4.e(); // p4_energy
    interaction.secondary_momenta[nu_index][1] = p4.px(); // p4_x
    interaction.secondary_momenta[nu_index][2] = p4.py(); // p4_y
    interaction.secondary_momenta[nu_index][3] = p4.pz(); // p4_z
    interaction.secondary_masses[nu_index] = p4.m();

    interaction.secondary_helicity[nu_index] = interaction.primary_helicity;
}

std::vector<LI::dataclasses::Particle::ParticleType> ElasticScattering::GetPossibleTargets() const {
    std::vector<LI::dataclasses::Particle::ParticleType> res;
    res.push_back(LI::dataclasses::Particle::ParticleType::EMinus);
    return res;
}

std::vector<LI::dataclasses::Particle::ParticleType> ElasticScattering::GetPossibleTargetsFromPrimary(LI::dataclasses::Particle::ParticleType primary_type) const {
    if(not primary_types.count(primary_type)) {
        return std::vector<LI::dataclasses::Particle::ParticleType>();
    }
    return GetPossibleTargets();
}

std::vector<LI::dataclasses::Particle::ParticleType> ElasticScattering::GetPossiblePrimaries() const {
    return std::vector<LI::dataclasses::Particle::ParticleType>(primary_types.begin(), primary_types.end());
}

std::vector<dataclasses::InteractionSignature> ElasticScattering::GetPossibleSignatures() const {
    std::vector<LI::dataclasses::Particle::ParticleType> targets = GetPossibleTargets();
    std::vector<dataclasses::InteractionSignature> signatures;
    dataclasses::InteractionSignature signature;
    signature.secondary_types.resize(2);

    for(auto primary : primary_types) {
        signature.primary_type = primary;
        signature.secondary_types[0] = primary;
        for(auto target : targets) {
            signature.target_type = target;
            signature.secondary_types[1] = target;
            signatures.push_back(signature);
        }
    }
    return signatures;
}

std::vector<dataclasses::InteractionSignature> ElasticScattering::GetPossibleSignaturesFromParents(LI::dataclasses::Particle::ParticleType primary_type, LI::dataclasses::Particle::ParticleType target_type) const {
    std::vector<LI::dataclasses::Particle::ParticleType> targets = GetPossibleTargets();
    if(primary_types.count(primary_type) > 0 and std::find(targets.begin(), targets.end(), target_type) != targets.end()) {
        dataclasses::InteractionSignature signature;
        signature.secondary_types.resize(2);
        signature.primary_type = primary_type;
        signature.target_type = target_type;
        signature.secondary_types[1] = target_type;
        if(primary_types.count(primary_type))
            signature.secondary_types[0] = primary_type;
        else
            throw std::runtime_error("Primary type not in primary_types!");
        return std::vector<dataclasses::InteractionSignature>{signature};
    } else {
        return std::vector<dataclasses::InteractionSignature>();
    }
}

double ElasticScattering::FinalStateProbability(dataclasses::InteractionRecord const & interaction) const {
    double dxs = DifferentialCrossSection(interaction);
    double txs = TotalCrossSection(interaction);
    if(dxs == 0) {
        return 0.0;
    } else if (txs == 0) {
        return 0.0;
    } else {
        return dxs / txs;
    }
}

std::vector<std::string> ElasticScattering::DensityVariables() const {
    return std::vector<std::string>{"Bjorken y"};
}

} // namespace crosssections
} // namespace LI

