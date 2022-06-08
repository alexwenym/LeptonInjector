
#include "earthmodel-service/Path.h"
#include "earthmodel-service/Vector3D.h"
#include "earthmodel-service/EarthModel.h"
#include "earthmodel-service/EarthModelCalculator.h"

#include "phys-services/CrossSection.h"

#include "LeptonInjector/Random.h"
#include "LeptonInjector/Particle.h"

#include "LeptonInjector/Distributions.h"

namespace LeptonInjector {

namespace {
    double log_one_minus_exp_of_negative(double x) {
        if(x < 1e-1) {
            return std::log(x) - x/2.0 + x*x/24.0 - x*x*x*x/2880.0;
        } else if(x > 3) {
            double ex = std::exp(-x);
            double ex2 = ex * ex;
            double ex3 = ex2 * ex;
            double ex4 = ex3 * ex;
            double ex5 = ex4 * ex;
            double ex6 = ex5 * ex;
            return -(ex + ex2 / 2.0 + ex3 / 3.0 + ex4 / 4.0 + ex5 / 5.0 + ex6 / 6.0);
        } else {
            return std::log(1.0 - std::exp(-x));
        }
    }
    bool fexists(const std::string filename)
    {
            std::ifstream ifile(filename.c_str());
            return (bool)ifile;
    }

}

//---------------
// PhysicallyNormalizedDistribution
//---------------
PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution() {
    SetNormalization(1.0);
}

PhysicallyNormalizedDistribution::PhysicallyNormalizedDistribution(double norm) {
    SetNormalization(norm);
}

void PhysicallyNormalizedDistribution::SetNormalization(double norm) {
    normalization = norm;
}

double PhysicallyNormalizedDistribution::GetNormalization() const {
    return normalization;
}

bool PhysicallyNormalizedDistribution::IsNormalizationSet() const {
    return normalization != 1.0;
}


//---------------
// class WeightableDistribution
//---------------
std::vector<std::string> WeightableDistribution::DensityVariables() const {
    return {};
}

bool WeightableDistribution::operator==(WeightableDistribution const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool WeightableDistribution::operator<(WeightableDistribution const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

bool WeightableDistribution::AreEquivalent(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<earthmodel::EarthModel const> second_earth_model, std::shared_ptr<CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

//---------------
// class NormalizationConstant : WeightableDistribution, PhysicallyNormalizedDistribution
//---------------
//
NormalizationConstant::NormalizationConstant() {}

NormalizationConstant::NormalizationConstant(double norm) {
    SetNormalization(norm);
}

double NormalizationConstant::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return 1.0;
}

std::string NormalizationConstant::Name() const {
    return "NormalizationConstant";
}

bool NormalizationConstant::equal(WeightableDistribution const & distribution) const {
    return false;
}

bool NormalizationConstant::less(WeightableDistribution const & distribution) const {
    return false;
}

//---------------
// class InjectionDistribution
//---------------

void InjectionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
}


//---------------
// class PrimaryInjector : InjectionDistribution
//---------------

PrimaryInjector::PrimaryInjector(LeptonInjector::Particle::ParticleType primary_type, double primary_mass) :
    primary_type(primary_type),
    primary_mass(primary_mass)
{}

LeptonInjector::Particle::ParticleType PrimaryInjector::PrimaryType() const {
    return primary_type;
}

double PrimaryInjector::PrimaryMass() const {
    return primary_mass;
}

void PrimaryInjector::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    record.signature.primary_type = primary_type;
    record.primary_mass = primary_mass;
}
double PrimaryInjector::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    if(record.signature.primary_type != primary_type)
        return 0.0;
    if(2.0 * abs(record.primary_mass - primary_mass) / (record.primary_mass + primary_mass) > 1e-9) {
        std::cerr << "Event primary mass does not match injector primary mass!" << std::endl;
        std::cerr << "Event primary_mass: " << record.primary_mass << std::endl;
        std::cerr << "Injector primary_mass: " << primary_mass << std::endl;
        std::cerr << "Particle mass definitions should be consistent." << std::endl;
        std::cerr << "Are you using the wrong simulation?" << std::endl;
        return 0.0;
    }
    return 1.0;
}

std::vector<std::string> PrimaryInjector::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryInjector::Name() const {
    return "PrimaryInjector";
}

std::shared_ptr<InjectionDistribution> PrimaryInjector::clone() const {
    return std::shared_ptr<InjectionDistribution>(new PrimaryInjector(*this));
}

bool PrimaryInjector::equal(WeightableDistribution const & other) const {
    const PrimaryInjector* x = dynamic_cast<const PrimaryInjector*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(primary_type, primary_mass)
            ==
            std::tie(x->primary_type, x->primary_mass);
}

bool PrimaryInjector::less(WeightableDistribution const & other) const {
    const PrimaryInjector* x = dynamic_cast<const PrimaryInjector*>(&other);
    return
        std::tie(primary_type, primary_mass)
        <
        std::tie(x->primary_type, x->primary_mass);
}


//---------------
// class TargetMomentumDistribution : InjectionDistribution
//---------------

void TargetMomentumDistribution::Sample(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel const> earth_model,
        std::shared_ptr<CrossSectionCollection const> cross_sections,
        InteractionRecord & record) const {
    record.target_momentum = SampleMomentum(rand, earth_model, cross_sections, record);
}

std::vector<std::string> TargetMomentumDistribution::DensityVariables() const {
    return std::vector<std::string>{"TargetMomentum"};
}

//---------------
// class TargetAtRest : TargetMomentumDistribution : InjectionDistribution
//---------------
std::array<double, 4> TargetAtRest::SampleMomentum(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel const> earth_model,
        std::shared_ptr<CrossSectionCollection const> cross_sections,
        InteractionRecord const & record) const {
    return std::array<double, 4>{record.target_mass, 0, 0, 0};
}

double TargetAtRest::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return 1.0;
}

std::vector<std::string> TargetAtRest::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> TargetAtRest::clone() const {
    return std::shared_ptr<InjectionDistribution>(new TargetAtRest(*this));
}

std::string TargetAtRest::Name() const {
    return "TargetAtRest";
}

bool TargetAtRest::equal(WeightableDistribution const & other) const {
    const TargetAtRest* x = dynamic_cast<const TargetAtRest*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool TargetAtRest::less(WeightableDistribution const & other) const {
    return false;
}

//---------------
// class PrimaryEnergyDistribution : InjectionDistribution
//---------------
void PrimaryEnergyDistribution::Sample(
        std::shared_ptr<LI_random> rand,
        std::shared_ptr<earthmodel::EarthModel const> earth_model,
        std::shared_ptr<CrossSectionCollection const> cross_sections,
        InteractionRecord & record) const {
    record.primary_momentum[0] = SampleEnergy(rand, earth_model, cross_sections, record);
}

std::vector<std::string> PrimaryEnergyDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryEnergy"};
}

//---------------
// class PowerLaw : PrimaryEnergyDistribution
//---------------
PowerLaw::PowerLaw(double powerLawIndex, double energyMin, double energyMax)
    : powerLawIndex(powerLawIndex)
    , energyMin(energyMin)
    , energyMax(energyMax)
{}

double PowerLaw::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    if(energyMin == energyMax)
        return energyMin; //return the only allowed energy

    if(powerLawIndex == 1.0) //sample uniformly in log space
        return pow(10.0, rand->Uniform(log10(energyMin), log10(energyMax)));
    else {
        double u = rand->Uniform();
        double energyP = (1 - u) * pow(energyMin, 1 - powerLawIndex) + u * pow(energyMax, 1 - powerLawIndex);
        return pow(energyP, 1 / (1 - powerLawIndex));
    }
}

double PowerLaw::pdf(double energy) const {
    if(energyMin == energyMax)
        return 1.0; // only one allowed energy

    if(powerLawIndex == 1.0)
        return 1.0 / (energy * log(energyMax / energyMin));
    else {
        return pow(energy, -powerLawIndex) * (powerLawIndex - 1.0) * (pow(energyMin, powerLawIndex - 1.0) - pow(energyMax, powerLawIndex - 1.0));
    }
}

double PowerLaw::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return pdf(record.primary_momentum[0]);
}

std::string PowerLaw::Name() const {
    return "PowerLaw";
}

std::shared_ptr<InjectionDistribution> PowerLaw::clone() const {
    return std::shared_ptr<InjectionDistribution>(new PowerLaw(*this));
}

bool PowerLaw::equal(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, powerLawIndex)
            ==
            std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

bool PowerLaw::less(WeightableDistribution const & other) const {
    const PowerLaw* x = dynamic_cast<const PowerLaw*>(&other);
    return
        std::tie(energyMin, energyMax, powerLawIndex)
        <
        std::tie(x->energyMin, x->energyMax, x->powerLawIndex);
}

void PowerLaw::SetNormalizationAtEnergy(double norm, double energy) {
    SetNormalization(norm / pdf(energy));
}

//---------------
// class ModifiedMoyalPlusExponentialEnergyDistribution : PrimaryEnergyDistribution
//---------------

double ModifiedMoyalPlusExponentialEnergyDistribution::unnormed_pdf(double energy) const {
    double x = (energy - mu) / sigma;
    double moyal = (A / sigma) * std::exp(-(x + std::exp(-x))/2) / std::sqrt(2.0 * M_PI);
    double exponential = (B / l) * std::exp(-energy / l);
    return moyal + exponential;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

ModifiedMoyalPlusExponentialEnergyDistribution::ModifiedMoyalPlusExponentialEnergyDistribution(double energyMin, double energyMax, double mu, double sigma, double A, double l, double B, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , mu(mu)
    , sigma(sigma)
    , A(A)
    , l(l)
    , B(B)
{
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = earthmodel::Integration::rombergIntegrate(integrand, energyMin, energyMax);
    if(has_physical_normalization) {
        SetNormalization(integral);
    }
}

double ModifiedMoyalPlusExponentialEnergyDistribution::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    // Metropolis-Hastings algorithm to sample from PDF.
    // Pass in a function pointer for the PDF

    double energy, density, test_energy, test_density, odds;
    bool accept;

    // sample an initial point uniformly
    energy = rand->Uniform(energyMin, energyMax);
    density = pdf(energy);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        test_energy = rand->Uniform(energyMin, energyMax);
        test_density = pdf(test_energy);
        odds = test_density / density;
        accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
        if(accept) {
            energy = test_energy;
            density = test_density;
        }
    }

    return energy;
}

double ModifiedMoyalPlusExponentialEnergyDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string ModifiedMoyalPlusExponentialEnergyDistribution::Name() const {
    return "ModifiedMoyalPlusExponentialEnergyDistribution";
}

std::shared_ptr<InjectionDistribution> ModifiedMoyalPlusExponentialEnergyDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ModifiedMoyalPlusExponentialEnergyDistribution(*this));
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::equal(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, mu, sigma, A, l, B)
            ==
            std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

bool ModifiedMoyalPlusExponentialEnergyDistribution::less(WeightableDistribution const & other) const {
    const ModifiedMoyalPlusExponentialEnergyDistribution* x = dynamic_cast<const ModifiedMoyalPlusExponentialEnergyDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, mu, sigma, A, l, B)
        <
        std::tie(x->energyMin, x->energyMax, x->mu, x->sigma, x->A, x->l, x->B);
}

//---------------
// class TabulatedFluxDistribution : PrimaryEnergyDistribution
//---------------
TabulatedFluxDistribution::TabulatedFluxDistribution() {}

void TabulatedFluxDistribution::ComputeIntegral() {
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    integral = earthmodel::Integration::rombergIntegrate(integrand, energyMin, energyMax);
}

void TabulatedFluxDistribution::LoadFluxTable() {
    if(fexists(fluxTableFilename)) {
        std::ifstream in(fluxTableFilename.c_str());
        std::string buf;
        std::string::size_type pos;
        TableData1D<double> table_data;
        while(std::getline(in, buf)) {
            // Ignore comments and blank lines
            if((pos = buf.find('#')) != std::string::npos)
                buf.erase(pos);
            const char* whitespace=" \n\r\t\v";
            if((pos=buf.find_first_not_of(whitespace))!=0)
                buf.erase(0,pos);
            if(!buf.empty() && (pos=buf.find_last_not_of(whitespace))!=buf.size()-1)
                buf.erase(pos+1);
            if(buf.empty())
                continue;

            std::stringstream ss(buf);
            double x, f;
            ss >> x >> f;
            table_data.x.push_back(x);
            table_data.f.push_back(f);
        }
        // If no physical are manually set, use first/last entry of table
        if(not bounds_set) {
            energyMin = table_data.x[0];
            energyMax = table_data.x[table_data.x.size()-1];
        }
        fluxTable = Interpolator1D<double>(table_data);
    } else {
        throw std::runtime_error("Failed to open flux table file!");
    }
}

double TabulatedFluxDistribution::unnormed_pdf(double energy) const {
    return fluxTable(energy);
}

double TabulatedFluxDistribution::pdf(double energy) const {
    return unnormed_pdf(energy) / integral;
}

void TabulatedFluxDistribution::SetEnergyBounds(double eMin, double eMax) {
    energyMin = eMin;
    energyMax = eMax;
    bounds_set = true;
    ComputeIntegral();
}

TabulatedFluxDistribution::TabulatedFluxDistribution(std::string fluxTableFilename, bool has_physical_normalization)
    : bounds_set(false)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

TabulatedFluxDistribution::TabulatedFluxDistribution(double energyMin, double energyMax, std::string fluxTableFilename, bool has_physical_normalization)
    : energyMin(energyMin)
    , energyMax(energyMax)
    , bounds_set(true)
    , fluxTableFilename(fluxTableFilename)
{
    LoadFluxTable();
    std::function<double(double)> integrand = [&] (double x) -> double {
        return unnormed_pdf(x);
    };
    ComputeIntegral();
    if(has_physical_normalization)
        SetNormalization(integral);
}

double TabulatedFluxDistribution::SampleEnergy(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    // Metropolis-Hastings algorithm to sample from PDF.
    // Pass in a function pointer for the PDF

    double energy, density, test_energy, test_density, odds;
    bool accept;

    // sample an initial point uniformly
    energy = rand->Uniform(energyMin, energyMax);
    density = pdf(energy);

    // Metropolis Hastings loop
    for (size_t j = 0; j <= burnin; ++j) {
        test_energy = rand->Uniform(energyMin, energyMax);
        test_density = pdf(test_energy);
        odds = test_density / density;
        accept = (odds > 1.) or (rand->Uniform(0,1) < odds);
        if(accept) {
            energy = test_energy;
            density = test_density;
        }
    }

    return energy;
}

double TabulatedFluxDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    double const & energy = record.primary_momentum[0];
    if(energy < energyMin or energy > energyMax)
        return 0.0;
    else
        return pdf(energy);
}

std::string TabulatedFluxDistribution::Name() const {
    return "TabulatedFluxDistribution";
}

std::shared_ptr<InjectionDistribution> TabulatedFluxDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new TabulatedFluxDistribution(*this));
}

bool TabulatedFluxDistribution::equal(WeightableDistribution const & other) const {
    const TabulatedFluxDistribution* x = dynamic_cast<const TabulatedFluxDistribution*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(energyMin, energyMax, fluxTableFilename)
            ==
            std::tie(x->energyMin, x->energyMax, x->fluxTableFilename);
}

bool TabulatedFluxDistribution::less(WeightableDistribution const & other) const {
    const TabulatedFluxDistribution* x = dynamic_cast<const TabulatedFluxDistribution*>(&other);
    return
        std::tie(energyMin, energyMax, fluxTable)
        <
        std::tie(x->energyMin, x->energyMax, x->fluxTable);
}

//---------------
// class PrimaryDirectionDistribution : InjectionDistribution
//---------------
void PrimaryDirectionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir = SampleDirection(rand, earth_model, cross_sections, record);
    double energy = record.primary_momentum[0];
    double mass = record.primary_mass;
    double momentum = std::sqrt(energy*energy - mass*mass);
    record.primary_momentum[1] = momentum * dir.GetX();
    record.primary_momentum[2] = momentum * dir.GetY();
    record.primary_momentum[3] = momentum * dir.GetZ();
}

std::vector<std::string> PrimaryDirectionDistribution::DensityVariables() const {
    return std::vector<std::string>{"PrimaryDirection"};
}

//---------------
// class IsotropicDirection : PrimaryDirectionDistribution
//---------------
earthmodel::Vector3D IsotropicDirection::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    double nx = rand->Uniform(0, 1);
    double ny = rand->Uniform(0, 1);
    double nz = rand->Uniform(0, 1);
    earthmodel::Vector3D res(nx, ny, nz);
    res.normalize();
    return res;
}

double IsotropicDirection::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return 1.0 / (4.0 * M_PI);
}

std::shared_ptr<InjectionDistribution> IsotropicDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new IsotropicDirection(*this));
}

std::string IsotropicDirection::Name() const {
    return "IsotropicDirection";
}

bool IsotropicDirection::equal(WeightableDistribution const & other) const {
    const IsotropicDirection* x = dynamic_cast<const IsotropicDirection*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool IsotropicDirection::less(WeightableDistribution const & other) const {
    return false;
}

//---------------
// class FixedDirection : PrimaryDirectionDistribution
//---------------
earthmodel::Vector3D FixedDirection::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return dir;
}

double FixedDirection::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    if(abs(1.0 - earthmodel::scalar_product(dir, event_dir)) < 1e-9)
        return 1.0;
    else
        return 0.0;
}

std::vector<std::string> FixedDirection::DensityVariables() const {
    return std::vector<std::string>();
}

std::shared_ptr<InjectionDistribution> FixedDirection::clone() const {
    return std::shared_ptr<InjectionDistribution>(new FixedDirection(*this));
}

std::string FixedDirection::Name() const {
    return "FixedDirection";
}

bool FixedDirection::equal(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9);
}

bool FixedDirection::less(WeightableDistribution const & other) const {
    const FixedDirection* x = dynamic_cast<const FixedDirection*>(&other);
    if(abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        double X = dir.GetX();
        double Y = dir.GetY();
        double Z = dir.GetZ();
        double other_X = dir.GetX();
        double other_Y = dir.GetY();
        double other_Z = dir.GetZ();
        return
            std::tie(X, Y, Z)
            <
            std::tie(other_X, other_Y, other_Z);
    }
}

//---------------
// class Cone : PrimaryDirectionDistribution
//---------------
Cone::Cone(earthmodel::Vector3D dir, double min_angle, double max_angle, double min_azimuth, double max_azimuth) : dir(dir), min_angle(min_angle), max_angle(max_angle), min_azimuth(min_azimuth), max_azimuth(max_azimuth) {
    this->dir.normalize();
    if(this->dir == earthmodel::Vector3D(0,0,1)) {
        rotation = earthmodel::Quaternion(0,0,0,1);
    } else {
        earthmodel::Vector3D r = cross_product(earthmodel::Vector3D(0, 0, 1), dir);
        r.normalize();
        rotation = earthmodel::Quaternion(r);
        rotation.SetW(1.0 + dir.GetZ());
    }
}

earthmodel::Vector3D Cone::SampleDirection(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const{
    double theta = cos(rand->Uniform(acos(max_angle), acos(min_angle)));
    double phi = rand->Uniform(min_azimuth, max_azimuth);
    earthmodel::Quaternion q;
    q.SetEulerAnglesZXZr(phi, theta, 0.0);
    return rotation.rotate(q.rotate(earthmodel::Vector3D(0,0,1), false), false);
}

double Cone::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D event_dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    event_dir.normalize();
    double theta = acos(earthmodel::scalar_product(dir, event_dir));
    earthmodel::Vector3D cone_frame_dir = rotation.rotate(event_dir, true);
    double phi = std::atan2(cone_frame_dir.GetY(). cone_frame_dir.GetX());
    if(theta >= min_angle and theta <= max_angle and phi >= min_azimuth and phi <= max_azimuth)
        return 1.0 / ((max_azimuth - min_azimuth) * (acos(min_angle) - acos(max_angle)));
    else
        return 0.0;
}

std::shared_ptr<InjectionDistribution> Cone::clone() const {
    return std::shared_ptr<InjectionDistribution>(new Cone(*this));
}

std::string Cone::Name() const {
    return "Cone";
}

bool Cone::equal(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);

    if(!x)
        return false;
    else
        return (abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9
            and min_angle == x->min_angle)
            and max_angle == x-> max_angle;
}

bool Cone::less(WeightableDistribution const & other) const {
    const Cone* x = dynamic_cast<const Cone*>(&other);
    if(abs(1.0 - earthmodel::scalar_product(dir, x->dir)) < 1e-9) {
        return false;
    } else {
        return std::tie(dir, min_angle, max_angle)
            <
            std::tie(x->dir, x->min_angle, x->max_angle);
    }
}

//---------------
// class VertexPositionDistribution : InjectionDistribution
//---------------
void VertexPositionDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D pos = SamplePosition(rand, earth_model, cross_sections, record);
    record.interaction_vertex[0] = pos.GetX();
    record.interaction_vertex[1] = pos.GetY();
    record.interaction_vertex[2] = pos.GetZ();
}

std::vector<std::string> VertexPositionDistribution::DensityVariables() const {
    return std::vector<std::string>{"InteractionVertexPosition"};
}

bool VertexPositionDistribution::AreEquivalent(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<earthmodel::EarthModel const> second_earth_model, std::shared_ptr<CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution) and earth_model->operator==(*second_earth_model) and cross_sections->operator==(*second_cross_sections);
}

//---------------
// class OrientedCylinderPositionDistribution : VertexPositionDistribution
//---------------
//
earthmodel::Vector3D OrientedCylinderPositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D OrientedCylinderPositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    /*
    std::pair<earthmodel::Vector3D, earthmodel::Vector3D> GetBounds(earth_model, cross_sections, pca);

    earthmodel::Vector3D p0;
    earthmodel::Vector3D p1;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    */
    return earthmodel::Vector3D();
}

double OrientedCylinderPositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    return 0.0;
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> OrientedCylinderPositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & interaction) const {
    return std::make_pair(earthmodel::Vector3D(), earthmodel::Vector3D());
}

bool OrientedCylinderPositionDistribution::AreEquivalent(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<earthmodel::EarthModel const> second_earth_model, std::shared_ptr<CrossSectionCollection const> second_cross_sections) const {
    return false;
}

//---------------
// class CylinderVolumePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D CylinderVolumePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    double t = rand->Uniform(0, 2 * M_PI);
    const double outer_radius = cylinder.GetRadius();
    const double inner_radius = cylinder.GetInnerRadius();
    const double height = cylinder.GetZ();
    double r = std::sqrt(rand->Uniform(inner_radius*inner_radius, outer_radius*outer_radius));
    double z = rand->Uniform(-height/2.0, height/2.0);
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), z);
    return cylinder.LocalToGlobalPosition(pos);
}

double CylinderVolumePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D pos(record.interaction_vertex);
    double z = pos.GetZ();
    double r = sqrt(pos.GetX() * pos.GetX() + pos.GetY() * pos.GetY());
    if(abs(z) >= 0.5 * cylinder.GetZ()
            or r <= cylinder.GetInnerRadius()
            or r >= cylinder.GetRadius()) {
        return 0.0;
    } else {
        return 1.0 / ((cylinder.GetRadius() * cylinder.GetRadius() - cylinder.GetInnerRadius() * cylinder.GetInnerRadius()) * cylinder.GetZ());
    }
}


CylinderVolumePositionDistribution::CylinderVolumePositionDistribution(earthmodel::Cylinder cylinder) : cylinder(cylinder) {}

std::string CylinderVolumePositionDistribution::Name() const {
    return "CylinderVolumePositionDistribution";
}

std::shared_ptr<InjectionDistribution> CylinderVolumePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new CylinderVolumePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> CylinderVolumePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & interaction) const {
    earthmodel::Vector3D dir(interaction.primary_momentum[1], interaction.primary_momentum[2], interaction.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pos(interaction.interaction_vertex);
    std::vector<earthmodel::Geometry::Intersection> intersections = cylinder.Intersections(pos, dir);
    earthmodel::EarthModel::SortIntersections(intersections);
    if(intersections.size() == 0) {
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));
    } else if(intersections.size() >= 2) {
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(intersections.front().position, intersections.back().position);
    } else {
        throw std::runtime_error("Only found one cylinder intersection!");
    }
}

bool CylinderVolumePositionDistribution::equal(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (cylinder == x->cylinder);
}

bool CylinderVolumePositionDistribution::less(WeightableDistribution const & other) const {
    const CylinderVolumePositionDistribution* x = dynamic_cast<const CylinderVolumePositionDistribution*>(&other);
    return cylinder < x->cylinder;
}

bool CylinderVolumePositionDistribution::AreEquivalent(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<earthmodel::EarthModel const> second_earth_model, std::shared_ptr<CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

//---------------
// class DepthFunction
//---------------
DepthFunction::DepthFunction() {}

double DepthFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool DepthFunction::operator==(DepthFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool DepthFunction::operator<(DepthFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class LeptonDepthFunction : DepthFunction
//---------------
LeptonDepthFunction::LeptonDepthFunction() {}

void LeptonDepthFunction::SetMuParams(double mu_alpha, double mu_beta) {
    this->mu_alpha = mu_alpha;
    this->mu_beta = mu_beta;
}

void LeptonDepthFunction::SetTauParams(double tau_alpha, double tau_beta) {
    this->tau_alpha = tau_alpha;
    this->tau_beta = tau_beta;

}

void LeptonDepthFunction::SetScale(double scale) {
    this->scale = scale;
}

void LeptonDepthFunction::SetMaxDepth(double max_depth) {
    this->max_depth = max_depth;
}

double LeptonDepthFunction::GetMuAlpha() const {
    return mu_alpha;
}

double LeptonDepthFunction::GetMuBeta() const {
    return mu_beta;
}

double LeptonDepthFunction::GetTauAlpha() const {
    return tau_alpha;
}

double LeptonDepthFunction::GetTauBeta() const {
    return tau_beta;
}

double LeptonDepthFunction::GetScale() const {
    return scale;
}

double LeptonDepthFunction::GetMaxDepth() const {
    return max_depth;
}

double LeptonDepthFunction::operator()(InteractionSignature const & signature, double energy) const {
    double range = log(1.0 + energy * mu_beta / mu_alpha) / mu_beta;
    if(tau_primaries.count(signature.primary_type) > 0)
        range += log(1.0 + energy * tau_beta / tau_alpha) / tau_beta;
    return std::min(range, max_depth);
}

bool LeptonDepthFunction::equal(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        ==
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

bool LeptonDepthFunction::less(DepthFunction const & other) const {
    const LeptonDepthFunction* x = dynamic_cast<const LeptonDepthFunction*>(&other);

    if(not x)
        return false;

    return
        std::tie(mu_alpha, mu_beta, tau_alpha, tau_beta, scale, max_depth, tau_primaries)
        <
        std::tie(x->mu_alpha, x->mu_beta, x->tau_alpha, x->tau_beta, x->scale, x->max_depth, x->tau_primaries);
}

//---------------
// class RangeFunction
//---------------
RangeFunction::RangeFunction() {}

double RangeFunction::operator()(InteractionSignature const & signature, double energy) const {
    return 0.0;
}

bool RangeFunction::operator==(RangeFunction const & distribution) const {
    if(this == &distribution)
        return true;
    else
        return this->equal(distribution);
}

bool RangeFunction::operator<(RangeFunction const & distribution) const {
    if(typeid(this) == typeid(&distribution))
        return this->less(distribution);
    else
        return std::type_index(typeid(this)) < std::type_index(typeid(&distribution));
}

//---------------
// class DecayRangeFunction
//---------------
//
//
double DecayRangeFunction::DecayLength(double particle_mass, double decay_width, double energy) {
    double beta = sqrt(energy*energy - particle_mass*particle_mass) / energy;
    double gamma = energy / particle_mass;
    double time_in_rest_frame = 1.0 / decay_width; // inverse GeV
    double time_in_lab_frame = time_in_rest_frame * gamma; // inverse GeV
    constexpr double iGeV_in_m = 1.973269804593025e-16; // meters per inverse GeV
    double length = time_in_lab_frame * beta * iGeV_in_m; // meters = ((inverse GeV * dimensionless) * (meters per inverse GeV))
    return length; // meters
}

double DecayRangeFunction::DecayLength(InteractionSignature const & signature, double energy) const {
    return DecayRangeFunction::DecayLength(particle_mass, decay_width, energy);
}

double DecayRangeFunction::Range(InteractionSignature const & signature, double energy) const {
    return std::min(DecayLength(signature, energy) * multiplier, max_distance);
}

double DecayRangeFunction::operator()(InteractionSignature const & signature, double energy) const {
    return Range(signature, energy);
}

double DecayRangeFunction::Multiplier() const {
    return multiplier;
}

double DecayRangeFunction::ParticleMass() const {
    return particle_mass;
}

double DecayRangeFunction::DecayWidth() const {
    return decay_width;
}

double DecayRangeFunction::MaxDistance() const {
    return max_distance;
}

DecayRangeFunction::DecayRangeFunction(double particle_mass, double decay_width, double multiplier, double max_distance) : particle_mass(particle_mass), decay_width(decay_width), multiplier(multiplier), max_distance(max_distance) {}

bool DecayRangeFunction::equal(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(particle_mass, decay_width, multiplier, max_distance)
            ==
            std::tie(x->particle_mass, x->decay_width, x->multiplier, x->max_distance);
}

bool DecayRangeFunction::less(RangeFunction const & other) const {
    const DecayRangeFunction* x = dynamic_cast<const DecayRangeFunction*>(&other);

    return
        std::tie(particle_mass, decay_width, multiplier, max_distance)
        <
        std::tie(x->particle_mass, x->decay_width, x->multiplier, x->max_distance);
}

//---------------
// class ColumnDepthPositionDistribution : VertexPositionDistribution
//---------------
earthmodel::Vector3D ColumnDepthPositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D ColumnDepthPositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);
    if(total_interaction_depth == 0) {
        throw(InjectionFailure("No available interactions along path!"));
    }

    double traversed_interaction_depth;
    if(total_interaction_depth < 1e-6) {
        traversed_interaction_depth = rand->Uniform() * total_interaction_depth;
    } else {
        double exp_m_total_interaction_depth = exp(-total_interaction_depth);

        double y = rand->Uniform();
        traversed_interaction_depth = -log(y * exp_m_total_interaction_depth + (1.0 - y));
    }

    double dist = path.GetDistanceFromStartAlongPath(traversed_interaction_depth, targets, total_cross_sections);
    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double ColumnDepthPositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    double interaction_density = earth_model->GetInteractionDensity(path.GetIntersections(), earth_model->GetEarthCoordPosFromDetCoordPos(vertex), targets, total_cross_sections);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

ColumnDepthPositionDistribution::ColumnDepthPositionDistribution(double radius, double endcap_length, std::shared_ptr<DepthFunction> depth_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), depth_function(depth_function), target_types(target_types) {}

std::string ColumnDepthPositionDistribution::Name() const {
    return "ColumnDepthPositionDistribution";
}

std::shared_ptr<InjectionDistribution> ColumnDepthPositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new ColumnDepthPositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> ColumnDepthPositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double lepton_depth = (*depth_function)(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByColumnDepth(lepton_depth);
    path.ClipToOuterBounds();
    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool ColumnDepthPositionDistribution::equal(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (depth_function and x->depth_function and *depth_function == *x->depth_function)
                    or (!depth_function and !x->depth_function)
                )
            and target_types == x->target_types);
}

bool ColumnDepthPositionDistribution::less(WeightableDistribution const & other) const {
    const ColumnDepthPositionDistribution* x = dynamic_cast<const ColumnDepthPositionDistribution*>(&other);
    bool depth_less =
        (!depth_function and x->depth_function) // this->NULL and other->(not NULL)
        or (depth_function and x->depth_function // both not NULL
                and *depth_function < *x->depth_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, depth_less, x->target_types);
}

//---------------
// class RangePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D RangePositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D RangePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);
    if(total_interaction_depth == 0) {
        throw(InjectionFailure("No available interactions along path!"));
    }
    double traversed_interaction_depth;
    if(total_interaction_depth < 1e-6) {
        traversed_interaction_depth = rand->Uniform() * total_interaction_depth;
    } else {
        double exp_m_total_interaction_depth = exp(-total_interaction_depth);

        double y = rand->Uniform();
        traversed_interaction_depth = -log(y * exp_m_total_interaction_depth + (1.0 - y));
    }

    double dist = path.GetDistanceFromStartAlongPath(traversed_interaction_depth, targets, total_cross_sections);
    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double RangePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)))
        return 0.0;

    std::set<Particle::ParticleType> const & possible_targets = cross_sections->TargetTypes();

    std::vector<LeptonInjector::Particle::ParticleType> targets(possible_targets.begin(), possible_targets.end());
    std::vector<double> total_cross_sections(targets.size(), 0.0);
    InteractionRecord fake_record = record;
    for(unsigned int i=0; i<targets.size(); ++i) {
        LeptonInjector::Particle::ParticleType const & target = targets[i];
        fake_record.signature.target_type = target;
        fake_record.target_mass = earth_model->GetTargetMass(target);
        fake_record.target_momentum = {fake_record.target_mass,0,0,0};
        for(auto const & cross_section : cross_sections->GetCrossSectionsForTarget(target)) {
            total_cross_sections[i] += cross_section->TotalCrossSection(fake_record);
        }
    }
    double total_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    path.SetPointsWithRay(path.GetFirstPoint(), path.GetDirection(), path.GetDistanceFromStartInBounds(earth_model->GetEarthCoordPosFromDetCoordPos(vertex)));

    double traversed_interaction_depth = path.GetInteractionDepthInBounds(targets, total_cross_sections);

    double interaction_density = earth_model->GetInteractionDensity(path.GetIntersections(), earth_model->GetEarthCoordPosFromDetCoordPos(vertex), targets, total_cross_sections);

    double prob_density;
    if(total_interaction_depth < 1e-6) {
        prob_density = interaction_density / total_interaction_depth;
    } else {
        prob_density = interaction_density * exp(-log_one_minus_exp_of_negative(total_interaction_depth) - traversed_interaction_depth);
    }
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3

    return prob_density;
}

RangePositionDistribution::RangePositionDistribution() {}

RangePositionDistribution::RangePositionDistribution(double radius, double endcap_length, std::shared_ptr<RangeFunction> range_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {}

std::string RangePositionDistribution::Name() const {
    return "RangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> RangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new RangePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> RangePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double lepton_range = range_function->operator()(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(lepton_range);
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));
    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool RangePositionDistribution::equal(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            and target_types == x->target_types);
}

bool RangePositionDistribution::less(WeightableDistribution const & other) const {
    const RangePositionDistribution* x = dynamic_cast<const RangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, range_less, x->target_types);
}

//---------------
// class DecayRangePositionDistribution : public VertexPositionDistribution {
//---------------
earthmodel::Vector3D DecayRangePositionDistribution::SampleFromDisk(std::shared_ptr<LI_random> rand, earthmodel::Vector3D const & dir) const {
    double t = rand->Uniform(0, 2 * M_PI);
    double r = radius * std::sqrt(rand->Uniform());
    earthmodel::Vector3D pos(r * cos(t), r * sin(t), 0.0);
    earthmodel::Quaternion q = rotation_between(earthmodel::Vector3D(0,0,1), dir);
    return q.rotate(pos, false);
}

earthmodel::Vector3D DecayRangePositionDistribution::SamplePosition(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D pca = SampleFromDisk(rand, dir);

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    double y = rand->Uniform();
    double total_distance = path.GetDistance();
    double dist = -decay_length * log(y * (exp(-total_distance/decay_length) - 1) + 1);

    earthmodel::Vector3D vertex = earth_model->GetDetCoordPosFromEarthCoordPos(path.GetFirstPoint() + dist * path.GetDirection());

    return vertex;
}

double DecayRangePositionDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return 0.0;

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return 0.0;

    double total_distance = path.GetDistance();
    double dist = earthmodel::scalar_product(path.GetDirection(), vertex - path.GetFirstPoint());

    double prob_density = exp(-dist / decay_length) / (decay_length * (1.0 - exp(-total_distance / decay_length))); // m^-1
    prob_density /= (M_PI * radius * radius); // (m^-1 * m^-2) -> m^-3
    return prob_density;
}

DecayRangePositionDistribution::DecayRangePositionDistribution() {}

DecayRangePositionDistribution::DecayRangePositionDistribution(double radius, double endcap_length, std::shared_ptr<DecayRangeFunction> range_function, std::set<Particle::ParticleType> target_types) : radius(radius), endcap_length(endcap_length), range_function(range_function), target_types(target_types) {}

std::string DecayRangePositionDistribution::Name() const {
    return "DecayRangePositionDistribution";
}

std::shared_ptr<InjectionDistribution> DecayRangePositionDistribution::clone() const {
    return std::shared_ptr<InjectionDistribution>(new DecayRangePositionDistribution(*this));
}

std::pair<earthmodel::Vector3D, earthmodel::Vector3D> DecayRangePositionDistribution::InjectionBounds(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    earthmodel::Vector3D dir(record.primary_momentum[1], record.primary_momentum[2], record.primary_momentum[3]);
    dir.normalize();
    earthmodel::Vector3D vertex(record.interaction_vertex); // m
    earthmodel::Vector3D pca = vertex - dir * earthmodel::scalar_product(dir, vertex);

    if(pca.magnitude() >= radius)
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    double decay_length = range_function->DecayLength(record.signature, record.primary_momentum[0]);

    earthmodel::Vector3D endcap_0 = pca - endcap_length * dir;
    earthmodel::Vector3D endcap_1 = pca + endcap_length * dir;

    earthmodel::Path path(earth_model, earth_model->GetEarthCoordPosFromDetCoordPos(endcap_0), earth_model->GetEarthCoordDirFromDetCoordDir(dir), endcap_length*2);
    path.ExtendFromStartByDistance(decay_length * range_function->Multiplier());
    path.ClipToOuterBounds();

    if(not path.IsWithinBounds(vertex))
        return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(earthmodel::Vector3D(0, 0, 0), earthmodel::Vector3D(0, 0, 0));

    return std::pair<earthmodel::Vector3D, earthmodel::Vector3D>(path.GetFirstPoint(), path.GetLastPoint());
}

bool DecayRangePositionDistribution::equal(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);

    if(!x)
        return false;
    else
        return (radius == x->radius
            and endcap_length == x->endcap_length
            and (
                    (range_function and x->range_function and *range_function == *x->range_function)
                    or (!range_function and !x->range_function)
                )
            and target_types == x->target_types);
}

bool DecayRangePositionDistribution::less(WeightableDistribution const & other) const {
    const DecayRangePositionDistribution* x = dynamic_cast<const DecayRangePositionDistribution*>(&other);
    bool range_less =
        (!range_function and x->range_function) // this->NULL and other->(not NULL)
        or (range_function and x->range_function // both not NULL
                and *range_function < *x->range_function); // Less than
    bool f = false;
    return
        std::tie(radius, endcap_length, f, target_types)
        <
        std::tie(radius, x->endcap_length, range_less, x->target_types);
}

bool DecayRangePositionDistribution::AreEquivalent(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, std::shared_ptr<WeightableDistribution const> distribution, std::shared_ptr<earthmodel::EarthModel const> second_earth_model, std::shared_ptr<CrossSectionCollection const> second_cross_sections) const {
    return this->operator==(*distribution);
}

//---------------
// class PrimaryNeutrinoHelicityDistribution : InjectionDistribution
//---------------
void PrimaryNeutrinoHelicityDistribution::Sample(std::shared_ptr<LI_random> rand, std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord & record) const {
    Particle::ParticleType & t = record.signature.primary_type;
    if(t > 0) // Particles are left handed, anti-particles are right handed
        record.primary_helicity = -0.5;
    else
        record.primary_helicity = 0.5;
}

double PrimaryNeutrinoHelicityDistribution::GenerationProbability(std::shared_ptr<earthmodel::EarthModel const> earth_model, std::shared_ptr<CrossSectionCollection const> cross_sections, InteractionRecord const & record) const {
    std::array<double, 4> const & mom = record.primary_momentum;
    earthmodel::Vector3D dir(mom[1], mom[2], mom[3]);
    dir.normalize();

    if(abs(0.5 - abs(record.primary_helicity)) > 1e-9) // Helicity magnitude must be 0.5
        return 0.0;

    Particle::ParticleType const & t = record.signature.primary_type;
    // Particles are left handed, anti-particles are right handed
    if(t > 0) {
        if(record.primary_helicity < 0) // expect opposite direction
            return 1.0;
        else
            return 0.0;
    } else {
        if(record.primary_helicity > 0) // expect same direction
            return 1.0;
        else
            return 0.0;
    }
}

PrimaryNeutrinoHelicityDistribution::PrimaryNeutrinoHelicityDistribution() {}

std::vector<std::string> PrimaryNeutrinoHelicityDistribution::DensityVariables() const {
    return std::vector<std::string>{};
}

std::string PrimaryNeutrinoHelicityDistribution::Name() const {
    return "PrimaryNeutrinoHelicityDistribution";
}

std::shared_ptr<InjectionDistribution> PrimaryNeutrinoHelicityDistribution::clone() const {
    return std::shared_ptr<PrimaryNeutrinoHelicityDistribution>(new PrimaryNeutrinoHelicityDistribution(*this));
}

bool PrimaryNeutrinoHelicityDistribution::equal(WeightableDistribution const & other) const {
    const PrimaryNeutrinoHelicityDistribution* x = dynamic_cast<const PrimaryNeutrinoHelicityDistribution*>(&other);

    if(!x)
        return false;
    else
        return true;
}

bool PrimaryNeutrinoHelicityDistribution::less(WeightableDistribution const & other) const {
    return false;
}

} // namespace LeptonInjector
