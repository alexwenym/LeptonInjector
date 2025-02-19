#include "LeptonInjector/dataclasses/InteractionRecord.h"
#include "LeptonInjector/utilities/Random.h"

#include "LeptonInjector/distributions/Distributions.h"
#include "LeptonInjector/distributions/primary/energy/PrimaryEnergyDistribution.h"
#include "LeptonInjector/distributions/primary/energy/Monoenergetic.h"

namespace LI {
namespace distributions {

//---------------
// class Monoenergetic : PrimaryEnergyDistribution
//---------------
Monoenergetic::Monoenergetic(double gen_energy)
    : gen_energy(gen_energy)
{}

double Monoenergetic::SampleEnergy(std::shared_ptr<LI::utilities::LI_random> rand, std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    return gen_energy;
}

double Monoenergetic::pdf(double energy) const {
    if(std::abs(energy - gen_energy) < 1e-6 * gen_energy)
        return 1.0; // only one allowed energy
    return 0.0;
}

double Monoenergetic::GenerationProbability(std::shared_ptr<LI::detector::EarthModel const> earth_model, std::shared_ptr<LI::crosssections::CrossSectionCollection const> cross_sections, LI::dataclasses::InteractionRecord const & record) const {
    return pdf(record.primary_momentum[0]);
}

std::string Monoenergetic::Name() const {
    return "Monoenergetic";
}

std::shared_ptr<InjectionDistribution> Monoenergetic::clone() const {
    return std::shared_ptr<InjectionDistribution>(new Monoenergetic(*this));
}

bool Monoenergetic::equal(WeightableDistribution const & other) const {
    const Monoenergetic* x = dynamic_cast<const Monoenergetic*>(&other);

    if(!x)
        return false;
    else
        return
            std::tie(gen_energy)
            ==
            std::tie(x->gen_energy);
}

bool Monoenergetic::less(WeightableDistribution const & other) const {
    const Monoenergetic* x = dynamic_cast<const Monoenergetic*>(&other);
    return
        std::tie(gen_energy)
        <
        std::tie(x->gen_energy);
}

} // namespace distributions
} // namespace LI

