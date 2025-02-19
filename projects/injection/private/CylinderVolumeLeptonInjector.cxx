#include <map>
#include <set>
#include <string>
#include <memory>
#include <vector>

#include "LeptonInjector/math/Vector3D.h"

#include "LeptonInjector/injection/InjectorBase.h"
#include "LeptonInjector/injection/CylinderVolumeLeptonInjector.h"

namespace LI {
namespace injection {

//---------------
// class CylinderVolumeLeptonInjector : InjectorBase
//---------------
CylinderVolumeLeptonInjector::CylinderVolumeLeptonInjector() {}

CylinderVolumeLeptonInjector::CylinderVolumeLeptonInjector(
        unsigned int events_to_inject,
        std::shared_ptr<LI::detector::EarthModel> earth_model,
        std::shared_ptr<injection::InjectionProcess> primary_process,
        std::vector<std::shared_ptr<injection::InjectionProcess>> secondary_processes,
        std::shared_ptr<LI::utilities::LI_random> random,
        LI::geometry::Cylinder cylinder) :
    InjectorBase(events_to_inject, earth_model, random),
    position_distribution(std::make_shared<LI::distributions::CylinderVolumePositionDistribution>(cylinder)) {
    cross_sections = primary_process->cross_sections;
    primary_process->injection_distributions.push_back(position_distribution);
    SetPrimaryProcess(primary_process);
    for(auto & sec_process : secondary_processes) {
      AddSecondaryProcess(sec_process);
      // Assume each secondary process already has a position distribution
      // Otherwise uncomment below
      /*
      sec_process->injection_distributions.push_back(position_distribution);
      */
    }
}

std::string CylinderVolumeLeptonInjector::Name() const {
    return("VolumeInjector");
}

std::pair<LI::math::Vector3D, LI::math::Vector3D> CylinderVolumeLeptonInjector::InjectionBounds(LI::dataclasses::InteractionRecord const & interaction) const {
    return position_distribution->InjectionBounds(earth_model, cross_sections, interaction);
}

} // namespace injection
} // namespace LI
