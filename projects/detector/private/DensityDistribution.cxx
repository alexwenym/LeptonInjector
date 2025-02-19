#include <cmath>
#include <iostream>
#include "LeptonInjector/detector/DensityDistribution.h"

namespace LI {
namespace detector {

DensityDistribution::DensityDistribution() {}

DensityDistribution::DensityDistribution(const DensityDistribution& density_distr) {}

bool DensityDistribution::operator==(const DensityDistribution& dens_distr) const
{
    if (!this->compare(dens_distr) )
        return false;
    return true;
}


bool DensityDistribution::operator!=(const DensityDistribution& dens_distr) const {
    return !(*this == dens_distr);
}

} // namespace LI
} // namespace detector
