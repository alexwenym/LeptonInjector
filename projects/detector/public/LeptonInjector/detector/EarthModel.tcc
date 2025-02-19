#pragma once
#ifndef LI_EarthModel_TCC
#define LI_EarthModel_TCC

#include <numeric>

#include "LeptonInjector/detector/EarthModel.h"

namespace LI {
namespace detector {

template<typename Iterator, class>
double EarthModel::GetMassDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, Iterator begin, Iterator end) const {
    math::Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();

    std::function<bool(std::vector<geometry::Geometry::Intersection>::const_iterator, std::vector<geometry::Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<geometry::Geometry::Intersection>::const_iterator current_intersection, std::vector<geometry::Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            EarthSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            std::vector<double> mass_fractions = materials_.GetTargetMassFraction(sector.material_id, begin, end);
            density *= std::accumulate(mass_fractions.begin(), mass_fractions.end(), 0.0);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    assert(density >= 0);

    return density;
}

template<typename Iterator, class>
double EarthModel::GetMassDensity(math::Vector3D const & p0, Iterator begin, Iterator end) const {
    math::Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    geometry::Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetMassDensity(intersections, p0, begin, end);
}


template<typename Iterator, class>
std::vector<double> EarthModel::GetParticleDensity(geometry::Geometry::IntersectionList const & intersections, math::Vector3D const & p0, Iterator begin, Iterator end) const {
    math::Vector3D direction = p0 - intersections.position;
    if(direction.magnitude() == 0) {
        direction = intersections.direction;
    } else {
        direction.normalize();
    }
    double dot = direction * intersections.direction;
    assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
    double offset = (intersections.position - p0) * direction;

    if(dot < 0) {
        dot = -1;
    } else {
        dot = 1;
    }
    double density = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> particle_fractions;

    std::function<bool(std::vector<geometry::Geometry::Intersection>::const_iterator, std::vector<geometry::Geometry::Intersection>::const_iterator, double)> callback =
        [&] (std::vector<geometry::Geometry::Intersection>::const_iterator current_intersection, std::vector<geometry::Geometry::Intersection>::const_iterator intersection, double last_point) {
        // The local integration is bounded on the upper end by the intersection
        double end_point = offset + dot * intersection->distance;
        // whereas the lower end is bounded by the end of the last line segment, and the entry into the sector
        double start_point = std::max(offset + dot * current_intersection->distance, offset + dot * last_point);
        if(start_point <= 0 and end_point >= 0) {
            EarthSector sector = GetSector(current_intersection->hierarchy);
            density = sector.density->Evaluate(p0);
            particle_fractions = materials_.GetTargetParticleFraction(sector.material_id, begin, end);
            return true;
        } else {
            return false;
        }
    };

    SectorLoop(callback, intersections, dot < 0);

    for(unsigned int i=0; i<particle_fractions.size(); ++i) {
        particle_fractions[i] *= density;
    }

    assert(density >= 0);

    return particle_fractions;
}

template<typename Iterator, typename>
std::vector<double> EarthModel::GetParticleDensity(math::Vector3D const & p0, Iterator begin, Iterator end) const {
    math::Vector3D direction(1,0,0); // Any direction will work for determining the sector heirarchy
    geometry::Geometry::IntersectionList intersections = GetIntersections(p0, direction);
    return GetParticleDensity(intersections, p0, begin, end);
}


} // namespace detector
} // namespace LI

#endif
