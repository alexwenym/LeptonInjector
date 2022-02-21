#include <map>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "earthmodel-service/Path.h"
#include "earthmodel-service/Geometry.h"

using namespace earthmodel;

Path::Path() {

}

Path::Path(std::shared_ptr<const EarthModel> earth_model) {
    SetEarthModel(earth_model);
}

Path::Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & last_point) {
    SetEarthModel(earth_model);
    SetPoints(first_point, last_point);
}

Path::Path(std::shared_ptr<const EarthModel> earth_model, Vector3D const & first_point, Vector3D const & direction, double distance) {
    SetEarthModel(earth_model);
    SetPointsWithRay(first_point, direction, distance);
}

bool Path::HasEarthModel() {
    return set_earth_model_;
}

bool Path::HasPoints() {
    return set_points_;
}

bool Path::HasIntersections() {
    return set_intersections_;
}

bool Path::HasTargetColumnDepth(std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    return column_depth_cache_.count(targets) > 0;
}

bool Path::HasColumnDepth() {
    return set_column_depth_;
}

std::shared_ptr<const EarthModel> Path::GetEarthModel() {
    return earth_model_;
}

Vector3D Path::GetFirstPoint() {
    return first_point_;
}

Vector3D Path::GetLastPoint() {
    return last_point_;
}

Vector3D Path::GetDirection() {
    return direction_;
}

double Path::GetDistance() {
    return distance_;
}

Geometry::IntersectionList Path::GetIntersections() {
    return intersections_;
}

void Path::SetEarthModel(std::shared_ptr<const EarthModel> earth_model) {
    earth_model_ = earth_model;
    set_earth_model_ = true;
}

void Path::EnsureEarthModel() {
    if(not set_earth_model_) {
        throw(std::runtime_error("Earth model not set!"));
    }
}

void Path::SetPoints(Vector3D first_point, Vector3D last_point) {
    first_point_ = first_point;
    last_point_ = last_point;
    direction_ = last_point_ - first_point_;
    distance_ = direction_.magnitude();
    direction_.normalize();
    set_points_ = true;
    set_intersections_ = false;
    column_depth_cache_.clear();
    set_column_depth_ = false;
}

void Path::SetPointsWithRay(Vector3D first_point, Vector3D direction, double distance) {
    first_point_ = first_point;
    direction_ = direction;
    direction_.normalize();
    //double dif = std::abs(direction_.magnitude() - direction.magnitude()) / std::max(direction_.magnitude(), direction.magnitude());
    //if(not std::isnan(dif)) assert(dif < 1e-12);
    distance_ = distance;
    last_point_ = first_point + direction * distance;
    set_points_ = true;
    set_intersections_ = false;
    column_depth_cache_.clear();
    set_column_depth_ = false;
}

void Path::EnsurePoints() {
    if(not set_points_) {
        throw(std::runtime_error("Points not set!"));
    }
}

void Path::SetIntersections(Geometry::IntersectionList const & intersections) {
    intersections_ = intersections;
    set_intersections_ = true;
}

void Path::ComputeIntersections() {
    EnsureEarthModel();
    EnsurePoints();
    intersections_ = earth_model_->GetIntersections(first_point_, direction_);
    set_intersections_ = true;
}

void Path::EnsureIntersections() {
    if(not set_intersections_) {
        ComputeIntersections();
    }
}

void Path::ClipToOuterBounds() {
    EnsureIntersections();
    EnsurePoints();
    Geometry::IntersectionList bounds = EarthModel::GetOuterBounds(intersections_);
    if(bounds.intersections.size() > 0) {
        assert(bounds.intersections.size() == 2);
        Vector3D p0 = bounds.intersections[0].position;
        Vector3D p1 = bounds.intersections[1].position;
        Vector3D direction = p1 - p0;
        double distance = direction.magnitude();
        direction.normalize();
        double dot = direction_ * direction;
        assert(std::abs(1.0 - std::abs(dot)) < 1e-6);
        if(dot < 0) {
            p0.swap(p1);
        }
        bool clip_0 = (p0 - first_point_) * direction_ > 0;
        bool clip_1 = (p1 - last_point_) * direction_ < 0;
        bool clip = clip_0 or clip_1;
        if(clip_0) {
            first_point_ = p0;
        }
        if(clip_1) {
            last_point_ = p1;
        }
        if(clip) {
            distance_ = (last_point_ - first_point_).magnitude();
            column_depth_cache_.clear();
            set_column_depth_ = false;
        }
    } else {
        return;
    }
}

void Path::Flip() {
    EnsurePoints();
    std::swap(first_point_, last_point_);
    direction_ *= -1;
}

void Path::ExtendFromEndByDistance(double distance) {
    distance_ += distance;
    last_point_ += direction_ * distance;
    if(distance_ < 0) {
        distance_ = 0;
        last_point_ = first_point_;
    }
    column_depth_cache_.clear();
    set_column_depth_ = false;
}

void Path::ExtendFromEndByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double distance = GetDistanceFromEndAlongPath(column_depth, targets);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromEndByColumnDepth(double column_depth) {
    double distance = GetDistanceFromEndAlongPath(column_depth);
    ExtendFromEndByDistance(distance);
}

void Path::ExtendFromStartByDistance(double distance) {
    distance_ += distance;
    first_point_ += direction_ * -distance;
    if(distance_ < 0) {
        distance_ = 0;
        first_point_ = last_point_;
    }
    column_depth_cache_.clear();
    set_column_depth_ = false;
}

void Path::ExtendFromStartByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double distance = GetDistanceFromStartInReverse(column_depth, targets);
    ExtendFromStartByDistance(distance);
}

void Path::ExtendFromStartByColumnDepth(double column_depth) {
    double distance = GetDistanceFromStartInReverse(column_depth);
    ExtendFromStartByDistance(distance);
}

void Path::ShrinkFromEndByDistance(double distance) {
    ExtendFromEndByDistance(-distance);
}

void Path::ShrinkFromEndByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double distance = GetDistanceFromEndInReverse(column_depth, targets);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromEndByColumnDepth(double column_depth) {
    double distance = GetDistanceFromEndInReverse(column_depth);
    ShrinkFromEndByDistance(distance);
}

void Path::ShrinkFromStartByDistance(double distance) {
    ExtendFromStartByDistance(-distance);
}

void Path::ShrinkFromStartByColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double distance = GetDistanceFromStartAlongPath(column_depth, targets);
    ShrinkFromStartByDistance(distance);
}

void Path::ShrinkFromStartByColumnDepth(double column_depth) {
    double distance = GetDistanceFromStartAlongPath(column_depth);
    ShrinkFromStartByDistance(distance);
}

void Path::ExtendFromEndToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromEndByDistance(shift);
    }
}

void Path::ExtendFromEndToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double shift = column_depth - GetColumnDepthInBounds(targets);
    if(shift > 0) {
        ExtendFromEndByColumnDepth(shift);
    }
}

void Path::ExtendFromEndToColumnDepth(double column_depth) {
    double shift = column_depth - GetColumnDepthInBounds();
    if(shift > 0) {
        ExtendFromEndByColumnDepth(shift);
    }
}

void Path::ExtendFromStartToDistance(double distance) {
    double shift = distance - distance_;
    if(shift > 0) {
        ExtendFromStartByDistance(shift);
    }
}

void Path::ExtendFromStartToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double shift = column_depth - GetColumnDepthInBounds(targets);
    if(shift > 0) {
        ExtendFromStartByColumnDepth(shift);
    }
}

void Path::ExtendFromStartToColumnDepth(double column_depth) {
    double shift = column_depth - GetColumnDepthInBounds();
    if(shift > 0) {
        ExtendFromStartByColumnDepth(shift);
    }
}

void Path::ShrinkFromEndToDistance(double distance) {
    double shift = distance_ - distance;
    if(shift > 0) {
        ShrinkFromEndByDistance(shift);
    }
}

void Path::ShrinkFromEndToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double shift = GetColumnDepthInBounds(targets) - column_depth;
    if(shift > 0) {
        ShrinkFromEndByColumnDepth(shift);
    }
}

void Path::ShrinkFromEndToColumnDepth(double column_depth) {
    double shift = GetColumnDepthInBounds() - column_depth;
    if(shift > 0) {
        ShrinkFromEndByColumnDepth(shift);
    }
}

void Path::ShrinkFromStartToDistance(double distance) {
    double shift = distance_ - distance;
    if(shift > 0) {
        ShrinkFromStartByDistance(shift);
    }
}

void Path::ShrinkFromStartToColumnDepth(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    double shift = GetColumnDepthInBounds(targets) - column_depth;
    if(shift > 0) {
        ShrinkFromStartByColumnDepth(shift);
    }
}

void Path::ShrinkFromStartToColumnDepth(double column_depth) {
    double shift = GetColumnDepthInBounds() - column_depth;
    if(shift > 0) {
        ShrinkFromStartByColumnDepth(shift);
    }
}

double Path::GetColumnDepthInBounds(std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    if(HasTargetColumnDepth(targets)) {
        return column_depth_cache_[targets];
    } else {
        double column_depth = earth_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_, targets);
        column_depth_cache_[targets] = column_depth;
        return column_depth;
    }
}

double Path::GetColumnDepthInBounds() {
    EnsureIntersections();
    if(HasColumnDepth()) {
        return column_depth_cached_;
    } else {
        double column_depth = earth_model_->GetColumnDepthInCGS(intersections_, first_point_, last_point_);
        column_depth_cached_ = column_depth;
        return column_depth;
    }
}

double Path::GetColumnDepthFromStartInBounds(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, targets);
}

double Path::GetColumnDepthFromStartInBounds(double distance) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance);
}

double Path::GetColumnDepthFromEndInBounds(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, targets);
}

double Path::GetColumnDepthFromEndInBounds(double distance) {
    if(distance > distance_) {
        distance = distance_;
    } else if(distance <= 0) {
        return 0.0;
    }
    EnsureIntersections();
    return earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance);
}

double Path::GetColumnDepthFromStartAlongPath(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance, targets), distance);
}

double Path::GetColumnDepthFromStartAlongPath(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * distance), distance);
}

double Path::GetColumnDepthFromEndAlongPath(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * distance, targets), distance);
}

double Path::GetColumnDepthFromEndAlongPath(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * distance), distance);
}

double Path::GetColumnDepthFromStartInReverse(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * -distance, targets), distance);
}

double Path::GetColumnDepthFromStartInReverse(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, first_point_, first_point_ + direction_ * -distance), distance);
}

double Path::GetColumnDepthFromEndInReverse(double distance, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance, targets), distance);
}

double Path::GetColumnDepthFromEndInReverse(double distance) {
    EnsureIntersections();
    return std::copysign(earth_model_->GetColumnDepthInCGS(intersections_, last_point_, last_point_ + direction_ * -distance), distance);
}

double Path::GetDistanceFromStartInBounds(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth, targets);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartInBounds(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth, targets);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromEndInBounds(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    if(distance > distance_) {
        distance = distance_;
    } else if(column_depth <= 0) {
        return 0.0;
    }
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth, targets);
    return distance;
}

double Path::GetDistanceFromStartAlongPath(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, direction_, column_depth, targets);
    return distance;
}

double Path::GetDistanceFromEndAlongPath(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, -direction_, column_depth, targets);
    return distance;
}

double Path::GetDistanceFromStartInReverse(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, first_point_, -direction_, column_depth);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double column_depth, std::set<LeptonInjector::Particle::ParticleType> const & targets) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth, targets);
    return distance;
}

double Path::GetDistanceFromEndInReverse(double column_depth) {
    EnsureIntersections();
    double distance = earth_model_->DistanceForColumnDepthFromPoint(intersections_, last_point_, -direction_, column_depth);
    return distance;
}

bool Path::IsWithinBounds(Vector3D point) {
    EnsurePoints();
    double d0 = earthmodel::scalar_product(direction_, first_point_ - point);
    double d1 = earthmodel::scalar_product(direction_, last_point_ - point);
    return d0 <= 0 and d1 >= 0;
}

