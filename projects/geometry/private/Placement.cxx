#include <cmath>
#include <tuple>
#include <iostream>

#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/math/EulerQuaternionConversions.h"

using namespace LI::geometry;

//----------------------------------------------------------------------//
//------------------------- Constructors -------------------------------//
//----------------------------------------------------------------------//

Placement::Placement() :
    position_(LI::math::Vector3D(0,0,0)),
    quaternion_(LI::math::QFromZXZr(0,0,0))
{
    quaternion_.normalize();
}

Placement::Placement(LI::math::Vector3D const & position) :
    position_(position)
{
}

Placement::Placement(LI::math::Quaternion const & quaternion) :
    quaternion_(quaternion)
{
    quaternion_.normalize();
}

Placement::Placement(LI::math::Vector3D const & position, LI::math::Quaternion const & quaternion) :
    position_(position),
    quaternion_(quaternion)
{
    quaternion_.normalize();
}

Placement::Placement(const Placement& placement) :
    position_(placement.position_),
    quaternion_(placement.quaternion_)
{
}

Placement::Placement(Placement&& other) :
    position_(std::move(other.position_)),
    quaternion_(std::move(other.quaternion_))
{
}

//destructor
Placement::~Placement() {};

//----------------------------------------------------------------------//
//-----------------------operator functions and swap--------------------//
//----------------------------------------------------------------------//


Placement& Placement::operator=(Placement const & placement) {
    if (this != &placement)
    {
        Placement tmp(placement);
        swap(tmp);
    }
    return *this;
}

Placement& Placement::operator=(Placement && other) {
    position_ = std::move(other.position_);
    quaternion_ = std::move(other.quaternion_);
    return *this;
}

Placement& Placement::operator=(Placement const && other) {
    position_ = std::move(other.position_);
    quaternion_ = std::move(other.quaternion_);
    return *this;
}


bool Placement::operator==(const Placement& placement) const
{
    return (this == &placement) or (
            position_ == placement.position_ and
            quaternion_ == placement.quaternion_);
}

bool Placement::operator!=(const Placement& placement) const
{
    return !(*this == placement);
}

bool Placement::operator<(const Placement& placement) const
{
    return (this != &placement) and
        std::tie(position_, quaternion_)
        <
        std::tie(placement.position_, placement.quaternion_);
}

void Placement::swap(Placement& placement)
{
    using std::swap;

    swap(position_, placement.position_);
    swap(quaternion_, placement.quaternion_);
}

namespace LI {
namespace geometry {
    std::ostream& operator<<(std::ostream& os, Placement const& placement) {
        os << "Placement (" << &placement << ")" << std::endl;
        os << placement.position_ << std::endl;
        os << placement.quaternion_ << std::endl;
        return os;
    }
} // namespace geometry
} // namespace LI

std::shared_ptr<const Placement> Placement::create() const { return std::shared_ptr<const Placement>( new Placement(*this) ); }

//-------------------------------------//
// getter and setter functions

LI::math::Vector3D Placement::GetPosition() const { return position_; }

LI::math::Quaternion Placement::GetQuaternion() const { return quaternion_; }

void Placement::SetPosition(LI::math::Vector3D const & p) { position_ = p; }

void Placement::SetQuaternion(LI::math::Quaternion const & q) {
    quaternion_ = q;
    quaternion_.normalize();
}

LI::math::Vector3D Placement::Rotate(LI::math::Vector3D const & p0, bool inv) const
{
    return quaternion_.rotate(p0,inv);
}

LI::math::Vector3D Placement::LocalToGlobalPosition(LI::math::Vector3D const & p0) const
{
    LI::math::Vector3D p1 = quaternion_.rotate(p0, false); // Rotate about local origin to get orientation in the global system
    LI::math::Vector3D p2 = p1 + position_; // Translate local origin to global origin
    return p2;
}

LI::math::Vector3D Placement::LocalToGlobalDirection(LI::math::Vector3D const & p0) const
{
    LI::math::Vector3D p1 = quaternion_.rotate(p0, false); // Rotate about local origin to get orientation in the global system
    return p1;
}

LI::math::Vector3D Placement::GlobalToLocalPosition(LI::math::Vector3D const & p0) const
{
    LI::math::Vector3D p1 = p0 - position_; // Translate global origin to local origin
    LI::math::Vector3D p2 = quaternion_.rotate(p1, true); // Inverse rotatation about local origin to get orientation in the local system
    return p2;
}

LI::math::Vector3D Placement::GlobalToLocalDirection(LI::math::Vector3D const & p0) const
{
    LI::math::Vector3D p1 = quaternion_.rotate(p0, false); // Inverse rotatation about local origin to get orientation in the local system
    return p1;
}

