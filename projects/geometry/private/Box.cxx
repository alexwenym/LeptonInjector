#include <map>
#include <cmath>
#include <tuple>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "LeptonInjector/math/Vector3D.h"
#include "LeptonInjector/geometry/Geometry.h"
#include "LeptonInjector/geometry/Placement.h"
#include "LeptonInjector/geometry/Box.h"

namespace LI {
namespace geometry {

Box::Box()
    : Geometry((std::string)("Box"))
    , x_(0.0)
    , y_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Box::Box(double x, double y, double z)
    : Geometry("Box")
    , x_(x)
    , y_(y)
      , z_(z)
{
    // Do nothing here
}

Box::Box(Placement const & placement)
    : Geometry((std::string)("Box"), placement)
    , x_(0.0)
    , y_(0.0)
      , z_(0.0)
{
    // Do nothing here
}

Box::Box(Placement const & placement, double x, double y, double z)
    : Geometry((std::string)("Box"), placement)
    , x_(x)
    , y_(y)
      , z_(z)
{
    // Do nothing here
}

Box::Box(const Box& box)
    : Geometry(box)
    , x_(box.x_)
    , y_(box.y_)
      , z_(box.z_)
{
    // Nothing to do here
}

/*Box::Box(const nlohmann::json& config)
  : Geometry(config)
  {
  if(not config.at("length").is_number())
  throw std::invalid_argument("Length is not a number.");
  if(not config.at("width").is_number())
  throw std::invalid_argument("Width is not a number.");
  if(not config.at("height").is_number())
  throw std::invalid_argument("Height is not a number.");

  x_ = config["length"].get<double>();
  y_ = config["width"].get<double>();
  z_ = config["height"].get<double>();

  if(x_ < 0) throw std::logic_error("lenght must be > 0");
  if(y_ < 0) throw std::logic_error("width must be > 0");
  if(z_ < 0) throw std::logic_error("height must be > 0");
  }*/

// ------------------------------------------------------------------------- //
void Box::swap(Geometry& geometry)
{
    Box* box = dynamic_cast<Box*>(&geometry);
    if (!box)
    {
        //log_warn("Cannot swap Box!");
        return;
    }

    Geometry::swap(*box);

    std::swap(x_, box->x_);
    std::swap(y_, box->y_);
    std::swap(z_, box->z_);
}

//------------------------------------------------------------------------- //
Box& Box::operator=(const Geometry& geometry)
{
    if (this != &geometry)
    {
        const Box* box = dynamic_cast<const Box*>(&geometry);
        if (!box)
        {
            //log_warn("Cannot assign Sphere!");
            return *this;
        }

        Box tmp(*box);
        swap(tmp);
    }
    return *this;
}

// ------------------------------------------------------------------------- //
bool Box::equal(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    if (!box)
        return false;
    else if (x_ != box->x_)
        return false;
    else if (y_ != box->y_)
        return false;
    else if (z_ != box->z_)
        return false;
    else
        return true;
}

// ------------------------------------------------------------------------- //
bool Box::less(const Geometry& geometry) const
{
    const Box* box = dynamic_cast<const Box*>(&geometry);

    return
        std::tie(x_, y_, z_)
        <
        std::tie(box->x_, box->y_, box->z_);
}

// ------------------------------------------------------------------------- //
void Box::print(std::ostream& os) const
{
    os << "Width_x: " << x_ << "\tWidth_y " << y_ << "\tHeight: " << z_ << '\n';
}

// ------------------------------------------------------------------------- //
std::vector<Geometry::Intersection> Box::ComputeIntersections(LI::math::Vector3D const & position, LI::math::Vector3D const & direction) const {
    // Calculate intersection of particle trajectory and the box
    // Surface of the box is defined by six planes:
    // E1: x1   =   position.GetX() + 0.5*x
    // E2: x1   =   position.GetX() - 0.5*x
    // E3: x2   =   position.GetY() + 0.5*y
    // E4: x2   =   position.GetY() - 0.5*y
    // E5: x3   =   position.GetZ() + 0.5*z
    // E6: x3   =   position.GetZ() - 0.5*z
    // straight line (particle trajectory) g = vec(x,y,z) + t * dir_vec( cosph
    // *sinth, sinph *sinth , costh)
    // We are only interested in postive values of t
    // ( we want to find the intersection in direction of the particle
    // trajectory)

    double dir_vec_x = direction.GetX();
    double dir_vec_y = direction.GetY();
    double dir_vec_z = direction.GetZ();

    double t;
    double intersection_x;
    double intersection_y;
    double intersection_z;
    bool entering;

    std::vector<Intersection> dist;

    std::function<void()> save = [&](){
        Intersection i;
        i.position = LI::math::Vector3D(intersection_x,intersection_y,intersection_z);
        i.distance = t;
        i.hierarchy = 0;
        i.entering = entering;
        dist.push_back(i);
    };

    double x_calc_pos =   0.5 * x_;
    double x_calc_neg = - 0.5 * x_;
    double y_calc_pos =   0.5 * y_;
    double y_calc_neg = - 0.5 * y_;
    double z_calc_pos =   0.5 * z_;
    double z_calc_neg = - 0.5 * z_;

    // intersection with E1
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E1
    {
        t = (x_calc_pos - position.GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_y = position.GetY() + t * dir_vec_y;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_x = position.GetX() + t * dir_vec_x;
            entering = direction.GetX() < 0;
            save();
        }
    }

    // intersection with E2
    if (dir_vec_x != 0) // if dir_vec == 0 particle trajectory is parallel to E2
    {
        t = (x_calc_neg - position.GetX()) / dir_vec_x;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_y = position.GetY() + t * dir_vec_y;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_y >= y_calc_neg && intersection_y <= y_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_x = position.GetX() + t * dir_vec_x;
            entering = direction.GetX() > 0;
            save();
        }
    }

    // intersection with E3
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E3
    {
        t = (y_calc_pos - position.GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_y = position.GetY() + t * dir_vec_y;
            entering = direction.GetY() < 0;
            save();
        }
    }

    // intersection with E4
    if (dir_vec_y != 0) // if dir_vec == 0 particle trajectory is parallel to E4
    {
        t = (y_calc_neg - position.GetY()) / dir_vec_y;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_z = position.GetZ() + t * dir_vec_z;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_z >= z_calc_neg &&
                intersection_z <= z_calc_pos)
        {
            intersection_y = position.GetY() + t * dir_vec_y;
            entering = direction.GetY() > 0;
            save();
        }
    }

    // intersection with E5
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E5
    {
        t = (z_calc_pos - position.GetZ()) / dir_vec_z;

        // Computer precision controll
        if (std::fabs(t) < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            entering = direction.GetZ() < 0;
            save();
        }
    }

    // intersection with E6
    if (dir_vec_z != 0) // if dir_vec == 0 particle trajectory is parallel to E6
    {
        t = (z_calc_neg - position.GetZ()) / dir_vec_z;

        // Computer precision controll
        if (t > 0 && t < GEOMETRY_PRECISION)
            t = 0;

        // Check if intersection is inside the box borders
        intersection_x = position.GetX() + t * dir_vec_x;
        intersection_y = position.GetY() + t * dir_vec_y;
        if (intersection_x >= x_calc_neg && intersection_x <= x_calc_pos && intersection_y >= y_calc_neg &&
                intersection_y <= y_calc_pos)
        {
            intersection_z = position.GetZ() + t * dir_vec_z;
            entering = direction.GetZ() > 0;
            save();
        }
    }

    std::function<bool(Intersection const &, Intersection const &)> comp = [](Intersection const & a, Intersection const & b){
        return a.distance < b.distance;
    };

    std::sort(dist.begin(), dist.end(), comp);
    return dist;
}

// ------------------------------------------------------------------------- //
std::pair<double, double> Box::ComputeDistanceToBorder(const LI::math::Vector3D& position, const LI::math::Vector3D& direction) const
{
    // Compute the surface intersections
    std::vector<Intersection> intersections = Intersections(position, direction);
    std::vector<double> dist;
    for(unsigned int i=0; i<intersections.size(); ++i) {
        if(intersections[i].distance > 0) {
            dist.push_back(intersections[i].distance);
        }
    }

    std::pair<double, double> distance;

    if (dist.size() < 1) // No intersection with the box
    {
        distance.first  = -1;
        distance.second = -1;
    } else if (dist.size() == 1) // Particle is inside the box and we have one
        // intersection in direction of the particle
        // trajectory
    {
        distance.first  = dist.at(0);
        distance.second = -1;
    } else if (dist.size() == 2) // Particle is outside and the box is infront
        // of the particle trajectory ( two
        // intersections).
    {
        distance.first  = dist.at(0);
        distance.second = dist.at(1);
        if (distance.second < distance.first)
        {
            std::swap(distance.first, distance.second);
        }

    } else
    {
        //log_error("This point should nerver be reached... (-1/-1) is returned");

        distance.first  = -1;
        distance.second = -1;
    }

    // Make a computer precision controll!
    // This is necessary cause due to numerical effects it meight be happen
    // that a particle which is located on a gemoetry border is treated as
    // inside
    // or outside
    if (distance.first < GEOMETRY_PRECISION)
        distance.first = -1;
    if (distance.second < GEOMETRY_PRECISION)
        distance.second = -1;
    if (distance.first < 0)
        std::swap(distance.first, distance.second);

    return distance;
}

} // namespace geometry
} // namespace LI
