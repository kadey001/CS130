#include "plane.h"
#include "ray.h"
#include <cfloat>
#include <limits>

// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
Hit Plane::Intersection(const Ray& ray, int part) const
{
    Hit hit{0, 0, -1};

    //Check if ray is perpendicular to plane (check == 0). If it is then return empty hit.
    double check = dot(ray.direction, normal);
    if(!check) {
        return hit;
    }

    //If distance t is less than the threshold of small_t, return empty hit.
    double t = dot((x1 - ray.endpoint), normal) / dot(ray.direction, normal);
    if(t < small_t) {
        return hit;
    }

    //Otherwise set object as hit with distance t
    hit.object = this;
    hit.dist = t;
    return hit;
}

vec3 Plane::Normal(const vec3& point, int part) const
{
    return normal;
}

// There is not a good answer for the bounding box of an infinite object.
// The safe thing to do is to return a box that contains everything.
Box Plane::Bounding_Box(int part) const
{
    Box b;
    b.hi.fill(std::numeric_limits<double>::max());
    b.lo=-b.hi;
    return b;
}
