#include "plane.h"
#include "ray.h"
#include <cfloat>
#include <limits>

// Intersect with the half space defined by the plane.  The plane's normal
// points outside.  If the ray starts on the "inside" side of the plane, be sure
// to record a hit with t=0 as the first entry in hits.
Hit Plane::Intersection(const Ray& ray, int part) const
{
    //TODO;
    // float denominator = dot(normal, ray.direction);
    // if(denominator > 1e-5) {
    //     vec3 pointToPlane = x1 - ray.endpoint;
    // }
    
    return {0,0,0};
    // Hit hit;
    // double t;
    // double U;
    // vec3 W = ray.endpoint - x1; // E - x1
    // U = dot(ray.direction, normal);
    // if(U == 0){
    //     hit = {0, 0, 0}; //no intersection
    // }
    // else{
    //     t = (-1)*dot(W, normal)/U;
    //     if(t > small_t){
    //         hit = {this, t, 1};
    //     }
    //     else{
    //         hit = {0, 0, 0};
    //     }
    // }
    // return hit;
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
