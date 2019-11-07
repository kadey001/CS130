#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    //vector from center to ray endpoint
    vec3 L = center - ray.endpoint;
    //product of L and ray direction
    double tc = dot(L, ray.direction);
    //Hit
    Hit hit{0,0,0};

    //If tc < 0 then the ray is pointing away from the sphere
    if(tc < 0)
        return hit;//Empty Hit

    //Get the length of the line segment from center to point p that is perpendicular to the ray
    vec3 p = ray.Point(tc);
    double d = (center - p).magnitude();

    if(d > radius)//Not touching the sphere if the d is greater than the raidus
        return hit;//Empty Hit
    
    //Calculate the points intersecting with the sphere
    double tc1 = sqrt((radius * radius) - (d * d));
    double t1 = tc - tc1;
    double t2 = tc + tc1;

    if(t2 > t1) 
        hit = {this, t1, 1};
    else
        hit = {this, t2, 1};
    
    return hit;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal = point - center;
    normal = normal.normalized();
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    return box;
}
