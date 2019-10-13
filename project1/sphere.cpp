#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    vec3 L = center - ray.endpoint;//Get vector L from center to ray endpoint
    double tc = dot(L, ray.direction);//gets dot product of L and ray direction

    if(tc < 0)//If tc < 0 then the ray is pointing away from the sphere
        return {0,0,0};
    
    //Get the length of the line segment from center that is perpendicular to the ray
    double d = sqrt((L.magnitude * L.magnitude) - (tc * tc));

    if(d > radius)//Not touching the sphere if the d is greater than the raidus
        return {0,0,0};
    
    //Calculate the points intersecting with the sphere
    double tc1 = sqrt((radius * radius) - (d * d));
    double t1 = tc - tc1;
    double t2 = tc + tc1;

    Ray p1Ray = Ray(ray.endpoint, ray.direction);
    Ray p2Ray = Ray(ray.endpoint, ray.direction);

    //Get vectors to the intersections on the sphere
    vec3 p1 = p1Ray.Point(t1);
    vec3 p2 = p2Ray.Point(t2);
    
    Hit sphereIntersection;
    sphereIntersection.object = this;
    sphereIntersection.dist = (p1.magnitude > p2.magnitude) ? p1.magnitude : p2.magnitude;
    sphereIntersection.part = part;
    return sphereIntersection;
}

vec3 Sphere::Normal(const vec3& point, int part) const
{
    vec3 normal;
    TODO; // compute the normal direction
    return normal;
}

Box Sphere::Bounding_Box(int part) const
{
    Box box;
    TODO; // calculate bounding box
    return box;
}
