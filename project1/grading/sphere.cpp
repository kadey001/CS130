#include "sphere.h"
#include "ray.h"

// Determine if the ray intersects with the sphere
Hit Sphere::Intersection(const Ray& ray, int part) const
{
    vec3 L = center - ray.endpoint;//Get vector L from center to ray endpoint
    double tc = dot(L, ray.direction);//gets dot product of L and ray direction
    if(tc < 0)//If tc < 0 then the ray is pointing away from the sphere
        return {0,0,0};

    //Get the length of the line segment from center to point p that is perpendicular to the ray
    vec3 p = ray.Point(tc);
    double d = (center - p).magnitude();

    if(d > radius)//Not touching the sphere if the d is greater than the raidus
        return {0,0,0};
    
    //Calculate the points intersecting with the sphere
    double tc1 = sqrt((radius * radius) - (d * d));
    double t1 = tc - tc1;
    double t2 = tc + tc1;
    Hit hit;

    if(t2 > t1) {
        hit = {this, t1, 1};
    }
    else {
        hit = {this, t2, 1};
    }

    return hit;

    //Ray p1Ray = Ray(ray.endpoint, ray.direction);
    //Ray p2Ray = Ray(ray.endpoint, ray.direction);

    //Get vectors to the intersections on the sphere
    //vec3 p1 = p1Ray.Point(t1);
    //vec3 p2 = p2Ray.Point(t2);
    
    // Hit sphereIntersection;
    // sphereIntersection.object = this;
    // sphereIntersection.dist = (p1.magnitude() > p2.magnitude()) 
    //     ? p1.magnitude() : p2.magnitude();
    // sphereIntersection.part = 0;//part;
    // return sphereIntersection;

////////////////////////////////////////////////////////////////////
    
    // Hit hit;
    // vec3 W = ray.endpoint - center;
    // double a, b, c, t1, t2, root;
    // a = ray.direction.magnitude_squared();
    // b = 2*dot(ray.direction, W);
    // c = W.magnitude_squared() - (radius*radius);
    // root = (b*b) - (4*a*c);
    // if(root < 0){
    //     hit = {0, 0, 0};
    // }
    // else{
    //     t1 = ((-1)*b + sqrt(root))/(2*a);
    //     t2 = ((-1)*b - sqrt(root))/(2*a);
    //     if(t1 < t2 && t1 >= small_t){ //t1 is closer to camera
    //         hit = {this, t1, 1};
    //     }
    //     else if(t2 <= t1 && t2 >= small_t){ //t2 is closer to camera
    //         hit = {this, t2, 1};
    //     }
    //     else{ //no intersect
    //         hit = {0, 0, 0};
    //     }
    // }    
    // return hit;
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
