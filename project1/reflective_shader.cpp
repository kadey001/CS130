#include "reflective_shader.h"
#include "ray.h"
#include "render_world.h"

vec3 Reflective_Shader::
Shade_Surface(const Ray& ray,const vec3& intersection_point,
    const vec3& normal,int recursion_depth) const
{
    //Initilize color to surface color (determined by phong shader)
    vec3 color = shader->Shade_Surface(ray, intersection_point, normal, recursion_depth);
    
    //Calculate reflection direction and build ray from it and the intersection point
    vec3 reflection_direction = ray.direction + 2 * dot((-1.0) * ray.direction, normal) * normal;
    Ray reflection_ray(intersection_point, reflection_direction);

    //Our color is determined by (1-reflectivity) * color + reflectivity * cast_ray(reflected ray, recurrsion_depth + 1)
    if(recursion_depth != world.recursion_depth_limit) {
        return color = (1 - reflectivity) * color + reflectivity * world.Cast_Ray(reflection_ray, recursion_depth + 1);
    } else {
        return (1 - reflectivity) * color;
    }
}
