#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

vec3 Phong_Shader::Shade_Surface(const Ray& ray,const vec3& intersection_point, 
    const vec3& normal,int recursion_depth) const
{
    //TODO Ambient lighting
    vec3 color, diffuse, specular;
    color = specular = world.ambient_color * world.ambient_intensity * color_ambient;

    //Loop through each ligh source in the world
    for(unsigned i = 0; i < world.lights.size(); i++) {
        //Vector to light
        vec3 L = world.lights[i]->position - intersection_point;
        vec3 R = (-L + 2 * dot(L, normal) * normal).normalized();

        //If shadows are enabled
        if(world.enable_shadows) {
            //Get ray for the shadow
            Ray Shadow_ray(intersection_point, L);
            //Get closest object intersected by shadow ray
            Hit Obj = world.Closest_Intersection(Shadow_ray);

            //If there is an object and the light source is not between the starting point
            //and the object, return the color without specular and diffuse modifiers.
            if(Obj.object && Obj.dist < L.magnitude()) {
                return color;
            }
        }
        //Diffuse vector
        diffuse = world.lights[i]->Emitted_Light(L) * color_diffuse * std::max(dot(normal, L.normalized()), 0.0);
        //Specular vector
        specular = world.lights[i]->Emitted_Light(L) * color_specular * std::pow(std::max(dot(R, -(ray.direction)), 0.0), specular_power);
        //Update color based on current light
        color += diffuse + specular;
    }
    return color;
}
