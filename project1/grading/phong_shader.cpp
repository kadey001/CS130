#include "light.h"
#include "phong_shader.h"
#include "ray.h"
#include "render_world.h"
#include "object.h"

using namespace std;

vec3 Phong_Shader::Shade_Surface(const Ray& ray,const vec3& intersection_point, 
    const vec3& normal,int recursion_depth) const
{
    vec3 color, ambient, diffuse, specular;

    //Get ambient color and initilize color/specular with it
    ambient = world.ambient_color * world.ambient_intensity * color_ambient;
    color = ambient;
    specular = ambient;

    //Loop through each ligh source in the world
    for(unsigned i = 0; i < world.lights.size(); i++) {
        //Vector to light, reflection
        vec3 L = world.lights[i]->position - intersection_point;
        vec3 R = (2 * dot(L, normal) * normal - L).normalized();
        vec3 emitted_light = world.lights[i]->Emitted_Light(L);

        //Check if shadows are enabled
        if(world.enable_shadows) {
            //Fire shadow ray to check for an object inbetween the light.
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
        diffuse = emitted_light * color_diffuse * max(dot(normal, L.normalized()), (double) 0);
        //Specular vector
        specular = emitted_light * color_specular * pow(max(dot(R, -ray.direction), (double) 0), specular_power);
        //Update color with diffuse and specular lighting
        color += (diffuse + specular);
    }
    return color;
}
