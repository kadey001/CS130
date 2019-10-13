#include "render_world.h"
#include "flat_shader.h"
#include "object.h"
#include "light.h"
#include "ray.h"

extern bool disable_hierarchy;

Render_World::Render_World()
    :background_shader(0),ambient_intensity(0),enable_shadows(true),
    recursion_depth_limit(3)
{}

Render_World::~Render_World()
{
    delete background_shader;
    for(size_t i=0;i<objects.size();i++) delete objects[i];
    for(size_t i=0;i<lights.size();i++) delete lights[i];
}

// Find and return the Hit structure for the closest intersection.  Be careful
// to ensure that hit.dist>=small_t.
Hit Render_World::Closest_Intersection(const Ray& ray)
{
    Hit closestHit;
    closestHit.dist = __DBL_MAX__;
    double min_t = 1e+200;//Large value set to min_t
    for(unsigned int i = 0; i < objects.size(); i++) {
        Hit hit = objects[i]->Intersection(ray, 0);
        if(hit.dist < min_t && hit.dist > small_t) {
            closestHit = hit;
            min_t = hit.dist;
        }
    }
    return closestHit;
}

// set up the initial view ray and call
void Render_World::Render_Pixel(const ivec2& pixel_index)
{
    Ray ray;
    ray.direction = (camera.position - camera.World_Position(pixel_index)).normalized();
    ray.endpoint = camera.position;

    vec3 color=Cast_Ray(ray,1);
    camera.Set_Pixel(pixel_index,Pixel_Color(color));
}

void Render_World::Render()
{
    if(!disable_hierarchy)
        Initialize_Hierarchy();

    for(int j=0;j<camera.number_pixels[1];j++)
        for(int i=0;i<camera.number_pixels[0];i++)
            Render_Pixel(ivec2(i,j));
}

// cast ray and return the color of the closest intersected surface point,
// or the background color if there is no object intersection
vec3 Render_World::Cast_Ray(const Ray& ray,int recursion_depth)
{
    vec3 color;
    // determine the color here TODO
    Hit closestHit = Closest_Intersection(ray);
    vec3 intersectionPoint = ray.endpoint;
    vec3 norm = closestHit.object->Normal(intersectionPoint, closestHit.part);

    if(closestHit.dist < __DBL_MAX__) {//If thee is an intersection
        color = closestHit.object->material_shader->Shade_Surface(ray, intersectionPoint, norm, recursion_depth);
    } 
    else {
        color = background_shader->Shade_Surface(ray, intersectionPoint, norm, recursion_depth);
    }
    return color;
}

void Render_World::Initialize_Hierarchy()
{
    TODO; // Fill in hierarchy.entries; there should be one entry for
    // each part of each object.

    hierarchy.Reorder_Entries();
    hierarchy.Build_Tree();
}
