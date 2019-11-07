#include <limits>
#include "box.h"

// Return whether the ray intersects this box.
bool Box::Intersection(const Ray& ray) const
{
    //TODO
    int x = 0;
    int y = 1;
    int z = 2;

    //X plane
    double x_min = (lo[x] - ray.endpoint[x]) / ray.direction[x];
    double x_max = (hi[x] - ray.endpoint[x]) / ray.direction[x];
    if(x_max < x_min) {
        double temp = x_max;
        x_max = x_min;
        x_min = temp;
    }

    //Y plane
    double y_min = (lo[y] - ray.endpoint[y]) / ray.direction[y];
    double y_max = (hi[y] - ray.endpoint[y]) / ray.direction[y];
    if (y_max < y_min)
    {
        double temp = y_max;
        y_max = y_min;
        y_min = temp;
    }

    //Z plane
    double z_min = (lo[z] - ray.endpoint[z]) / ray.direction[z];
    double z_max = (hi[z] - ray.endpoint[z]) / ray.direction[z];
    if (z_max < z_min)
    {
        double temp = z_max;
        z_max = z_min;
        z_min = temp;
    }

    double min = 0.0;
    double max = 0.0;
    if(x_min > y_min)//Greatest minimum x/y
        min = x_min;
    if(x_max < y_max)//Smallest maximum x/y
        max = x_max;

    if(x_min > y_max || y_min > x_max)
        return false;
    if(min > z_max || z_min > max)
        return false;
        
    return true;
}

// Compute the smallest box that contains both *this and bb.
Box Box::Union(const Box& bb) const
{
    Box box;
    TODO;
    return box;
}

// Enlarge this box (if necessary) so that pt also lies inside it.
void Box::Include_Point(const vec3& pt)
{
    //Find if point is inside of box
    vec3 box_difference = hi - lo;
    vec3 point_difference = hi - pt;

    int x = 0;
    int y = 1;
    int z = 2;
    bool too_small = true;

    while(too_small) {
        bool x_diff = point_difference[x] <= box_difference[x];
        bool y_diff = point_difference[y] <= box_difference[y];
        bool z_diff = point_difference[z] <= box_difference[z];
        
        if (x_diff && y_diff && z_diff) {
            too_small = false;
        } else {
            this->lo[x] += 1.0;
            this->lo[y] += 1.0;
            this->lo[z] += 1.0;
            this->hi[x] += 1.0;
            this->hi[y] += 1.0;
            this->hi[z] += 1.0;
        }
    }
    return;
}

// Create a box to which points can be correctly added using Include_Point.
void Box::Make_Empty()
{
    lo.fill(std::numeric_limits<double>::infinity());
    hi=-lo;
}
