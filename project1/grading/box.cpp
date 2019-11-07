#include <limits>
#include "box.h"

using namespace std;

// Return whether the ray intersects this box.
bool Box::Intersection(const Ray& ray) const
{
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
    int x = 0;
    int y = 1;
    int z = 2;

    //Get box hi
    if(lo[x] > bb.lo[x])
        box.hi[x] = lo[x];
    else 
        box.hi[x] = bb.lo[x];
    if(lo[y] > bb.lo[y])
        box.hi[y] = lo[y];
    else
        box.hi[y] = bb.lo[y];
    if (lo[z] > bb.lo[z])
        box.hi[z] = lo[z];
    else
        box.hi[z] = bb.lo[z];

    //Get box lo
    if (lo[x] < bb.lo[x])
        box.lo[x] = lo[x];
    else
        box.lo[x] = bb.lo[x];
    if (lo[y] < bb.lo[y])
        box.lo[y] = lo[y];
    else
        box.lo[y] = bb.lo[y];
    if (lo[z] < bb.lo[z])
        box.lo[z] = lo[z];
    else
        box.lo[z] = bb.lo[z];

    return box;
}

// Enlarge this box (if necessary) so that pt also lies inside it.
void Box::Include_Point(const vec3& pt)
{
    int x = 0;
    int y = 1;
    int z = 2;
    
    //X coordinate
    if(pt[x] > hi[x])
        hi[x] = pt[x];
    if(pt[x] < lo[x])
        lo[x] = pt[x];

    //Y coordinate
    if(pt[y] > hi[y])
        hi[y] = pt[y];
    if(pt[y] < lo[y])
        lo[y] = pt[y];

    //Z coordinate
    if(pt[z] > hi[z])
        hi[z] = pt[z];
    if(pt[z] < lo[z])
        lo[z] = pt[z];
}

// Create a box to which points can be correctly added using Include_Point.
void Box::Make_Empty()
{
    lo.fill(std::numeric_limits<double>::infinity());
    hi=-lo;
}
