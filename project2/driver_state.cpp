#include "driver_state.h"
#include "common.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color= new pixel[width * height];
    state.image_depth=0;
    for(size_t i = 0; i < width * height; i++) {
        state.image_color[i] = make_pixel(0,0,0);
    }
    std::cout<<"TODO: allocate and initialize state.image_depth."<<std::endl;
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    switch(type) {
        case render_type::triangle: {
            data_geometry* triangle_geometry = new data_geometry[3];
            int num_of_triangles = state.num_vertices / 3;
            int index = 0;
            for(size_t i = 0; i < num_of_triangles; i++) {
                for(size_t j = 0; j < 3; j++) {
                    std::cout << *(state.vertex_data + index) << std::endl;
                    triangle_geometry[j].data = state.vertex_data + index;
                    std::cout << *triangle_geometry[j].data << std::endl;
                    index += state.floats_per_vertex;
                }
            }
            rasterize_triangle(state, (const data_geometry**)&triangle_geometry);
            break;
        }
        case render_type::indexed: {
            break;
        }
        case render_type::fan: {
            break;
        }
        case render_type::strip: {
            break;
        }
        default: {
            std::cout << "Invalid Type" << std::endl;
        }
    }
    std::cout<<"TODO: implement rendering."<<std::endl;
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state,in,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout << "Raster:" << *(*in)[0].data << std::endl;
    // for(size_t i = 0; i < 3; i++) {
    //     std::cout << *(*in)[i].data/state.image_width << std::endl;
    //     //state.vertex_shader((*in)[i].data, )
    // }
    std::cout<<"TODO: implement rasterization"<<std::endl;
}

