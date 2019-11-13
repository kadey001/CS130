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
    state.image_width = width;
    state.image_height = height;
    state.image_color = new pixel[width * height];
    state.image_depth = 0;
    for(size_t i = 0; i < width * height; i++) {
        state.image_color[i] = make_pixel(0,0,0);
    }
    //std::cout<<"TODO: allocate and initialize state.image_depth."<<std::endl;
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
            //If the render is a triangle, the number of triangles should be
            //The number of verticies divided by 3.
            state.num_triangles = state.num_vertices / 3;
            int index = 0;
            //std::cout << state.num_triangles << std::endl;
            //For each triangle, copy each vertex by setting data_geometry data pointer 
            //to the first entry of the vertex data, allowing the x,y coordinates to be 
            //found for each vertex by looking at next location in memory.
            for(size_t triangle = 0; triangle < state.num_triangles; triangle++) {
                for(size_t i = 0; i < 3; i++) {
                    //std::cout << *(state.vertex_data + index) << std::endl;
                    triangle_geometry[i].data = state.vertex_data + index;
                    //std::cout << *triangle_geometry[j].data << std::endl;
                    index += state.floats_per_vertex;
                    
                }
                rasterize_triangle(state, (const data_geometry**)&triangle_geometry);
            }
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
    //std::cout<<"TODO: implement rendering."<<std::endl;
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
    //How to access vertex data:
    //Find pointer to data from: (*in)[vertex].data + (0 for x axis, 1 for y axis)
    //Must dereference address to get values
    data_vertex vertex_d[3];
    data_geometry geometry_d[3];
    for(size_t i = 0; i < 3; i++) {
        vertex_d[i].data = (*in)[i].data;
        geometry_d[i] = (*in)[i];
        state.vertex_shader(vertex_d[i], geometry_d[i], state.uniform_data);
    }
    //Get the w values for each vertex from geometry_d
    float w0 = geometry_d[0].gl_Position[3];
    float w1 = geometry_d[1].gl_Position[3];
    float w2 = geometry_d[2].gl_Position[3];
    //std::cout << w0 << w1 << w2 << std::endl;

    //Get each vertex coordinates divided by respective w
    float x0 = geometry_d[0].gl_Position[0] / w0;
    float x1 = geometry_d[1].gl_Position[0] / w1;
    float x2 = geometry_d[2].gl_Position[0] / w2;
    float x[] = {x0, x1, x2};

    float y0 = geometry_d[0].gl_Position[1] / w0;
    float y1 = geometry_d[1].gl_Position[1] / w1;
    float y2 = geometry_d[2].gl_Position[1] / w2;
    float y[] = {y0, y1, y2};

    float z0 = geometry_d[0].gl_Position[2] / w0;
    float z1 = geometry_d[0].gl_Position[2] / w1;
    float z2 = geometry_d[0].gl_Position[2] / w2;
    float z[] = {z0, z1, z2};

    //std::cout << "X:" << x0 << " " << x1 << " " << x2 << std::endl;
    //std::cout << "Y:" << y0 << " " << y1 << " " << y2 << std::endl;
    //std::cout << "Z:" << z0 << " " << z1 << " " << z2 << std::endl;

    int w = state.image_width;
    int h = state.image_height;
    int pixel_x[3], pixel_y[3];
    int i, j;
    int middle_w = (int)w/2;
    int middle_h = (int)h/2;
    
    for(size_t iter = 0; iter < 3; iter++) {
        i = (int)(middle_w * x[iter] + middle_w - .5);
        j = (int)(middle_h * y[iter] + middle_h - .5);
        pixel_x[iter] = i;
        pixel_y[iter] = j;
        //state.image_color[i + j * w] = make_pixel(255, 255, 255);
    }

    //Find the mins and maxes of triangle to check if it's in bounds and to reduce
    //pixels to check when rendering triangle.
    float minimum_x = std::min(pixel_x[0], std::min(pixel_x[1], pixel_x[2]));
    float maximum_x = std::max(pixel_x[0], std::max(pixel_x[1], pixel_x[2]));
    float minimum_y = std::min(pixel_y[0], std::min(pixel_y[1], pixel_y[2]));
    float maximum_y = std::max(pixel_y[0], std::max(pixel_y[1], pixel_y[2]));
    
    //Check to make sure that triangle is in bounds using max and min verticies
    //If the minimum vertex is less than zero or the maximum vertex is larger than
    //it's coresponding dimention then it is out of bounds and we need to set it to
    //zero;
    if(minimum_x < 0){
        minimum_x = 0;
    }
    if(maximum_x > w){
        maximum_x = w;
    }
    if(minimum_y < 0){
        minimum_y = 0;
    }
    if(maximum_y > h){
        maximum_y = h;
    }

    //Calculate area of triangle to use barycentric coordinates
    //Area(abc) = .5 * ((BxCy - CxBy)-(AxCy - CxAy)-(AxBy - BxAy))
    //Broke up into parts and used enum for readability.
    //(Parts also useful in barycentric coordinate calculations)
    enum {a, b, c};
    float part1 = (pixel_x[b] * pixel_y[c]) - (pixel_x[c] * pixel_y[b]);
    float part2 = (pixel_x[a] * pixel_y[c]) - (pixel_x[c] * pixel_y[a]);
    float part3 = (pixel_x[a] * pixel_y[b]) - (pixel_x[b] * pixel_y[a]);
    float triangle_area = .5 * (part1 - part2 - part3);
    //std::cout << triangle_area << std::endl;

    float alpha, beta, gamma;

    //std::cout << "X:" << minimum_x << maximum_x << std::endl;
    //std::cout << "Y:" << minimum_y << maximum_y << std::endl;

    //Find the barycentric coords for each pixel in the area of the triangle. If the coordinate is inside
    //the triangle then we update the color, otherwise leave it black.
    //We only need to check from the minimum x,y pixels to maximum x,y pixels since we know the triangle is in this area.
    for(int i = minimum_x; i < maximum_x; i++) {
        for(int j = minimum_y; j < maximum_y; j++) {
            alpha = (.5 * (part1 + (pixel_y[b] - pixel_y[c]) * i + (pixel_x[c] - pixel_x[b]) * j)/triangle_area);
            beta = (.5 * (((pixel_x[c] * pixel_y[a]) - (pixel_x[a] * pixel_y[c])) + (pixel_y[c] - pixel_y[a]) * i + (pixel_x[a] - pixel_x[c]) * j)/triangle_area);
            gamma = (.5 * (part3 + (pixel_y[a] - pixel_y[b]) * i + (pixel_x[b] - pixel_x[a]) * j)/triangle_area);

            //std::cout << alpha << beta << gamma << std::endl;

            if(alpha >= 0 && beta >= 0 && gamma >= 0) {
                state.image_color[i + j * w] = make_pixel(255,255,255);
            }
        }
    }

    //std::cout<<"TODO: implement rasterization"<<std::endl;
}

