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
    state.image_depth = new float[width * height];
    for(size_t i = 0; i < width * height; i++) {
        state.image_color[i] = make_pixel(0,0,0);
        state.image_depth[i] = 1;
    }
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
            const data_geometry* triangle_geometry[3];
            data_geometry geometry_d[3];
            data_vertex vertex_d[3];
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
                    vertex_d[i].data = state.vertex_data + index;
                    geometry_d[i].data = vertex_d[i].data;
                    //std::cout << *triangle_geometry[i].data << std::endl;

                    //Fill gl_Position for the triangle geometry allowing it to be accessed while clipping
                    state.vertex_shader(vertex_d[i], geometry_d[i], state.uniform_data);

                    triangle_geometry[i] = &geometry_d[i];
                    //std::cout << triangle_geometry[i].gl_Position << std::endl;

                    index += state.floats_per_vertex;
                }
                //rasterize_triangle(state, (const data_geometry**)&triangle_geometry);
                clip_triangle(state, triangle_geometry, 0);
            }
            break;
        }
        case render_type::indexed: {
            const data_geometry* triangle_geometry[3];
            data_geometry geometry_d[3];
            data_vertex vertex_d[3];

            for(size_t index = 0; index < state.num_triangles * 3; index += 3) {
                for(size_t i = 0; i < 3; i++) {
                    vertex_d[i].data = &state.vertex_data[state.index_data[index + i] * state.floats_per_vertex];
                    geometry_d[i].data = vertex_d[i].data;
                    state.vertex_shader(vertex_d[i], geometry_d[i], state.uniform_data);
                    triangle_geometry[i] = &geometry_d[i];
                }
                clip_triangle(state, triangle_geometry, 0);
            }
            break;
        }
        case render_type::fan: {
            const data_geometry* triangle_geometry[3];
            data_geometry geometry_d[3];
            data_vertex vertex_d[3];
            int temp = 0;

            for(size_t fan = 0; fan < state.num_vertices; fan++) {
                for(size_t i = 0; i < 3; i++) {
                    if(i == 0) {
                        temp = 0;
                    } else {
                        temp = fan + i;
                    }

                    vertex_d[i].data = &state.vertex_data[temp * state.floats_per_vertex];
                    geometry_d[i].data = vertex_d[i].data;
                    state.vertex_shader(vertex_d[i], geometry_d[i], state.uniform_data);
                    triangle_geometry[i] = &geometry_d[i];
                }
                clip_triangle(state, triangle_geometry, 0);
            }
            break;
        }
        case render_type::strip: {
            const data_geometry* triangle_geometry[3];
            data_geometry geometry_d[3];
            data_vertex vertex_d[3];

            for(size_t strip = 0; strip < state.num_vertices - 2; strip++) {
                for(size_t i = 0; i < 3; i++) {
                    vertex_d[i].data = &state.vertex_data[(strip + i) * state.floats_per_vertex];
                    geometry_d[i].data = vertex_d[i].data;
                    state.vertex_shader(vertex_d[i], geometry_d[i], state.uniform_data);
                    triangle_geometry[i] = &geometry_d[i];
                }
                clip_triangle(state, triangle_geometry, 0);
            }
            break;
        }
        default: {
            std::cout << "Invalid Type" << std::endl;
        }
    }
}

void interpolate_clipping(driver_state& state, const data_geometry* in[3], data_geometry(&geometry_d) [3], float& lerp, int& interp_case) {
    enum {x, y, z, w};
    enum {a, b, c};

    switch(interp_case) {
        case 0:
            for(size_t i = 0; i < state.floats_per_vertex; i++) {
                //Update A<->C based on interpolation type
                switch(state.interp_rules[i]) {
                    case interp_type::flat:
                        geometry_d[a].data[i] = in[a]->data[i];
                        break;
                    case interp_type::noperspective:{
                        float temp = lerp * in[c]->gl_Position[w] / (lerp * in[c]->gl_Position[w] + (1 - lerp) * in[a]->gl_Position[w]);
                        //std::cout << temp << std::endl;
                        geometry_d[a].data[i] = temp * in[c]->data[i] + (1 - temp) * in[a]->data[i];
                        //std::cout << geometry_d1[a].data[i] << std::endl;
                        break;
                    }
                    case interp_type::smooth:
                        geometry_d[a].data[i] = lerp * in[c]->data[i] + (1 - lerp) * in[a]->data[i];
                        //std::cout << geometry_d1[a].data[i] << std::endl;
                    default:
                        break;
                }
            }
            break;
        case 1:
            for(size_t i = 0; i < state.floats_per_vertex; i++) {
            //Update A<->B based on interpolation type
                switch(state.interp_rules[i]) {
                    case interp_type::flat:
                        geometry_d[a].data[i] = in[a]->data[i];
                        break;
                    case interp_type::noperspective:{
                        float temp = lerp * in[a]->gl_Position[w] / (lerp * in[a]->gl_Position[w] + (1 - lerp) * in[b]->gl_Position[w]);
                        geometry_d[a].data[i] = temp * in[a]->data[i] + (1 - temp) * in[b]->data[i];
                        break;
                    }
                    case interp_type::smooth:
                        geometry_d[a].data[i] = lerp * in[a]->data[i] + (1 - lerp) * in[b]->data[i];
                    default:
                        break;
                }
            }
            break;
        default:
            break;
    }
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3], int face)
{
    if(face == 6) {
        rasterize_triangle(state, in);
        return;
    }

    enum {x, y, z, w};
    enum {a, b, c};

    const data_geometry* triangle_geometry[3] = {in[a], in[b], in[c]};

    vec4 A = in[a]->gl_Position;
    vec4 B = in[b]->gl_Position;
    vec4 C = in[c]->gl_Position;

    int interp_case = 0;
    float lerp_1, lerp_2;
    vec4 P1, P2;

    data_geometry geometry_d1[3], geometry_d2[3];
    geometry_d1[a].data = new float[state.floats_per_vertex];
    geometry_d1[b] = *in[b];
    geometry_d1[c] = *in[c];

    geometry_d2[a].data = new float[state.floats_per_vertex];
    geometry_d2[b] = *in[b];
    geometry_d2[c] = *in[c];

    //Check if clipping is needed
    if(A[z] < -A[w] && B[z] < -B[w] && C[z] < -C[w]) {
        //All points outside so no need to clip
        return;
    } else if(A[z] < -A[w] && B[z] >= -B[w] && C[z] >= -C[w]) {
        //Clip triangle and update interpolation
        //A is outside
        //Find lerp factors
        lerp_1 = (-A[z] - A[w]) / (C[z] + C[w] - A[z] - A[w]);
        lerp_2 = (-B[z] - B[w]) / (A[z] + A[w] - B[z] - B[w]);

        //std::cout << lerp_1 << " " << lerp_2 << std::endl;

        P1 = lerp_1 * C + (1 - lerp_1) * A;
        P2 = lerp_2 * A + (1 - lerp_2) * B;

        //std::cout << P1 << " " << P2 << std::endl;

        //Interpolate geometry data
        interpolate_clipping(state, in, geometry_d1, lerp_1, interp_case);

        geometry_d1[a].gl_Position = P1;
        triangle_geometry[a] = &geometry_d1[a];
        triangle_geometry[b] = &geometry_d1[b];
        triangle_geometry[c] = &geometry_d1[c];

        clip_triangle(state, triangle_geometry, face + 1);

        interp_case = 1;

        //Interpolate geometry data
        interpolate_clipping(state, in, geometry_d2, lerp_2, interp_case);

        geometry_d2[a].gl_Position = P2;
        triangle_geometry[a] = &geometry_d2[a];
        triangle_geometry[b] = &geometry_d1[b];
        triangle_geometry[c] = &geometry_d1[a];
    }

    clip_triangle(state, triangle_geometry, face + 1);
}

void interpolate_pixels(driver_state& state, const data_geometry* in[3], data_output& output_d, data_fragment& fragment_d, float& alpha, float& beta, float& gamma, float(&w) [3]) {
    enum {a, b, c};
    const float ALPHA = alpha;
    const float BETA = beta;
    const float GAMMA = gamma;
    
    //Interpolation
    for (int k = 0; k < state.floats_per_vertex; k++) {
        float temp;
        switch (state.interp_rules[k]) {
            case interp_type::flat:
                //Flat fragment shader
                fragment_d.data[k] = in[a]->data[k];
                break;

            case interp_type::noperspective:
                //interpolate using image-space barycentric coordinates
                fragment_d.data[k] = (
                    alpha * in[a]->data[k] + 
                    beta * in[b]->data[k] + 
                    gamma * in[c]->data[k]
                );
                break;

            case interp_type::smooth:
                //Perspective interpolation
                temp = (ALPHA / w[0]) + (BETA / w[1]) + (GAMMA / w[2]);
                
                alpha = ALPHA / temp / w[0];
                beta = BETA / temp / w[1];
                gamma = GAMMA / temp / w[2];

                fragment_d.data[k] = (
                    alpha * in[a]->data[k] + 
                    beta * in[b]->data[k] + 
                    gamma * in[c]->data[k]
                );
                break;

            default:
                break;
        }
    }
}


// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    enum {a, b, c};
    //How to access vertex data:
    //Find pointer to data from: (*in)[vertex].data + (0 for x axis, 1 for y axis)
    //Or use in[vertex]->data[axis]
    //Or in[vertex]->gl_Position[axis]
    //Must dereference address to get values
    //std::cout << *(*in)[1].data << std::endl;
    //std::cout << in[1]->data << std::endl;
    //std::cout << in[1]->gl_Position[0];

    //Get the w values for each vertex from geometry_d
    float w0 = in[0]->gl_Position[3];
    float w1 = in[1]->gl_Position[3];
    float w2 = in[2]->gl_Position[3];
    float w[] = {w0, w1, w2};
    //std::cout << w0 << w1 << w2 << std::endl;

    //Get each vertex coordinates divided by respective w
    float x0 = in[0]->gl_Position[0] / w0;
    float x1 = in[1]->gl_Position[0] / w1;
    float x2 = in[2]->gl_Position[0] / w2;
    float x[] = {x0, x1, x2};

    float y0 = in[0]->gl_Position[1] / w0;
    float y1 = in[1]->gl_Position[1] / w1;
    float y2 = in[2]->gl_Position[1] / w2;
    float y[] = {y0, y1, y2};

    float z0 = in[0]->gl_Position[2] / w0;
    float z1 = in[1]->gl_Position[2] / w1;
    float z2 = in[2]->gl_Position[2] / w2;

    //std::cout << "X:" << x0 << " " << x1 << " " << x2 << std::endl;
    //std::cout << "Y:" << y0 << " " << y1 << " " << y2 << std::endl;
    //std::cout << "Z:" << z0 << " " << z1 << " " << z2 << std::endl;

    int width = state.image_width;
    int height = state.image_height;
    float pixel_x[3], pixel_y[3];
    float i, j;
    float middle_w = width / 2.0;
    float middle_h = height / 2.0;
    
    for(size_t iter = 0; iter < 3; iter++) {
        i = middle_w * x[iter] + middle_w - .5;
        j = middle_h * y[iter] + middle_h - .5;
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
    //zero or the corresponding maximum image width/height.
    if(minimum_x < 0){
        minimum_x = 0;
    }
    if(maximum_x > width){
        maximum_x = width;
    }
    if(minimum_y < 0){
        minimum_y = 0;
    }
    if(maximum_y > height){
        maximum_y = height;
    }

    //Calculate area of triangle to use barycentric coordinates
    //Area(abc) = .5 * ((BxCy - CxBy)-(AxCy - CxAy)-(AxBy - BxAy))
    //Broke up into parts and used enum for readability.
    //(Parts also useful in barycentric coordinate calculations)
    float part1 = (pixel_x[b] * pixel_y[c]) - (pixel_x[c] * pixel_y[b]);
    float part2 = (pixel_x[a] * pixel_y[c]) - (pixel_x[c] * pixel_y[a]);
    float part3 = (pixel_x[a] * pixel_y[b]) - (pixel_x[b] * pixel_y[a]);
    float triangle_area = .5 * (part1 - part2 + part3);
    //std::cout << triangle_area << std::endl;
    // for(unsigned i = 0; i < 3; i++) {
    //     std::cout << pixel_x[i] << std::endl;
    //     std::cout << pixel_y[i] << std::endl;
    // }

    float alpha, beta, gamma;

    //std::cout << "X:" << minimum_x << maximum_x << std::endl;
    //std::cout << "Y:" << minimum_y << maximum_y << std::endl;    

    //Find the barycentric coords for each pixel in the area of the triangle. If the coordinate is inside
    //the triangle then we update the color, otherwise leave it black.
    //We only need to check from the minimum x,y pixels to maximum x,y pixels since we know the triangle is in this area.
    //Interpolate points within the triangle based on interp_rules.
    for(int i = minimum_x; i < maximum_x; i++) {
        for(int j = minimum_y; j < maximum_y; j++) {
            alpha = (.5 * (part1 + (pixel_y[b] - pixel_y[c]) * i + (pixel_x[c] - pixel_x[b]) * j) / triangle_area);
            beta = (.5 * (((pixel_x[c] * pixel_y[a]) - (pixel_x[a] * pixel_y[c])) + (pixel_y[c] - pixel_y[a]) * i + (pixel_x[a] - pixel_x[c]) * j) / triangle_area);
            gamma = (.5 * (part3 + (pixel_y[a] - pixel_y[b]) * i + (pixel_x[b] - pixel_x[a]) * j) / triangle_area); 

            //Find depth for Z-Buffering
            float depth = alpha * z0 + beta * z1 + gamma * z2;

            //Only render pixel when the point is inside of the triangle and 
            if(alpha >= 0 && beta >= 0 && gamma >= 0 && state.image_depth[i + j * width] > depth) {
                data_output output_d;
                data_fragment fragment_d;
                fragment_d.data = new float[state.floats_per_vertex];
                
                //Interpolation
                interpolate_pixels(state, in, output_d, fragment_d, alpha, beta, gamma, w);
                
                state.fragment_shader(fragment_d, output_d, state.uniform_data);
                //std::cout << output_d.output_color << std::endl;

                state.image_depth[i + j * width] = depth;
                state.image_color[i + j * width] = make_pixel(output_d.output_color[0] * 255, output_d.output_color[1] * 255, output_d.output_color[2] * 255);
            }
        }
    }
}
