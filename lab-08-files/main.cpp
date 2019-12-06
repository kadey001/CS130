// Name:
// Quarter, Year:
// Lab:
//
// This file is to be modified by the student.
// main.cpp
////////////////////////////////////////////////////////////
#ifndef __APPLE__
#include <GL/glut.h>
#else
#include <GLUT/glut.h>
#endif

#include <vector>
#include <cstdio>
#include <math.h>
#include "vec.h"
#include <iostream>

using namespace std;
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 800;

vector<vec2> controlPoints;

vec2 binomial(int n, int k, float t);

void GL_render()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glutSwapBuffers();

    glBegin(GL_LINE_STRIP);
    glColor3f(1.0f,0.0f,0.0f);
    // just for example, remove if desired
    size_t n = controlPoints.size();
    for(float t = 0; t <= 1; t += .01) {
        vec2 tempVec;
        tempVec = binomial(n, 0, t);
        glVertex2f(tempVec.x[0], tempVec.x[1]);
    }
    //std::cout << controlPoints[0] << std::endl;
    //glVertex2f(binomial(controlPoints.size(), 1, 1), 1.0f);
    //glVertex2f(.5f, .5f);
    //glVertex2f(-.5f,-.5f);
    // glVertex2f(.5f,.5f);
    // glVertex2f(-.5f,.5f);
    glEnd();
    glFlush();
}

void GL_mouse(int button,int state,int x,int y)
{
    y=WINDOW_HEIGHT-y;
    GLdouble mv_mat[16];
    GLdouble proj_mat[16];
    GLint vp_mat[4];
    glGetDoublev(GL_MODELVIEW_MATRIX,mv_mat);
    glGetDoublev(GL_PROJECTION_MATRIX,proj_mat);
    glGetIntegerv(GL_VIEWPORT,vp_mat);

    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN){
        double px,py,dummy_z; // we don't care about the z-value but need something to pass in
        gluUnProject(x,y,0,mv_mat,proj_mat,vp_mat,&px,&py,&dummy_z);
        glutPostRedisplay();
        controlPoints.push_back(vec2(px, py));
    }
}

float factorial(int n) {
    if (n <= 1)
        return(1);
    else
        n = n * factorial(n - 1);
    return n;
}

float combination(int n, int k) {
    return (factorial(n) / (factorial(k) * factorial(n - k)));
}

vec2 binomial(int n, int k, float t) {
    vec2 result;
    for(int i = 0; i <= n; i++) {
        result[0] += (combination(n, i) * pow(t, i) * pow(1 - t, n - i) * controlPoints[i][0]);
        result[1] += (combination(n, i) * pow(t, i) * pow(1 - t, n - i) * controlPoints[i][1]);
    }

    return result;
}

//Initializes OpenGL attributes
void GLInit(int* argc, char** argv)
{
    glutInit(argc, argv);
    //glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
    glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);

    //glMatrixMode(GL_PROJECTION_MATRIX);
    //glOrtho(0, WINDOW_WIDTH, 0, WINDOW_HEIGHT, -1, 1);
    glutCreateWindow("CS 130 - Kelton");
    glutDisplayFunc(GL_render);
    glutMouseFunc(GL_mouse);
}

int main(int argc, char** argv)
{
    GLInit(&argc, argv);
    glutMainLoop();
    return 0;
}
