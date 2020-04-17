#ifndef INZ_DISPLAY_H
#define INZ_DISPLAY_H

#include <string>
#include <stdexcept>
#include <cstdio>
#include <freetype2/ft2build.h>
#include FT_FREETYPE_H

#include "OpenGL.h"
#include "MACgrid.h"

#define GLSL(src) "#version 330 core\n" #src

class Display {
    const int Nx, Ny, Nz, maxN, gridSize;
    GLuint render, renderText, renderSolids, renderPressure;  //vbo
    GLuint program, programText, programPressure; //shader programs
    FT_Library library;
    FT_Face face;
    FT_Error error;
public:
    Display(const int Nx, const int Ny, const int Nz);
    ~Display();
    void draw(const float *MVP, const MACgrid &grid, const float pointSize) const;
    void drawP(const float *MVP, const MACgrid &grid, const float pointSize, const float pMax) const;
    void print(int x, int y, const char *text, const int size, const int w, const int h);
private:
    GLuint createProgram(const char *vertex_code, const char *fragment_code) const;
    GLuint createProgram(const char *vertex_code, const char *fragment_code, const char *geometry_code) const;
    GLuint createVBO(GLuint program, const char *variable, int N, int dimensions) const;
};

#endif //INZ_DISPLAY_H
