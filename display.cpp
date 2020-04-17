#include "display.h"

Display::Display(const int Nx, const int Ny, const int Nz) :
        Nx(Nx),Ny(Ny),Nz(Nz),gridSize(Nx*Ny*Nz), maxN(Nz>(Ny>Nx?Ny:Nx)?Nz:(Ny>Nx?Ny:Nx)) {
    //render solid
    {
        const char *vertexShader = GLSL(
            in int gridCell;
            flat out int cell;

            uniform mat4 MVP;
            uniform int Nx;
            uniform int Ny;
            uniform int Nz;
            uniform int maxN;

            void main() {
                cell = gridCell;

                vec3 position = vec3(
                    2.0/maxN*(gl_VertexID % Nx - 0.5*(Nx-1.0)),
                    2.0/maxN*(gl_VertexID/Nx % Ny - 0.5*(Ny-1.0)),
                    2.0/maxN*(gl_VertexID/Nx/Ny % Nz - 0.5*(Nz-1.0))
                );
                gl_Position = MVP*vec4(position, 1.0);
            }
        );
        const char *fragmentShader = GLSL(
            out vec4 color;
            flat in int cell;

            void main() {
                if(cell==0) discard;
                else if(cell==2) color = vec4(0.1, 0.1, 0.1, 0.1);
                else color = vec4(0.0, 0.0, 1.0, 0.5);
            }
        );

        program = createProgram(vertexShader, fragmentShader);
        render = createVBO(program, "gridCell", 4*gridSize, 1);
    }
    //!render solid
    //render pressure
    {
        const char *vertexShader = GLSL(
            in int gridCell;
            flat out int cell;
            in float pressure;
            out float p;

            uniform mat4 MVP;
            uniform int Nx;
            uniform int Ny;
            uniform int Nz;
            uniform int maxN;

            void main() {
                cell = gridCell;
                p = pressure;

                vec3 position = vec3(
                    2.0/maxN*(gl_VertexID % Nx - 0.5*(Nx-1.0)),
                    2.0/maxN*(gl_VertexID/Nx % Ny - 0.5*(Ny-1.0)),
                    2.0/maxN*(gl_VertexID/Nx/Ny % Nz - 0.5*(Nz-1.0))
                );
                gl_Position = MVP*vec4(position, 1.0);
            }
        );
        const char *fragmentShader = GLSL(
            out vec4 color;
            flat in int cell;
            in float p;

            uniform float pMax;

            void main() {
                if(cell!=2) discard;
                else color = vec4(p/pMax, 0.0, (pMax-p)/pMax, p/pMax+0.1);
            }
        );

        programPressure = createProgram(vertexShader, fragmentShader);
        renderSolids = createVBO(programPressure, "gridCell", 4*gridSize, 1);
        renderPressure = createVBO(programPressure, "pressure", 4*gridSize, 1);
    }
    //!render pressure
    //render text
    {
        error = FT_Init_FreeType(&library);
        if(error) perror("Error: initialize FreeType library");
        error = FT_New_Face(library, "FreeSans.ttf", 0, &face);
        if(error) perror("Error: loading font");

        const char *vertexShader = GLSL(
            in int glyph;
            flat out int color;

            uniform int x;
            uniform int y;
            uniform int w;
            uniform int h;
            uniform int width;

            void main() {
                color = glyph;

                vec2 position = vec2(
                    2.0*(x + gl_VertexID%width)/w-1.0,
                    2.0*(y - gl_VertexID/width)/h-1.0
                );
                gl_Position = vec4(position, 1.0, 1.0);
            }
        );
        const char *fragmentShader = GLSL(
            flat in int color;
            out vec4 Color;

            void main() {
                if(color==0) discard;
                else Color = vec4(0.0, 0.0, 0.0, 1.0);
            }
        );

        programText = createProgram(vertexShader, fragmentShader);
        renderText = createVBO(programText, "glyph", 4*100*100, 1);
    }
    //render text
}

Display::~Display() {
    glDeleteProgram(program);
    glDeleteBuffers(1, &render);

    glDeleteProgram(programPressure);
    glDeleteBuffers(1, &renderSolids);
    glDeleteBuffers(1, &renderPressure);

    glDeleteProgram(programText);
    glDeleteBuffers(1, &renderText);
    FT_Done_Face(face);
    FT_Done_FreeType(library);
}

void Display::draw(const float *MVP, const MACgrid &grid, const float pointSize) const {
//    glEnable(GL_DEPTH_TEST);
//    glDepthFunc(GL_LESS);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1.0f,1.0f,1.0f,1.0f);

    glUseProgram(program);

    GLint position;
    position = glGetUniformLocation(program, "MVP");
    glUniformMatrix4fv(position, 1, GL_FALSE, MVP);
    position = glGetUniformLocation(program, "Nx");
    glUniform1i(position, Nx);
    position = glGetUniformLocation(program, "Ny");
    glUniform1i(position, Ny);
    position = glGetUniformLocation(program, "Nz");
    glUniform1i(position, Nz);
    position = glGetUniformLocation(program, "maxN");
    glUniform1i(position, maxN);

    glBindBuffer(GL_ARRAY_BUFFER, render);
    glVertexAttribPointer((GLuint)glGetAttribLocation(program, "gridCell"), 1, GL_FLOAT, GL_FALSE, 0, 0);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 4*gridSize, grid.solids);

    glPointSize(pointSize);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawArrays(GL_POINTS, 0, gridSize);
}

void Display::drawP(const float *MVP, const MACgrid &grid, const float pointSize, const float pMax) const {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glClearColor(1.0f,1.0f,1.0f,1.0f);

    glUseProgram(programPressure);

    GLint position;
    position = glGetUniformLocation(programPressure, "MVP");
    glUniformMatrix4fv(position, 1, GL_FALSE, MVP);
    position = glGetUniformLocation(programPressure, "Nx");
    glUniform1i(position, Nx);
    position = glGetUniformLocation(programPressure, "Ny");
    glUniform1i(position, Ny);
    position = glGetUniformLocation(programPressure, "Nz");
    glUniform1i(position, Nz);
    position = glGetUniformLocation(programPressure, "maxN");
    glUniform1i(position, maxN);
    position = glGetUniformLocation(programPressure, "pMax");
    glUniform1f(position, pMax);

    glBindBuffer(GL_ARRAY_BUFFER, renderSolids);
    glVertexAttribPointer((GLuint)glGetAttribLocation(programPressure, "gridCell"), 1, GL_FLOAT, GL_FALSE, 0, 0);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 4*gridSize, grid.solids);
    glBindBuffer(GL_ARRAY_BUFFER, renderPressure);
    glVertexAttribPointer((GLuint)glGetAttribLocation(programPressure, "pressure"), 1, GL_FLOAT, GL_FALSE, 0, 0);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 4*gridSize, grid.pressure);

    glPointSize(pointSize);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glDrawArrays(GL_POINTS, 0, gridSize);
}

void Display::print(int x, int y, const char *text, const int size, const int w, const int h) {
    FT_Set_Char_Size(face, size * 64, 0, 100, 0);

    for(const char *s=text; *s>0; s++) {
        error = FT_Load_Char(face, (FT_ULong)*s, FT_LOAD_RENDER);
        if(error) perror("Error: loading character");

        FT_Bitmap bitmap = face->glyph->bitmap;

        glUseProgram(programText);

        GLint position;
        position = glGetUniformLocation(programText, "x");
        glUniform1i(position, x + face->glyph->bitmap_left);
        position = glGetUniformLocation(programText, "y");
        glUniform1i(position, y + face->glyph->bitmap_top);
        position = glGetUniformLocation(programText, "w");
        glUniform1i(position, w);
        position = glGetUniformLocation(programText, "h");
        glUniform1i(position, h);
        position = glGetUniformLocation(programText, "width");
        glUniform1i(position, bitmap.width);

        int buffer[100*100];
        for(int i=0; i<bitmap.width*bitmap.rows; i++) buffer[i] = bitmap.buffer[i];

        glBindBuffer(GL_ARRAY_BUFFER, renderText);
        glVertexAttribPointer((GLuint)glGetAttribLocation(programText, "glyph"), 1, GL_FLOAT, GL_FALSE, 0, 0);
        glBufferSubData(GL_ARRAY_BUFFER, 0, 4*bitmap.width*bitmap.rows, buffer);

        glPointSize(1.0f);
        glDrawArrays(GL_POINTS, 0, bitmap.width*bitmap.rows);

        x += face->glyph->advance.x / 64;
        y += face->glyph->advance.y / 64;
    }
}

GLuint Display::createProgram(const char *vertex_code, const char *fragment_code) const {
    //Create shaders
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    GLint compileStatus = GL_FALSE, infoLogLength = 0, linkStatus = GL_FALSE;

    // Compile Vertex Shader
    glShaderSource(vertexShader, 1, &vertex_code , NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compileStatus);
    if(compileStatus==GL_FALSE) {
        glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string vertexShaderError(infoLogLength, 0);
        glGetShaderInfoLog(vertexShader, infoLogLength, NULL, &vertexShaderError[0]);
        glDeleteShader(vertexShader);
        throw std::runtime_error("createProgram: vertexShader: "+vertexShaderError);
    }

    // Compile Fragment Shader
    glShaderSource(fragmentShader, 1, &fragment_code , NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compileStatus);
    if(compileStatus==GL_FALSE) {
        glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string fragmentShaderError(infoLogLength, 0);
        glGetShaderInfoLog(fragmentShader, infoLogLength, NULL, &fragmentShaderError[0]);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        throw std::runtime_error("createProgram: fragmentShader: "+fragmentShaderError);
    }

    // Link the program
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glLinkProgram(program);

    glDetachShader(program, vertexShader);
    glDetachShader(program, fragmentShader);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);

    glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
    if(linkStatus==GL_FALSE) {
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string programError(infoLogLength, 0);
        glGetProgramInfoLog(program, infoLogLength, NULL, &programError[0]);
        glDeleteProgram(program);
        throw std::runtime_error("createProgram: link: "+programError);
    }

    return program;
}

GLuint Display::createProgram(const char *vertex_code, const char *fragment_code, const char *geometry_code) const {
    //Create shaders
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    GLuint geometryShader = glCreateShader(GL_GEOMETRY_SHADER);

    GLint compileStatus = GL_FALSE, infoLogLength = 0, linkStatus = GL_FALSE;

    // Compile Vertex Shader
    glShaderSource(vertexShader, 1, &vertex_code , NULL);
    glCompileShader(vertexShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &compileStatus);
    if(compileStatus==GL_FALSE) {
        glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string vertexShaderError(infoLogLength, 0);
        glGetShaderInfoLog(vertexShader, infoLogLength, NULL, &vertexShaderError[0]);
        glDeleteShader(vertexShader);
        throw std::runtime_error("createProgram: vertexShader: "+vertexShaderError);
    }

    // Compile Fragment Shader
    glShaderSource(fragmentShader, 1, &fragment_code , NULL);
    glCompileShader(fragmentShader);
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &compileStatus);
    if(compileStatus==GL_FALSE) {
        glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string fragmentShaderError(infoLogLength, 0);
        glGetShaderInfoLog(fragmentShader, infoLogLength, NULL, &fragmentShaderError[0]);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        throw std::runtime_error("createProgram: fragmentShader: "+fragmentShaderError);
    }

    //Compile Geometry Shader
    glShaderSource(geometryShader, 1, &geometry_code , NULL);
    glCompileShader(geometryShader);
    glGetShaderiv(geometryShader, GL_COMPILE_STATUS, &compileStatus);
    if(compileStatus==GL_FALSE) {
        glGetShaderiv(geometryShader, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string geometryShaderError(infoLogLength, 0);
        glGetShaderInfoLog(geometryShader, infoLogLength, NULL, &geometryShaderError[0]);
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        glDeleteShader(geometryShader);
        throw std::runtime_error("createProgram: geometryShader: "+geometryShaderError);
    }

    // Link the program
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    glAttachShader(program, geometryShader);
    glLinkProgram(program);

    glDetachShader(program, vertexShader);
    glDetachShader(program, fragmentShader);
    glDetachShader(program, geometryShader);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    glDeleteShader(geometryShader);

    glGetProgramiv(program, GL_LINK_STATUS, &linkStatus);
    if(linkStatus==GL_FALSE) {
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLogLength);
        std::string programError(infoLogLength, 0);
        glGetProgramInfoLog(program, infoLogLength, NULL, &programError[0]);
        glDeleteProgram(program);
        throw std::runtime_error("createProgram: link: "+programError);
    }

    return program;
}

GLuint Display::createVBO(GLuint program, const char *variable, int N, int dimensions) const {
    GLuint vbo;
    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, N, nullptr, GL_STATIC_DRAW);

    GLint position = glGetAttribLocation(program, variable);
    if(position<0) throw std::logic_error(std::string("createVBO: wrong variable name ")+variable);
    glEnableVertexAttribArray((GLuint)position);
    glVertexAttribPointer((GLuint)position, dimensions, GL_FLOAT, GL_FALSE, 0, 0);

    return vbo;
}
