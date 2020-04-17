#include <cstdio>
#include <ctime>
#include <string>

#define GLM_FORCE_RADIANS
#include <glm/gtx/transform.hpp>

#include "OpenGL.h"
#include "MACgrid.h"
#include "instance.h"
#include "display.h"
#include "navier.h"

//Global variables for callbacks
bool enableRotate = false;
bool pause = false;
bool pressure = false;
double xOrig, yOrig;
//MVP
float fov = 90.0f;
float width = 1024.0f, height = 768.0f;
glm::mat4 View = glm::lookAt(
    glm::vec3(0.0f,0.0f,3.0f),
    glm::vec3(0.0f,0.0f,0.0f),
    glm::vec3(0.0f,1.0f,0.0f)
);
glm::mat4 rotateX = glm::mat4(1.0f);
float rotateY = 0.0f;
//!Global variables

//Callbacks
void key_callback(GLFWwindow *window, int key, int, int action, int) {
    if(action == GLFW_PRESS || action == GLFW_REPEAT) {
        //ESC Close
        if(key == GLFW_KEY_ESCAPE) glfwSetWindowShouldClose(window, GL_TRUE);

        //Reset to default
        if(key == GLFW_KEY_DELETE) {
            fov = 90.0f;
            rotateX = glm::mat4(1.0f);
            rotateY = 0.0f;
        }

        //Pause
        if(key == GLFW_KEY_SPACE) pause ^= true;

        //Render pressure/solids
        if(key == GLFW_KEY_P) pressure ^= true;
    }
}

void window_size_callback(GLFWwindow*, int _width, int _height) { //Resize window
    glViewport(0, 0, _width, _height);
    width = _width;
    height = _height;
}

void mouse_button_callback(GLFWwindow *window, int button, int action, int) {
    //Enable rotate
    if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
        enableRotate = true;
        glfwGetCursorPos(window, &xOrig, &yOrig);
    }
    //Disable rotate
    if(button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
        enableRotate = false;
    }
}

void cursor_pos_callback(GLFWwindow*, double x, double y) {
    //Rotate
    if(enableRotate) {
        rotateX = glm::rotate(glm::radians((float)(x - xOrig)/2.0f), glm::vec3(0.0f, 1.0f, 0.0f)) * rotateX;
        if(rotateY + (float)(y - yOrig)/2.0f < 90.0f && rotateY + (float)(y - yOrig)/2.0f > -90.0f)
            rotateY += (float)(y - yOrig)/2.0f;

        xOrig = x;
        yOrig = y;
    }
}

void scroll_callback(GLFWwindow*, double, double y) {
    //Zoom
    if(y > 0.0) if(fov - 2.0f > 30.0f)  fov -= 2.0f;
    if(y < 0.0) if(fov + 2.0f < 180.0f) fov += 2.0f;
}
//!Callbacks

int main() {
    glfwSetErrorCallback([](int, const char *desc) {perror(desc);});

    if(!glfwInit()) {
        perror("Failed to initialize GLFW\n");
        return -1;
    }

    glfwWindowHint(GLFW_SAMPLES, 4);                               //4x antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);                 //We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);           //To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); //We don't want the old OpenGL

    GLFWwindow *window = glfwCreateWindow((int)width, (int)height, "Fluid simulation", NULL, NULL);
    if(!window) {
        perror("Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n" );
        glfwTerminate();
        return -1;
    }

    glfwMakeContextCurrent(window); //Initialize GLEW
    glfwSwapInterval(1);

    glfwSetKeyCallback(window, key_callback);
    glfwSetWindowSizeCallback(window, window_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_pos_callback);
    glfwSetScrollCallback(window, scroll_callback);

    glewExperimental = (GLboolean)true; //Needed in core profile
    if(glewInit() != GLEW_OK) {
        perror("Failed to initialize GLEW\n");
        return -1;
    }

    {
        GLuint vao;
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        //instance numbers:
        //pelna woda stojaca w miejscu = 0
        //woda lejaca z kranu = 1
        //jakis test = 2
        //polowa wody = 3
        //zwezenie = 4
        //wir = 5

        Instance instance(40,40,40, 1.0f,0.01f, 5); //grid size, dx, dt, instance number
        MACgrid grid(instance.Nx, instance.Ny, instance.Nz, instance.dx);
        Display display(instance.Nx, instance.Ny, instance.Nz);

        long long delta=0, frames=0;
        double fps = 0.0;
        float maxP = 0.0f;

        for(float time=0.0f, dt=instance.dt; !glfwWindowShouldClose(window);) {
            if(!pause) {
                delta -= clock();

                instance.update(grid, time);

                dt = Navier::step(grid, dt);
                time += (dt = instance.dt<=0.0f ? dt : instance.dt);

                delta += clock();
                frames++;
            }
            {
                if(!pressure) display.draw(                                             //Render fluid
                    &(glm::perspective(glm::radians(fov), width/height, 0.1f, 100.0f) * //Projection
                    View *                                                              //View
                    glm::rotate(glm::radians(rotateY), glm::vec3(1.0f,0.0f,0.0f)) *     //Model
                    rotateX)[0][0],
                    grid,                                                               //MACgrid
                    height/2.0f/(float)tan(glm::radians(fov/2))/instance.maxN           //Grid cell size
                );
                else display.drawP(                                                     //Render pressure
                    &(glm::perspective(glm::radians(fov), width/height, 0.1f, 100.0f) * //Projection
                    View *                                                              //View
                    glm::rotate(glm::radians(rotateY), glm::vec3(1.0f,0.0f,0.0f)) *     //Model
                    rotateX)[0][0],
                    grid,                                                               //MACgrid
                    height/2.0f/(float)tan(glm::radians(fov/2))/instance.maxN,          //Grid cell size
                    3800.0f/10.0f*(instance.dt<=0.0f ? dt : instance.dt)/instance.dx    //max pressure in hPa/(rho/100)*dt/dx
                );

                display.print(50, 700, ("FPS: "+std::to_string(fps)).c_str(), 20, (int)width, (int)height);
                display.print(50, 670, ("Czas: "+std::to_string(time)+" s").c_str(), 20, (int)width, (int)height);
                display.print(50, 640, ("Max P: "+std::to_string(maxP*instance.dx/instance.dt*10.0f)+" hPa").c_str(), 20, (int)width, (int)height);

                float gora = 0.0f, dol = 0.0f;
                for(int y=0; y<22; y++) dol += grid.velocity[3*((20) + (y)*(40+1) + (20)*(40+1)*(40+1)) + 1]; dol /= 22.0f;
                for(int y=22; y<instance.Ny; y++) gora += grid.velocity[3*((20) + (y)*(40+1) + (20)*(40+1)*(40+1)) + 1]; gora /= 17.0f;

//                display.print(50, 600, ("Avg vel at top: "+std::to_string(gora)).c_str(), 20, (int)width, (int)height); //For instance 4
//                display.print(50, 570, ("Avg vel at bot: "+std::to_string(dol)).c_str(), 20, (int)width, (int)height);
            }
            if(delta>CLOCKS_PER_SEC/4) {
                fps = frames/((double)delta/CLOCKS_PER_SEC);
                delta = frames = 0;

                maxP = 0.0f;
                for(int i=0; i<grid.Nx*grid.Ny*grid.Nz; i++) if(grid.solids[i]==MACgrid::SOLID_CELL && maxP<grid.pressure[i]) maxP = grid.pressure[i];
            }
            glfwSwapBuffers(window);
            glfwPollEvents();
        }

        glDeleteVertexArrays(1, &vao);
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
