#ifndef INZ_MACGRID_H
#define INZ_MACGRID_H

#include <vector>

class MACgrid {
public:
    static constexpr int AIR_CELL = 0;
    static constexpr int FLUID_CELL = 1;
    static constexpr int SOLID_CELL = 2;

    struct Marker {
        float x, y, z;
        Marker(float a, float b, float c) : x(a),y(b),z(c) {}
    };

    const float dx;
    const int Nx, Ny, Nz;
    const int scalarSize;
    const int vectorSize;

    std::vector<Marker> markers;
    std::vector<Marker> fluid;
    int *solids;

    float *velocity, *bodyForces;
    float *pressure;
    float *density, *densitySource;
    float *temperature, *temperatureSource;

    float *tmpV, *tmpS;

    MACgrid(const int Nx, const int Ny, const int Nz, const float dx);
    ~MACgrid();
};


#endif //INZ_MACGRID_H
