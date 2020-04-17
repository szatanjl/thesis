#ifndef INZ_NAVIER_H
#define INZ_NAVIER_H

#include <vector>
#include <algorithm>
#include "MACgrid.h"

class Navier {
    static int Nx, Ny, Nz;
    static inline int s(int x, int y, int z);
    static inline int u(int x, int y, int z);
    static inline int v(int x, int y, int z);
    static inline int w(int x, int y, int z);
public:
    static float step(MACgrid &grid, float dt); //return max optimal dt for next step
private:
    //advect markers + extrapolate velocity + update fluid
    static void advectMarkers(std::vector<MACgrid::Marker> &markers, float *velocity, std::vector<MACgrid::Marker> &fluid, const int *solids, const float dt, const float dx, float *tmpS);
    static void updateSolids(int *solids, const std::vector<MACgrid::Marker> &fluid);
    static void advectVelocity(float *&velocity, const float dt, const float dx, float *&tmpV);
    static void bodyForces(float *&velocity, const float *bodyForces, const int *solids, const float dt, float *&tmpV);
    static float project(float *velocity, float *pressure, const int *solids, const std::vector<MACgrid::Marker> &fluid, float *tmpV, float *tmpS); //return max velocity
};

#endif //INZ_NAVIER_H
