#ifndef INZ_INSTANCE_H
#define INZ_INSTANCE_H

#include "MACgrid.h"

class Instance {
    inline int s(int x, int y, int z);
    inline int u(int x, int y, int z);
    inline int v(int x, int y, int z);
    inline int w(int x, int y, int z);
public:
    const int Nx, Ny, Nz, maxN;
    const float dx, dt;
    float stepTime;

    Instance(const int x, const int y, const int z, const float dx, const float dt, const int instance);

    void update(MACgrid &grid, float time);
private:
    void(Instance::*f)(MACgrid&, float);
    void empty(MACgrid&, float);
    //Cases
    void fullWater(MACgrid &grid, float time);
    void fromTap(MACgrid &grid, float time);
    void fromTap2(MACgrid &grid, float time);
    void test(MACgrid &grid, float time);
    void half(MACgrid &grid, float time);
    void narrowing(MACgrid &grid, float time);
    void narrowing2(MACgrid &grid, float time);
    void wir(MACgrid &grid, float time);
    void wir2(MACgrid &grid, float time);
};

#endif //INZ_INSTANCE_H
