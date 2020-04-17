#include "instance.h"

inline int Instance::s(int x, int y, int z) {
    return (x) + (y)*Nx + (z)*Nx*Ny;
}
inline int Instance::u(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 0;
}
inline int Instance::v(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 1;
}
inline int Instance::w(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 2;
}

Instance::Instance(const int Nx, const int Ny, const int Nz, const float dx, const float dt, const int instance) :
        Nx(Nx),Ny(Ny),Nz(Nz),dx(dx),dt(dt), maxN(Nz>(Ny>Nx?Ny:Nx)?Nz:(Ny>Nx?Ny:Nx)) {
    if(instance==0) f = &Instance::fullWater;
    if(instance==1) f = &Instance::fromTap;
    if(instance==2) f = &Instance::test;
    if(instance==3) f = &Instance::half;
    if(instance==4) f = &Instance::narrowing;
    if(instance==5) f = &Instance::wir;
}

void Instance::update(MACgrid &grid, float time) {
    (this->*f)(grid, time);
}

void Instance::empty(MACgrid&, float) {}

//Cases

void Instance::fullWater(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -9.8f;
    }

    //walls
    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,0,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,Nz-1)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    //water
    for(int x=1; x<Nx-1; x++) for(int y=1; y<Ny-1; y++) for(int z=1; z<Nz-1; z++) {
        grid.markers.push_back(MACgrid::Marker(x + 0.01f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.99f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.01f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.99f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.01f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.99f,z + 0.5f));

        grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
    }

    f = &Instance::empty;
}

void Instance::fromTap(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -1.0f;
    }

    //walls
    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,0,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny/2; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,Nz-1)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny/2; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    stepTime = 0.0f;
    f = &Instance::fromTap2;
}

void Instance::fromTap2(MACgrid &grid, float time) {
    if(time < stepTime) return;

    stepTime += 1.0f;
    int x = 5, y = 10, z = 5;
//    grid.markers.push_back(MACgrid::Marker(x + 0.01f,y + 0.5f,z + 0.5f));
//    grid.markers.push_back(MACgrid::Marker(x + 0.99f,y + 0.5f,z + 0.5f));
//    grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.01f));
//    grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.99f));
//    grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.01f,z + 0.5f));
//    grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.99f,z + 0.5f));

    float ii, jj, kk;
    for(int i=0; i<3; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) {
        if(i==0) ii = 0.01f;
        else if(i==1) ii = 0.5f;
        else ii = 0.99f;
        if(j==0) jj = 0.01f;
        else if(j==1) jj = 0.5f;
        else jj = 0.99f;
        if(k==0) kk = 0.01f;
        else if(k==1) kk = 0.5f;
        else kk = 0.99f;

        grid.markers.push_back(MACgrid::Marker(x+ii,y+jj,z+kk));
    }

    grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
}

void Instance::test(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -9.8f;
    }

    //walls
    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,0,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,2)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    grid.velocity[u(2,1,1)] = grid.velocity[u(3,1,1)] = 50.0f;

    //water
    {
        int x = 2, y = 1, z = 1;
        grid.markers.push_back(MACgrid::Marker(x + 0.01f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.99f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.01f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.99f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.01f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.99f,z + 0.5f));

        grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
    }

    f = &Instance::empty;
}

void Instance::half(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -1.0f;
    }

    //walls
    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,0,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,Nz-1)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    //water
    for(int x=1; x<Nx/3-1; x++) for(int y=1; y<25; y++) for(int z=1; z<Nz-1; z++) {
        float ii, jj, kk;
        for(int i=0; i<3; i++) for(int j=0; j<3; j++) for(int k=0; k<3; k++) {
            if(i==0) ii = 0.01f;
            else if(i==1) ii = 0.5f;
            else ii = 0.99f;
            if(j==0) jj = 0.01f;
            else if(j==1) jj = 0.5f;
            else jj = 0.99f;
            if(k==0) kk = 0.01f;
            else if(k==1) kk = 0.5f;
            else kk = 0.99f;

            grid.markers.push_back(MACgrid::Marker(x+ii,y+jj,z+kk));
        }

        grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
    }

    f = &Instance::empty;
}

void Instance::narrowing(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=35; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -1.0f;
    }

    //walls
    for(int x=0; x<15; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,20,z)] = MACgrid::SOLID_CELL;
    for(int x=25; x<40; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,20,z)] = MACgrid::SOLID_CELL;


    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,39,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,Nz-1)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    stepTime = 0.0f;
    f = &Instance::narrowing2;
}

void Instance::narrowing2(MACgrid &grid, float time) {
    if(time < stepTime) return;

    stepTime += 1.0f;
    for(int x=1; x<Nx-1; x++) for(int y=37; y<38; y++) for(int z=1; z<Nz-1; z++) {
        grid.markers.push_back(MACgrid::Marker(x + 0.01f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.99f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.01f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.99f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.01f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.99f,z + 0.5f));
        grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
    }
}

void Instance::wir(MACgrid &grid, float) {
    //gravity
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        grid.bodyForces[v(x,y,z)] = -9.8f;
    }

    float C = 10.0f;
    //ruch wirowy
    for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
        if(x==20 && z==20) continue;
        grid.bodyForces[u(x,y,z)] = C*(z-20)*(z-20)/((x-20)*(x-20)+(z-20)*(z-20)) * (z-20>0 ? -1 : 1);
        grid.bodyForces[w(x,y,z)] = C*(x-20)*(x-20)/((x-20)*(x-20)+(z-20)*(z-20)) * (x-20>0 ? 1 : -1);
    }

    //walls
    for(int x=0; x<Nx; x++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(x,0,z)] = MACgrid::SOLID_CELL;

    for(int x=0; x<Nx; x++)
        for(int y=0; y<Ny; y++)
            grid.solids[s(x,y,0)] = grid.solids[s(x,y,Nz-1)] = MACgrid::SOLID_CELL;

    for(int y=0; y<Ny; y++)
        for(int z=0; z<Nz; z++)
            grid.solids[s(0,y,z)] = grid.solids[s(Nx-1,y,z)] = MACgrid::SOLID_CELL;

    //water
    for(int x=1; x<Nx-1; x++) for(int y=1; y<10; y++) for(int z=1; z<Nz-1; z++) {
        grid.markers.push_back(MACgrid::Marker(x + 0.01f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.99f,y + 0.5f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.01f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.5f,z + 0.99f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.01f,z + 0.5f));
        grid.markers.push_back(MACgrid::Marker(x + 0.5f,y + 0.99f,z + 0.5f));

        grid.solids[s(x,y,z)] = MACgrid::FLUID_CELL;
    }

    f = &Instance::wir2;
}

void Instance::wir2(MACgrid &grid, float time) {
    if(time>20.0f) {
        for(int x=0; x<=Nx; x++) for(int y=0; y<=Ny; y++) for(int z=0; z<=Nz; z++) {
            grid.bodyForces[u(x,y,z)] = grid.bodyForces[w(x,y,z)] = 0.0f;
        }

        f = &Instance::empty;
    }
}