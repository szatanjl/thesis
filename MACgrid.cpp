#include "MACgrid.h"

MACgrid::MACgrid(const int Nx, const int Ny, const int Nz, const float dx) :
        Nx(Nx),Ny(Ny),Nz(Nz),scalarSize(Nx*Ny*Nz),vectorSize(3*(Nx+1)*(Ny+1)*(Nz+1)),dx(dx) {
    solids            = new   int[scalarSize];

    velocity          = new float[vectorSize];
    bodyForces        = new float[vectorSize];

    pressure          = new float[scalarSize];

    density           = new float[scalarSize];
    densitySource     = new float[scalarSize];

    temperature       = new float[scalarSize];
    temperatureSource = new float[scalarSize];

    tmpV              = new float[vectorSize];
    tmpS              = new float[scalarSize];

    markers.clear();
    fluid.clear();
    for(int i=0; i<scalarSize; i++) solids[i] = MACgrid::AIR_CELL;
    for(int i=0; i<scalarSize; i++) pressure[i] = density[i] = densitySource[i] = temperature[i] = temperatureSource[i] = tmpS[i]= 0.0f;
    for(int i=0; i<vectorSize; i++) velocity[i] = bodyForces[i] = tmpV[i] = 0.0f;
}

MACgrid::~MACgrid() {
    delete[] solids;

    delete[] velocity;
    delete[] bodyForces;

    delete[] pressure;

    delete[] density;
    delete[] densitySource;

    delete[] temperature;
    delete[] temperatureSource;

    delete[] tmpV;
    delete[] tmpS;
}