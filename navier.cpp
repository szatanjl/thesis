#include "navier.h"

int Navier::Nx, Navier::Ny, Navier::Nz;

inline int Navier::s(int x, int y, int z) {
    return (x) + (y)*Nx + (z)*Nx*Ny;
}
inline int Navier::u(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 0;
}
inline int Navier::v(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 1;
}
inline int Navier::w(int x, int y, int z) {
    return 3*((x) + (y)*(Nx+1) + (z)*(Nx+1)*(Ny+1)) + 2;
}

float Navier::step(MACgrid &grid, float dt) {
    Nx = grid.Nx; Ny = grid.Ny; Nz = grid.Nz;
    if(dt<=0.0f) dt = 0.0f;
    advectMarkers(grid.markers, grid.velocity, grid.fluid, grid.solids, dt, grid.dx, grid.tmpS);
    advectVelocity(grid.velocity, dt, grid.dx, grid.tmpV);
    updateSolids(grid.solids, grid.markers);
    bodyForces(grid.velocity, grid.bodyForces, grid.solids, dt, grid.tmpV);
    float maxV = project(grid.velocity, grid.pressure, grid.solids, grid.fluid, grid.tmpV, grid.tmpS);

    return grid.dx/(maxV+10.0f);
}

void Navier::advectMarkers(std::vector<MACgrid::Marker> &M, float *U, std::vector<MACgrid::Marker> &fluid, const int *S, const float dt, const float dx, float *tmp) {
    bool remove = false;

    for(int i=0; i<Nx*Ny*Nz; i++) tmp[i] = 0.0f;

    fluid.clear();

    for(int xx=0; xx<Nx; xx++) for(int yy=0; yy<Ny; yy++) for(int zz=0; zz<Nz; zz++) {
        if(S[s(xx,yy,zz)]==MACgrid::AIR_CELL) {
            if(S[s(xx-1,yy,zz)]==MACgrid::AIR_CELL) U[u(xx,yy,zz)]   = 0.0f;
            if(S[s(xx+1,yy,zz)]==MACgrid::AIR_CELL) U[u(xx+1,yy,zz)] = 0.0f;
            if(S[s(xx,yy-1,zz)]==MACgrid::AIR_CELL) U[v(xx,yy,zz)]   = 0.0f;
            if(S[s(xx,yy+1,zz)]==MACgrid::AIR_CELL) U[v(xx,yy+1,zz)] = 0.0f;
            if(S[s(xx,yy,zz-1)]==MACgrid::AIR_CELL) U[w(xx,yy,zz)]   = 0.0f;
            if(S[s(xx,yy,zz+1)]==MACgrid::AIR_CELL) U[w(xx,yy,zz+1)] = 0.0f;
        }
    }

    for(auto &m : M) {
        int x = (int)m.x, y = (int)m.y, z = (int)m.z;
        //advect markers
        float a = m.x-x, b = m.y-y, c = m.z-z;
        m.x += ((1.0f-a)*U[u(x,y,z)] + a*U[u(x+1,y,z)])*dt/dx;
        m.y += ((1.0f-b)*U[v(x,y,z)] + b*U[v(x,y+1,z)])*dt/dx;
        m.z += ((1.0f-c)*U[w(x,y,z)] + c*U[w(x,y,z+1)])*dt/dx;
        //!advect
        //remove markers
        int xx = (int)m.x, yy = (int)m.y, zz = (int)m.z;
        if(xx<=0 || yy<=0 || zz<=0 || xx>=Nx-1 || yy>=Ny-1 || zz>=Nz-1) {
            m.x=-1.0f;
            remove = true;
        }
        //...remove
        else {
            //update fluid
            if(tmp[s(xx,yy,zz)]<0.5f) fluid.push_back(MACgrid::Marker(m.x,m.y,m.z));
            //!update
            //extrapolate velocity
            if(S[s(xx,yy,zz)]==MACgrid::AIR_CELL) {
                if(S[s(xx-1,yy,zz)]==MACgrid::AIR_CELL) U[u(xx,yy,zz)]   += U[u(x,y,z)];
                if(S[s(xx+1,yy,zz)]==MACgrid::AIR_CELL) U[u(xx+1,yy,zz)] += U[u(x+1,y,z)];
                if(S[s(xx,yy-1,zz)]==MACgrid::AIR_CELL) U[v(xx,yy,zz)]   += U[v(x,y,z)];
                if(S[s(xx,yy+1,zz)]==MACgrid::AIR_CELL) U[v(xx,yy+1,zz)] += U[v(x,y+1,z)];
                if(S[s(xx,yy,zz-1)]==MACgrid::AIR_CELL) U[w(xx,yy,zz)]   += U[w(x,y,z)];
                if(S[s(xx,yy,zz+1)]==MACgrid::AIR_CELL) U[w(xx,yy,zz+1)] += U[w(x,y,z+1)];
            }
            tmp[s(xx,yy,zz)] += 1.0f;
        }
    }

    for(int x=1; x<Nx; x++) for(int y=1; y<Ny; y++) for(int z=1; z<Nz; z++) {
        if(tmp[s(x,y,z)]>0.5f && S[s(x,y,z)]==MACgrid::AIR_CELL) {
            if(S[s(x-1,y,z)]==MACgrid::AIR_CELL) U[u(x,y,z)] /= tmp[s(x,y,z)]+tmp[s(x-1,y,z)];
            if(S[s(x+1,y,z)]==MACgrid::AIR_CELL && tmp[s(x+1,y,z)]<0.5f) U[u(x+1,y,z)] /= tmp[s(x,y,z)];
            if(S[s(x,y-1,z)]==MACgrid::AIR_CELL) U[v(x,y,z)] /= tmp[s(x,y,z)]+tmp[s(x,y-1,z)];
            if(S[s(x,y+1,z)]==MACgrid::AIR_CELL && tmp[s(x,y+1,z)]<0.5f) U[v(x,y+1,z)] /= tmp[s(x,y,z)];
            if(S[s(x,y,z-1)]==MACgrid::AIR_CELL) U[w(x,y,z)] /= tmp[s(x,y,z)]+tmp[s(x,y,z-1)];
            if(S[s(x,y,z+1)]==MACgrid::AIR_CELL && tmp[s(x,y,z+1)]<0.5f) U[w(x,y,z+1)] /= tmp[s(x,y,z)];
        }
    }
    //!extrapolate
    //...remove markers
    if(remove) M.erase(std::remove_if(M.begin(), M.end(), [](const MACgrid::Marker &m){return m.x<0.0f;}), M.end());
    //!remove
}

void Navier::updateSolids(int *S, const std::vector<MACgrid::Marker> &fluid) {
    for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
        if(S[s(x,y,z)]!=MACgrid::SOLID_CELL) S[s(x,y,z)] = MACgrid::AIR_CELL;
    }
    for(const auto &m : fluid) {
        int x = (int)m.x, y = (int)m.y, z = (int)m.z;
        S[s(x,y,z)] = MACgrid::FLUID_CELL;
    }
}

void Navier::advectVelocity(float *&U, const float dt, const float dx, float *&UU) {
    for(int x=1; x<Nx; x++) for(int y=1; y<Ny; y++) for(int z=1; z<Nz; z++) {
        {
            float a = x-U[u(x,y,z)]*dt/dx;
            float b = y-(U[v(x,y,z)]+U[v(x,y+1,z)]+U[v(x-1,y,z)]+U[v(x-1,y+1,z)])/4.0f*dt/dx;
            float c = z-(U[w(x,y,z)]+U[w(x,y,z+1)]+U[w(x-1,y,z)]+U[w(x-1,y,z+1)])/4.0f*dt/dx;

            int xx = (int)a, yy = (int)b, zz = (int)c;
            a -= xx; b -= yy; c -= zz;

            UU[u(x,y,z)] = (
                (1.0f-a)*(1.0f-b)*(1.0f-c) * U[u(xx,  yy,  zz  )] +
                (1.0f-a)*(1.0f-b)*(c)      * U[u(xx,  yy,  zz+1)] +
                (1.0f-a)*(b)*(1.0f-c)      * U[u(xx,  yy+1,zz  )] +
                (1.0f-a)*(b)*(c)           * U[u(xx,  yy+1,zz+1)] +
                (a)*(1.0f-b)*(1.0f-c)      * U[u(xx+1,yy,  zz  )] +
                (a)*(1.0f-b)*(c)           * U[u(xx+1,yy,  zz+1)] +
                (a)*(b)*(1.0f-c)           * U[u(xx+1,yy+1,zz  )] +
                (a)*(b)*(c)                * U[u(xx+1,yy+1,zz+1)]
            );
        }
        {
            float a = x-(U[u(x,y,z)]+U[u(x+1,y,z)]+U[u(x,y-1,z)]+U[u(x+1,y-1,z)])/4.0f*dt/dx;
            float b = y-U[v(x,y,z)]*dt/dx;
            float c = z-(U[w(x,y,z)]+U[w(x,y,z+1)]+U[w(x,y-1,z)]+U[w(x,y-1,z+1)])/4.0f*dt/dx;

            int xx = (int)a, yy = (int)b, zz = (int)c;
            a -= xx; b -= yy; c -= zz;

            UU[v(x,y,z)] = (
                (1.0f-a)*(1.0f-b)*(1.0f-c) * U[v(xx,  yy,  zz  )] +
                (1.0f-a)*(1.0f-b)*(c)      * U[v(xx,  yy,  zz+1)] +
                (1.0f-a)*(b)*(1.0f-c)      * U[v(xx,  yy+1,zz  )] +
                (1.0f-a)*(b)*(c)           * U[v(xx,  yy+1,zz+1)] +
                (a)*(1.0f-b)*(1.0f-c)      * U[v(xx+1,yy,  zz  )] +
                (a)*(1.0f-b)*(c)           * U[v(xx+1,yy,  zz+1)] +
                (a)*(b)*(1.0f-c)           * U[v(xx+1,yy+1,zz  )] +
                (a)*(b)*(c)                * U[v(xx+1,yy+1,zz+1)]
            );
        }
        {
            float a = x-(U[u(x,y,z)]+U[u(x+1,y,z)]+U[u(x,y,z-1)]+U[u(x+1,y,z-1)])/4.0f*dt/dx;
            float b = y-(U[v(x,y,z)]+U[v(x,y+1,z)]+U[v(x,y,z-1)]+U[v(x,y+1,z-1)])/4.0f*dt/dx;
            float c = z-U[w(x,y,z)]*dt/dx;

            int xx = (int)a, yy = (int)b, zz = (int)c;
            a -= xx; b -= yy; c -= zz;

            UU[w(x,y,z)] = (
                (1.0f-a)*(1.0f-b)*(1.0f-c) * U[w(xx,  yy,  zz  )] +
                (1.0f-a)*(1.0f-b)*(c)      * U[w(xx,  yy,  zz+1)] +
                (1.0f-a)*(b)*(1.0f-c)      * U[w(xx,  yy+1,zz  )] +
                (1.0f-a)*(b)*(c)           * U[w(xx,  yy+1,zz+1)] +
                (a)*(1.0f-b)*(1.0f-c)      * U[w(xx+1,yy,  zz  )] +
                (a)*(1.0f-b)*(c)           * U[w(xx+1,yy,  zz+1)] +
                (a)*(b)*(1.0f-c)           * U[w(xx+1,yy+1,zz  )] +
                (a)*(b)*(c)                * U[w(xx+1,yy+1,zz+1)]
            );
        }
    }
    float *tmp = U; U = UU; UU = tmp; //swap
}

void Navier::bodyForces(float *&U, const float *F, const int *S, const float dt, float *&UU) {
    for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
        UU[u(x,y,z)]   = U[u(x,y,z)];
        UU[u(x+1,y,z)] = U[u(x+1,y,z)];
        UU[v(x,y,z)]   = U[v(x,y,z)];
        UU[v(x,y+1,z)] = U[v(x,y+1,z)];
        UU[w(x,y,z)]   = U[w(x,y,z)];
        UU[w(x,y,z+1)] = U[w(x,y,z+1)];
    }
    for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
        if(S[s(x,y,z)]==MACgrid::FLUID_CELL) {
            UU[u(x,y,z)]   = U[u(x,y,z)]   + F[u(x,y,z)]*dt;
            UU[u(x+1,y,z)] = U[u(x+1,y,z)] + F[u(x+1,y,z)]*dt;
            UU[v(x,y,z)]   = U[v(x,y,z)]   + F[v(x,y,z)]*dt;
            UU[v(x,y+1,z)] = U[v(x,y+1,z)] + F[v(x,y+1,z)]*dt;
            UU[w(x,y,z)]   = U[w(x,y,z)]   + F[w(x,y,z)]*dt;
            UU[w(x,y,z+1)] = U[w(x,y,z+1)] + F[w(x,y,z+1)]*dt;
        }
    }
    float *tmp = U; U = UU; UU = tmp; //swap
}

float Navier::project(float *U, float *P, const int *S, const std::vector<MACgrid::Marker> &fluid, float *tmp, float *A) {
    float maxVelocity = 0.0f;

    //velocities on solids
    for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++)
        if(S[s(x,y,z)]==MACgrid::SOLID_CELL) {
            U[u(x,y,z)] = U[v(x,y,z)] = U[w(x,y,z)] = U[u(x+1,y,z)] = U[v(x,y+1,z)] = U[w(x,y,z+1)] = 0.0f;
        }

    {
        float *R  = tmp;
        float *W  = tmp + Nx*Ny*Nz;
        float *AW = tmp + 2*Nx*Ny*Nz;

        float rold = 0.0f, rnew, alpha;

        for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
            if(S[s(x,y,z)]==MACgrid::FLUID_CELL) {
                A[s(x,y,z)] = (6 -
                    (S[s(x-1,y,z)]==MACgrid::SOLID_CELL) -
                    (S[s(x+1,y,z)]==MACgrid::SOLID_CELL) -
                    (S[s(x,y-1,z)]==MACgrid::SOLID_CELL) -
                    (S[s(x,y+1,z)]==MACgrid::SOLID_CELL) -
                    (S[s(x,y,z-1)]==MACgrid::SOLID_CELL) -
                    (S[s(x,y,z+1)]==MACgrid::SOLID_CELL)
                );
                W[s(x,y,z)] = R[s(x,y,z)] = (U[u(x,y,z)]+U[v(x,y,z)]+U[w(x,y,z)]-U[u(x+1,y,z)]-U[v(x,y+1,z)]-U[w(x,y,z+1)]) -
                    A[s(x,y,z)]*P[s(x,y,z)] +
                    (S[s(x-1,y,z)]==MACgrid::FLUID_CELL)*P[s(x-1,y,z)] +
                    (S[s(x+1,y,z)]==MACgrid::FLUID_CELL)*P[s(x+1,y,z)] +
                    (S[s(x,y-1,z)]==MACgrid::FLUID_CELL)*P[s(x,y-1,z)] +
                    (S[s(x,y+1,z)]==MACgrid::FLUID_CELL)*P[s(x,y+1,z)] +
                    (S[s(x,y,z-1)]==MACgrid::FLUID_CELL)*P[s(x,y,z-1)] +
                    (S[s(x,y,z+1)]==MACgrid::FLUID_CELL)*P[s(x,y,z+1)]
                ;
                rold += R[s(x,y,z)]*R[s(x,y,z)];
            }
            else {
                A[s(x,y,z)] = 0.0f;
                W[s(x,y,z)] = R[s(x,y,z)] = 0.0f;
                AW[s(x,y,z)] = 0.0f;
                P[s(x,y,z)] = 0.0f;
            }
        }

        for(int iter=80; iter>0 && rold>1e-6; iter--) {
            alpha = 0.0f;
            for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
                if(S[s(x,y,z)]==MACgrid::FLUID_CELL) {
                    AW[s(x,y,z)] =
                        A[s(x,y,z)]*W[s(x,y,z)] -
                        (S[s(x-1,y,z)]==MACgrid::FLUID_CELL)*W[s(x-1,y,z)] -
                        (S[s(x+1,y,z)]==MACgrid::FLUID_CELL)*W[s(x+1,y,z)] -
                        (S[s(x,y-1,z)]==MACgrid::FLUID_CELL)*W[s(x,y-1,z)] -
                        (S[s(x,y+1,z)]==MACgrid::FLUID_CELL)*W[s(x,y+1,z)] -
                        (S[s(x,y,z-1)]==MACgrid::FLUID_CELL)*W[s(x,y,z-1)] -
                        (S[s(x,y,z+1)]==MACgrid::FLUID_CELL)*W[s(x,y,z+1)]
                    ;
                    alpha += W[s(x,y,z)]*AW[s(x,y,z)];
                }
            }
            alpha = rold/alpha;

            rnew = 0.0f;
            for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
                if(S[s(x,y,z)]==MACgrid::FLUID_CELL) {
                    P[s(x,y,z)] += alpha*W[s(x,y,z)];
                    R[s(x,y,z)] -= alpha*AW[s(x,y,z)];
                    rnew += R[s(x,y,z)]*R[s(x,y,z)];
                }
            }
            for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
                if(S[s(x,y,z)]==MACgrid::FLUID_CELL) {
                    W[s(x,y,z)] = R[s(x,y,z)] + (rnew/rold)*W[s(x,y,z)];
                }
            }
            rold = rnew;
        }
    }

    for(int x=1; x<Nx; x++) for(int y=1; y<Ny; y++) for(int z=1; z<Nz; z++) {
        U[u(x,y,z)] -= P[s(x,y,z)] - P[s(x-1,y,z)];
        U[v(x,y,z)] -= P[s(x,y,z)] - P[s(x,y-1,z)];
        U[w(x,y,z)] -= P[s(x,y,z)] - P[s(x,y,z-1)];
    }

    //velocities and pressure on solids
    for(int x=0; x<Nx; x++) for(int y=0; y<Ny; y++) for(int z=0; z<Nz; z++) {
        if(S[s(x,y,z)]==MACgrid::SOLID_CELL) {
            P[s(x,y,z)] = U[u(x,y,z)] + U[v(x,y,z)] + U[w(x,y,z)] - U[u(x+1,y,z)] - U[v(x,y+1,z)] - U[w(x,y,z+1)];
            U[u(x,y,z)] = U[v(x,y,z)] = U[w(x,y,z)] = U[u(x+1,y,z)] = U[v(x,y+1,z)] = U[w(x,y,z+1)] = 0.0f;
        }
        else {
            if(U[u(x,y,z)]>maxVelocity) maxVelocity=U[u(x,y,z)];
            if(U[v(x,y,z)]>maxVelocity) maxVelocity=U[v(x,y,z)];
            if(U[w(x,y,z)]>maxVelocity) maxVelocity=U[w(x,y,z)];
        }
    }

    return maxVelocity;
}
