#include <iostream>
#include <cmath>
#include <fstream>
#include <ctime>

using namespace std;

struct Cell{
    int i;
    int j;
};

// Model parameters
const double dn=0.001, dm=0.001;
const double gamma=0.005, alpha=0.1;
const int eta=10, bet=0;
const double epsilon=0.0025;

// Time and space discretization
const int nx=202, ny=202; // Add 2 for boundary conditions
const double hx=1.0/(nx-2), hy=1.0/(ny-2), hx2=hx*hx, hy2=hy*hy;
const int t_final=15;
const double k=0.001;
const int t=t_final/k;
const int saving_step=100;

// Gaining calculation time
const double m_center_coeff = 1 -4*(k*dm/hx2) -bet*k;
const double m_periph_coeff = k*dm/hx2;
const double m_alpha_coeff = alpha*k;
const double P0_coeff1 = 1-(4*dn*k/hx2);
const double P0_coeff2 = k*gamma/hx2;
const double P_coeff1 = dn*k/hx2;
const double P_coeff2 = k*gamma/(4*hx2);

int main(){
    // Model variables
    double *n, *f, *m, *n_new, *m_new;
    double P0, P1, P2, P3, P4;
    n = new double[nx*ny]();
    f = new double[nx*ny]();
    m = new double[nx*ny]();
    n_new = new double[nx*ny]();
    m_new = new double[nx*ny]();

    // Initial conditions
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            int ij = i + j*nx;
            double x=(i+0.5)*hx;
            double y=(j+0.5)*hy;
            double r2=(x-0.5)*(x-0.5)+(y-0.5)*(y-0.5);
            if(r2<0.01){
                n[ij] = exp(-r2/epsilon);
                f[ij] = 1-0.5*n[ij];
                m[ij] = 0.5*n[ij];
            }
            else{
                n[ij] = 0;
                f[ij] = 1;
                m[ij] = 0;
            }
            n_new[ij]=0;
            m_new[ij]=0;
        }
    }

    // Prepare binary file
    ofstream file("simulation_data_2d.bin", ios::binary);

    // Simulation loop
    double t1=clock();
    for(int t_step=0; t_step<t; t_step++){

        //Update f and m
        for(int j=1; j<ny-1; j++){
            for(int i=1; i<nx-1; i++){
                int ij = i + j*nx;

                P0 = P0_coeff1 - P0_coeff2*(f[ij-1]+f[ij+1]-4*f[ij]+f[ij-nx]+f[ij+nx]);
                P1 = P_coeff1 -P_coeff2*(f[ij+1]-f[ij-1]);
                P2 = P_coeff1 +P_coeff2*(f[ij+1]-f[ij-1]);
                P3 = P_coeff1 -P_coeff2*(f[ij+nx]-f[ij-nx]);
                P4 = P_coeff1 +P_coeff2*(f[ij+nx]-f[ij-nx]);

                if(P0<0) P0=0;
                if(P1<0) P1=0;
                if(P2<0) P2=0;
                if(P3<0) P3=0;
                if(P4<0) P4=0;

                n_new[ij] = n[ij]*P0 + n[ij+1]*P1 + n[ij-1]*P2 + n[ij+nx]*P3 + n[ij-nx]*P4;
                f[ij] = f[ij]*(1-k*eta*m[ij]);
                m_new[ij] = m[ij]*m_center_coeff
                        + m_periph_coeff*(m[ij+1]+m[ij-1]+m[ij+nx]+m[ij-nx]) + m_alpha_coeff*n[ij];
            }
        }
        swap(n,n_new);
        swap(m, m_new);

        // Saving
        if(t_step%saving_step==0){
            file.write(reinterpret_cast<char*>(n), nx * ny * sizeof(double));
            file.write(reinterpret_cast<char*>(f), nx * ny * sizeof(double));
            file.write(reinterpret_cast<char*>(m), nx * ny * sizeof(double));
        }
    }

    double t2=clock();
    double temps=(t2-t1)/CLOCKS_PER_SEC;
    std::cout << "Temps ecoule: " << temps << " secondes." << std::endl;

    delete [] n;
    delete [] f;
    delete [] m;
    delete [] n_new;
    delete [] m_new;

    return 0;
}