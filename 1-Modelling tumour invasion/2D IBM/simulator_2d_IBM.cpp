#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>

using namespace std;

struct Cell{
    int i;
    int j;
    int age;
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
    double *n, *f, *m, *m_new;
    double P0, P1, P2, P3, P4;
    n = new double[nx*ny]();
    f = new double[nx*ny]();
    m = new double[nx*ny]();
    m_new = new double[nx*ny]();
    vector<Cell> cells;
    cells.reserve(3000);

    // For random numbers
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<int> ic(1, 200);
    std::uniform_real_distribution<double> loop(0, 1);
    std::uniform_int_distribution<int> proliferate(1,4);

    // Initial conditions
    while (cells.size() < 500){
        int i_rand = ic(gen);
        int j_rand = ic(gen);
        double x_rand =(i_rand+0.5)*hx;
        double y_rand =(j_rand+0.5)*hy;
        double r2=(x_rand-0.5)*(x_rand-0.5)+(y_rand-0.5)*(y_rand-0.5);
        if(r2<0.01){
            cells.push_back({i_rand, j_rand, 0});
            n[i_rand + j_rand*nx] += 1;
        }
    }
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            int ij= i + j*nx;
            n[ij]<2 ? f[ij]=1-0.5*n[ij] : f[ij] = 0;
            m[ij] = 0.5*n[ij];
        }
    }

    // Prepare binary file
    ofstream file("simulation_data.bin", ios::binary);

    // Simulation loop
    double t1=clock();
    for(int t_step=0; t_step<t; t_step++){

        //Update f and m
        for(int j=1; j<ny-1; j++){
            for(int i=1; i<nx-1; i++){
                int ij = i + j*nx;
                f[ij] = f[ij]*(1-k*eta*m[ij]);
                m_new[ij] = m[ij]*m_center_coeff
                        + m_periph_coeff*(m[ij+1]+m[ij-1]+m[ij+nx]+m[ij-nx]) + m_alpha_coeff*n[ij];
            }
        }
        swap(m, m_new);

        for(auto it = cells.begin(); it != cells.end(); ++it){
            int ij = it->i + it->j*nx;

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

            // Moving cells
            double sum=P0+P1+P2+P3+P4;
            double nb = loop(gen)*sum;
            if(nb < P0){/*stationary*/}
            else if(nb < P0+P1){
                // Move left
                if(it->i>1){it->i -= 1;}
            }
            else if(nb < P0+P1+P2){
                // Move right
                if(it->i<nx-2){it->i += 1;}
            }
            else if(nb < P0+P1+P2+P3){
                // Move up
                if(it->j>1){
                    it->j -= 1;
                }
            }
            else{
                // Move down
                if(it->j<ny-2){
                    it->j += 1;
                }
            }
        }

        // Update n
        fill(n, n+nx*ny, 0);
        for(auto it = cells.begin(); it != cells.end(); ++it){
            int ij = it->i + it->j*nx;
            n[ij] += 1;
        }

        // Proliferation?
        vector<Cell> divided_cells;
        for(auto it = cells.begin(); it != cells.end(); ++it){
            it->age += 1;

            if(it->age%500==0 && (cells.size()+divided_cells.size()) < 3000){
                int i = it->i;
                int j = it->j;
                int ij = i + j*nx;

                if(n[ij-1]==0 || n[ij+1]==0 || n[ij-nx]==0 || n[ij+nx]==0){
                    int size_t = divided_cells.size();
                    while(divided_cells.size() == size_t){
                        int invasion = proliferate(gen);

                        if(invasion == 1 && n[ij-1] == 0 && i>1){
                            divided_cells.push_back({i-1, j, 0});
                            n[ij-1] += 1;
                        }
                        else if(invasion == 2 && n[ij+1] == 0 && i<nx-2){
                            divided_cells.push_back({i+1, j, 0});
                            n[ij+1] += 1;
                        }
                        else if(invasion == 3 && n[ij-nx] == 0 && j<ny-2){
                            divided_cells.push_back({i, j-1, 0});
                            n[ij-nx] += 1;
                        }
                        else if(invasion == 4 && n[ij+nx] == 0 && j>1){
                            divided_cells.push_back({i, j+1, 0});
                            n[ij+nx] += 1;
                        }
                    }
                }
            }
        }

        cells.insert(cells.end(), divided_cells.begin(), divided_cells.end());

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
    delete [] m_new;

    return 0;
}