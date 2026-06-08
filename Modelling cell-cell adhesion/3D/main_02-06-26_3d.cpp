#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <sstream>
#include <iomanip>
#include <string>
using namespace std;

struct cell {
    double x,y,z;
    double force_x, force_y,force_z;
    double vx, vy, vz, vx_new, vy_new, vz_new;
    double dt_Ec, dt_Eb, dt_b, dt_C;
    double Ec, Eb, b, C;
    bool is_boundary;
    bool is_invasive;
    bool is_recovering;
    double t_invasive;
};

// Time
const double t_max=50;
const double dt=0.0001;
const int t_save=0.5/dt;

// Viscosity
const double gam=0.4;//0.4;
const double gam_par=1e3, gam_perp=1e3;

// Other parameters
const double Et=100;
const double Pt=0.33514;
const double nu=100, kminus=19, km=0.01, k2=0.03, E=1;
const double cT=50,alpha=2;
const int nx=30, ny=25;
const double R=5;
const double d_eq=8.0;
const double R_eq=d_eq/2;
const int n_cells=nx*ny;

int main(){
    double *area = new double[n_cells*n_cells]();

    vector<cell> cells;
    cells.reserve(nx*ny);
    vector<vector<int>> neigh(n_cells);

    random_device rd{};
    mt19937 gen(rd());
    normal_distribution<double> dist(0.0, 1.0);
    double amplitude_bruit = 2e-4;

    // Initialisation
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            double x=i*2*R_eq+(j%2)*R_eq;
            double y=sqrt(3)*R_eq*j;
            double z=0;
            bool boundary=(i==0||i==nx-1||j==0||j==ny-1);
            cells.push_back({x,y,z,0,0,0,0,0,0,0,0,0,0,0,0,0,0,100,(km/k2*(kminus+k2))/(100*(Pt-km/k2)),km/k2,boundary, false, false,0});
        }
    }
    for(int i=0; i<n_cells; i++){
        for(int j=i+1; j<n_cells; j++){
            double dx=cells[j].x-cells[i].x;
            double dy=cells[j].y-cells[i].y;
            double dz=cells[j].z-cells[i].z;
            double d=sqrt(dx*dx+dy*dy+dz*dz);
            if(d<1.1*d_eq && i!=j){
                neigh[i].push_back(j);
                neigh[j].push_back(i);
            }
        }
    }

    int frame = 0;

    // Main loop
    double t1=omp_get_wtime();
    for(double t=0; t<=t_max/dt; t++){

        // 1. On pré-calcule le bruit séquentiellement (Thread-Safe)
        vector<double> bruit_x(n_cells), bruit_y(n_cells), bruit_z(n_cells);
        for(int i=0; i<n_cells; i++){
            bruit_x[i] = dist(gen) * amplitude_bruit;
            bruit_y[i] = dist(gen) * amplitude_bruit;
            bruit_z[i] = dist(gen) * amplitude_bruit;
        }

        #pragma omp parallel for
        for(int idx=0; idx<n_cells; idx++){

            if(cells[idx].is_boundary)continue;

            double c_i=0, d_i=0;
            double M_xx = gam;
            double M_yy = gam;
            double M_zz = gam;
            double M_xz = 0;
            double M_xy = 0;
            double M_yz = 0;
            double total_contact_area = 0;
            cells[idx].force_x=0; cells[idx].force_y=0; cells[idx].force_z=0;

            for(int k=0; k<neigh[idx].size(); k++){
                int jdx=neigh[idx][k];
                int ij=idx+jdx*n_cells;
                double dx=cells[jdx].x-cells[idx].x;
                double dy=cells[jdx].y-cells[idx].y;
                double dz=cells[jdx].z-cells[idx].z;
                double d=sqrt(dx*dx + dy*dy+dz*dz);
                if(d < 0.001) d = 0.001;
                double new_area=0;

                if(int(t)%t_save==0 && idx==nx*ny/2+1) cout << "d: "<< d <<" | ";

                // Biophysics
                if(d<2*R){ // Cells are touching

                    // Area
                    new_area=(2*R-d)/(4*R);
                    total_contact_area+=new_area;

                    // Repulsion
                    double F_rep=E*(10-d)*sqrt(10-d)*sqrt(2.5)*2./3.;
                    // Adhesion
                    double F_adh=2.981424*min(cells[idx].Eb/Et, cells[jdx].Eb/Et);

                    cells[idx].force_x+=(-F_rep+F_adh)*dx/d;
                    cells[idx].force_y+=(-F_rep+F_adh)*dy/d;
                    cells[idx].force_z+=(-F_rep+F_adh)*dz/d;

                    // Friction matrix
                    double dx_n = dx/d;
                    double dy_n = dy/d;
                    double dz_n = dz/d;

                    // Astuce mathématique : M_xy = (gam_par - gam_perp) * dx_n * dy_n
                    double diff = gam_par - gam_perp;

                    double M_local_xx = gam_par*dx_n*dx_n + gam_perp*(1.0 - dx_n*dx_n);
                    double M_local_yy = gam_par*dy_n*dy_n + gam_perp*(1.0 - dy_n*dy_n);
                    double M_local_zz = gam_par*dz_n*dz_n + gam_perp*(1.0 - dz_n*dz_n);

                    double M_local_xy = diff * dx_n * dy_n;
                    double M_local_xz = diff * dx_n * dz_n;
                    double M_local_yz = diff * dy_n * dz_n;

                    // Ajout à la matrice globale de la cellule
                    M_xx += M_local_xx;
                    M_yy += M_local_yy;
                    M_zz += M_local_zz;
                    M_xy += M_local_xy;
                    M_xz += M_local_xz;
                    M_yz += M_local_yz;

                    // Ajout de la force de friction dynamique due à la vitesse du voisin
                    cells[idx].force_x += M_local_xx * cells[jdx].vx + M_local_xy * cells[jdx].vy + M_local_xz * cells[jdx].vz;
                    cells[idx].force_y += M_local_xy * cells[jdx].vx + M_local_yy * cells[jdx].vy + M_local_yz * cells[jdx].vz;
                    cells[idx].force_z += M_local_xz * cells[jdx].vx + M_local_yz * cells[jdx].vy + M_local_zz * cells[jdx].vz;
                }

                // Biochemistry
                if(t>0){
                    double da=(new_area-area[ij])/dt;
                    if(da>0)c_i+=da*200;
                    else if(da<0)d_i+=abs(da)*200;
                }
                area[ij]=new_area;
            }

            if(int(t)%t_save==0 && idx==nx*ny/2+1) cout << "area: "<< total_contact_area <<" | ";

            if (c_i * dt > 1.0) c_i = 1.0 / dt;
            if (d_i * dt > 1.0) d_i = 1.0 / dt;

            // Inversion de la matrice de friction 3x3 (Règle de Cramer)
            double det = M_xx * (M_yy * M_zz - M_yz * M_yz)
                       - M_xy * (M_xy * M_zz - M_xz * M_yz)
                       + M_xz * (M_xy * M_yz - M_xz * M_yy);

            M_xx+=1;
            M_yy+=1;
            M_zz+=1;

            // Calcul des éléments de la matrice inverse
            double I_xx = (M_yy * M_zz - M_yz * M_yz) / det;
            double I_yy = (M_xx * M_zz - M_xz * M_xz) / det;
            double I_zz = (M_xx * M_yy - M_xy * M_xy) / det;
            double I_xy = (M_xz * M_yz - M_xy * M_zz) / det;
            double I_xz = (M_xy * M_yz - M_xz * M_yy) / det;
            double I_yz = (M_xy * M_xz - M_xx * M_yz) / det;

            if(cells[idx].is_invasive) {
                // 1. La force chimiotactique (À appliquer sur X ou Y, pas Z)
                // Commence avec une valeur modérée (par exemple 10.0 ou 20.0)
                // et ajuste-la. Elle doit rivaliser avec tes forces de répulsion (F_rep).
                double F_chemo = 500;
                cells[idx].force_z += F_chemo;
            }

            // Calcul des nouvelles vitesses (v = M_inverse * Force)
            cells[idx].vx_new = I_xx * cells[idx].force_x + I_xy * cells[idx].force_y + I_xz * cells[idx].force_z;
            cells[idx].vy_new = I_xy * cells[idx].force_x + I_yy * cells[idx].force_y + I_yz * cells[idx].force_z;
            cells[idx].vz_new = I_xz * cells[idx].force_x + I_yz * cells[idx].force_y + I_zz * cells[idx].force_z;

            // if(cells[idx].is_invasive /*&& cells[idx].Eb <5*/) {
            //     cells[idx].vz_new += 1e-1;
            //     bruit_x[idx]*=10;
            //     bruit_y[idx]*=10;
            //     bruit_z[idx]*=10;
            // }

            if(cells[idx].is_invasive) {
                bruit_z[idx] *= 10.0;
            }

            cells[idx].vx_new += bruit_x[idx];
            cells[idx].vy_new += bruit_y[idx];
            cells[idx].vz_new += bruit_z[idx];


            // We force beta > CT
            if(t==int(1/dt) && idx==nx*ny/2-5){
                cells[idx].b=100-(km/k2*(kminus+k2))/(100*(Pt-km/k2));
                cells[idx].is_invasive=true;
                cells[idx].t_invasive = t * dt;
            }

            if(cells[idx].b>cT && !cells[idx].is_recovering && !cells[idx].is_invasive){
                cells[idx].is_invasive=true;
                cells[idx].t_invasive = t * dt;
            }

            if (cells[idx].is_invasive && ((t * dt) - cells[idx].t_invasive) >= 5.0){
                cells[idx].is_invasive=false;
                cells[idx].is_recovering=true;

                for(int k=0; k<neigh[idx].size(); k++){
                    int jdx=neigh[idx][k];
                    int ij=idx+jdx*n_cells;
                    area[ij]=0;
                }
            }

            if(int(t)%t_save==0 && idx==nx*ny/2+1) cout << "ci: "<< c_i<< " | di: "<< d_i << " | b: "<< cells[idx].b << " | Eb: "<< cells[idx].Eb
            << " | C: "<< cells[idx].C << " | Ec: "<< cells[idx].Ec <<  endl;

            double kplus = 100;
            if (cells[idx].is_invasive) {
                kplus=0;
                cells[idx].dt_Eb = -(alpha + d_i) * cells[idx].Eb;
                cells[idx].dt_b  =  (alpha + d_i) * cells[idx].Eb - kplus * cells[idx].b * (Pt - cells[idx].C) + kminus * cells[idx].C + km;
                cells[idx].dt_Ec = -c_i * cells[idx].Ec + (alpha + d_i) * cells[idx].Eb;
            } else {
                // --- SÉCURITÉ OBLIGATOIRE ---
                double espace_libre = Et - cells[idx].Ec - cells[idx].Eb;
                if (espace_libre < 0.0) espace_libre = 0.0;

                double prod = nu * espace_libre * cells[idx].b;
                if (prod * dt > espace_libre) {
                    prod = espace_libre / dt;
                }

                cells[idx].dt_Eb =  prod - d_i * cells[idx].Eb;
                cells[idx].dt_b  = -prod + d_i * cells[idx].Eb - kplus * cells[idx].b * (Pt - cells[idx].C) + kminus * cells[idx].C + km;
                cells[idx].dt_Ec = -c_i * cells[idx].Ec + d_i * cells[idx].Eb;
            }
            cells[idx].dt_C  = kplus * cells[idx].b * (Pt - cells[idx].C) - kminus * cells[idx].C - k2 * cells[idx].C;
        }

        #pragma omp parallel for
        for(int idx=0; idx<n_cells; idx++){
            if(cells[idx].is_boundary)continue;

            // Update ODEs
            cells[idx].Ec+=dt*cells[idx].dt_Ec;
            cells[idx].Eb+=dt*cells[idx].dt_Eb;
            cells[idx].b+=dt*cells[idx].dt_b;
            cells[idx].C+=dt*cells[idx].dt_C;

            // Update velocity
            cells[idx].vx = cells[idx].vx_new;
            cells[idx].vy = cells[idx].vy_new;
            cells[idx].vz = cells[idx].vz_new;

            // Update position
            cells[idx].x += cells[idx].vx*dt;
            cells[idx].y += cells[idx].vy*dt;
            cells[idx].z += cells[idx].vz*dt;
        }

        // Save en format VTK pour ParaView
        if(int(t)%t_save==0 && t!=0){
            double t_actuel = t * dt;
            cout << " (t = " << t_actuel << ")" << endl;

            // Formatage du nom de fichier avec des zéros (ex: cells_0000.vtk)
            std::ostringstream filename;
            filename << "cells_" << std::setw(4) << std::setfill('0') << frame << ".vtk";
            std::ofstream vtk_file(filename.str());

            // --- En-tête VTK ---
            vtk_file << "# vtk DataFile Version 3.0\n";
            vtk_file << "Simulation cellules temps " << t_actuel << "\n";
            vtk_file << "ASCII\n";
            vtk_file << "DATASET POLYDATA\n";

            // --- Coordonnées des points (Cellules) ---
            vtk_file << "POINTS " << n_cells << " float\n";
            for(int idx = 0; idx < n_cells; idx++){
                // On ajoute un Z = 0.0 car ParaView travaille en 3D
                vtk_file << cells[idx].x << " " << cells[idx].y << " " << cells[idx].z << "\n";
            }

            // --- Valeurs attachées aux cellules (la beta-catenine) ---
            vtk_file << "POINT_DATA " << n_cells << "\n";
            vtk_file << "SCALARS b_adim float 1\n";
            vtk_file << "LOOKUP_TABLE default\n";
            for(int idx = 0; idx < n_cells; idx++){
                double b_adim = cells[idx].b / Et;
                vtk_file << b_adim << "\n";
            }

            vtk_file.close();
            frame++;
        }
    }
    double t2=omp_get_wtime();
    double temps=t2-t1;
    std::cout << "Fin de simulation | Temps : " << temps << " secondes." << std::endl;

    delete[] area;
    return 0;
}