#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
using namespace std;

struct cell {
    double x,y;
    double force_x, force_y;
    double vx, vy, vx_new, vy_new;
    double dt_Ec, dt_Eb, dt_b, dt_C;
    double Ec, Eb, b, C;
    bool is_boundary;
    bool is_invasive;
    bool is_recovering;
    double t_invasive;
};

// Time
const double t_max=40;
const double dt=0.0001;
const int t_save=t_max/(dt*40);

// Viscosity
const double gam=0.4;
const double gam_par=1e3, gam_perp=1e3;

// Other parameters
const double Et=100;
const double Pt=0.33514;
const double nu=100, kminus=19, km=0.01, k2=0.03, E=1;
const double cT=50,alpha=2;
const int nx=30, ny=30;
const double R=5;
const double d_eq=8.0;
const double R_eq=d_eq/2;
const int n_cells=nx*ny;

int main(){
    double *area = new double[n_cells*n_cells]();

    vector<cell> cells;
    cells.reserve(nx*ny);
    vector<vector<int>> neigh(n_cells);

    // Initialisation
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            double x=i*2*R_eq+(j%2)*R_eq;
            double y=sqrt(3)*R_eq*j;
            bool boundary=(i==0||i==nx-1||j==0||j==ny-1);
            cells.push_back({x,y,0,0,0,0,0,0,0,0,0,0,0,100,(km/k2*(kminus+k2))/(100*(Pt-km/k2)),km/k2,boundary, false, false,0});
        }
    }
    for(int i=0; i<n_cells; i++){
        for(int j=i+1; j<n_cells; j++){
            double dx=cells[j].x-cells[i].x;
            double dy=cells[j].y-cells[i].y;
            double d=sqrt(dx*dx+dy*dy);
            if(d<1.1*d_eq && i!=j){
                neigh[i].push_back(j);
                neigh[j].push_back(i);
            }
        }
    }

    ofstream data("data.bin", ios::binary);

    // Main loop
    double t1=omp_get_wtime();
    for(double t=0; t<=t_max/dt; t++){

        #pragma omp parallel for
        for(int idx=0; idx<n_cells; idx++){

            double c_i=0, d_i=0;
            double M_xx = gam;
            double M_yy = gam;
            double M_xy = 0;
            double total_contact_area = 0;
            cells[idx].force_x=0; cells[idx].force_y=0;

            for(int k=0; k<neigh[idx].size(); k++){
                int jdx=neigh[idx][k];
                int ij=idx+jdx*n_cells;
                double dx=cells[jdx].x-cells[idx].x;
                double dy=cells[jdx].y-cells[idx].y;
                double d=sqrt(dx*dx + dy*dy);
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

                    if(!cells[idx].is_boundary){
                        cells[idx].force_x+=(-F_rep+F_adh)*dx/d;
                        cells[idx].force_y+=(-F_rep+F_adh)*dy/d;
                    }

                    // Friction matrix
                    double dx_n = dx/d;
                    double dy_n = dy/d;
                    double M_par_xx = gam_par*dx_n*dx_n;
                    double M_par_xy = gam_par*dx_n*dy_n;
                    double M_par_yy = gam_par*dy_n*dy_n;
                    double M_perp_xx = gam_perp*(1-dx_n*dx_n);
                    double M_perp_xy = gam_perp*(-dx_n*dy_n);
                    double M_perp_yy = gam_perp*(1-dy_n*dy_n);

                    M_xx += M_par_xx+M_perp_xx;
                    M_yy += M_par_yy+M_perp_yy;
                    M_xy += M_par_xy+M_perp_xy;

                    if(!cells[idx].is_boundary){
                        cells[idx].force_x += (M_par_xx+M_perp_xx)*cells[jdx].vx+(M_par_xy+M_perp_xy)*cells[jdx].vy;
                        cells[idx].force_y += (M_par_xy+M_perp_xy)*cells[jdx].vx+(M_par_yy+M_perp_yy)*cells[jdx].vy;
                    }
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
            if(c_i>1/dt)c_i = 1/dt;

            double det=M_xx*M_yy-M_xy*M_xy;
            if(!cells[idx].is_boundary){
                cells[idx].vx_new=(M_yy*cells[idx].force_x-M_xy*cells[idx].force_y)/det;
                cells[idx].vy_new=(M_xx*cells[idx].force_y-M_xy*cells[idx].force_x)/det;
            }

            // We force beta > CT
            if(t==int(1/dt) && idx==nx*ny/2+14){
                cells[idx].b=100-(km/k2*(kminus+k2))/(100*(Pt-km/k2));
                cells[idx].is_invasive=true;
                cells[idx].t_invasive = t*dt;
            }

            if(cells[idx].b>cT && !cells[idx].is_recovering && !cells[idx].is_invasive){
                cells[idx].is_invasive=true;
                cells[idx].t_invasive = t*dt;
            }

            if (cells[idx].is_invasive && ((t*dt)-cells[idx].t_invasive)>=5.0){
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

            double kplus=100;
            if (cells[idx].is_invasive) {
                kplus = 0;
                cells[idx].dt_Eb=-(alpha+d_i)*cells[idx].Eb;
                cells[idx].dt_b =(alpha+d_i)*cells[idx].Eb- kplus*cells[idx].b*(Pt-cells[idx].C) +kminus*cells[idx].C +km;
            } else {
                cells[idx].dt_Eb =nu*(Et-cells[idx].Ec-cells[idx].Eb)*cells[idx].b -d_i*cells[idx].Eb;
                cells[idx].dt_b =-nu*(Et-cells[idx].Ec-cells[idx].Eb)*cells[idx].b +d_i*cells[idx].Eb -kplus*cells[idx].b*(Pt-cells[idx].C) + kminus*cells[idx].C +km;
            }
            cells[idx].dt_Ec= -c_i*cells[idx].Ec +d_i*cells[idx].Eb;
            cells[idx].dt_C =kplus*cells[idx].b*(Pt-cells[idx].C) -kminus*cells[idx].C -k2*cells[idx].C;
        }

        #pragma omp parallel for
        for(int idx=0; idx<n_cells; idx++){

            // Update ODEs
            cells[idx].Ec+=dt*cells[idx].dt_Ec;
            cells[idx].Eb+=dt*cells[idx].dt_Eb;
            cells[idx].b+=dt*cells[idx].dt_b;
            cells[idx].C+=dt*cells[idx].dt_C;

            // Update velocity
            if(!cells[idx].is_boundary){
                cells[idx].vx = cells[idx].vx_new;
                cells[idx].vy = cells[idx].vy_new;

                // Update position
                cells[idx].x += cells[idx].vx*dt;
                cells[idx].y += cells[idx].vy*dt;
            }
        }

        // Save
        if((int(t)%t_save==0 && t!=0) || int(t)==1/dt){
            double t_actuel=t*dt;
            cout << t_actuel <<endl;
            for(int idx=0; idx<n_cells-150; idx++){
                double b_adim=cells[idx].b/Et;
                data.write(reinterpret_cast<const char*>(&t_actuel), sizeof(double));
                data.write(reinterpret_cast<const char*>(&cells[idx].x), sizeof(double));
                data.write(reinterpret_cast<const char*>(&cells[idx].y), sizeof(double));
                data.write(reinterpret_cast<const char*>(&b_adim), sizeof(double));
            }
        }
    }
    double t2=omp_get_wtime();
    double temps=t2-t1;
    std::cout << "Fin de simulation | Temps : " << temps << " secondes." << std::endl;

    data.close();

    system("python -m plotter");

    delete[] area;
    return 0;
}