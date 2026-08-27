#include <iostream>
#include <cmath>
#include <random>
#include <fstream>
#include <ctime>
#include <vector>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <omp.h>

using namespace std;

struct Cell{
    int i;
    int j;
    int age;
    int id;
    int type; // 1, 2, 3, 4, 5=dead
};

struct Vessel{
    int i;
    int j;
    double age;
    bool active;
    int old_i;
    int old_j;
};

// Model parameters
const double dn=0.0005, dm=0.0005, dc=0.005, dcell=0.0025, dv=1e-5;
const double rho=0.001, nu=0, omega=0.57, phi=0.025;
const int eta=1, kappa=1, sigma=2;
const double pmut=0.1;
const double omega1=omega, omega2=4*omega/3, omega3=2*omega, omega4=4*omega, omega5=omega/2;
const double kappa1=kappa, kappa2=4*kappa/3, kappa3=2*kappa, kappa4=4*kappa;
const double rho1=rho, rho2=4*rho/3, rho3=2*rho, rho4=4*rho;
const double chi=0.6, rho5 = 0.34, alpha=1;

// Time and space discretization
const int nx=402, ny=402; // Add 2 for boundary conditions
const double hx=1.0/(nx-2), hy=1.0/(ny-2), hx2=hx*hx, hy2=hy*hy;
const int t_final=20;
const double k=0.0005;
const int t=t_final/k;
const int saving_step=2000; // t=1,2,3,...

// Gaining calculation time
const double m_center_coeff = 1 -4*(k*dm/hx2) -sigma*k;
const double m_periph_coeff = k*dm/hx2;
const double m_kappa_coeff = kappa*k;
const double P0_coeff1 = 1-(4*dn*k/hx2);
const double P0_coeff2 = k*rho/hx2;
const double P_coeff1 = dn*k/hx2;
const double P_coeff2 = k*rho/(4*hx2);
const double it_type1 = (16.0/16.0)*(1.0/k);
const double it_type2 = (14.0/16.0)*(1.0/k);
const double it_type3 = (12.0/16.0)*(1.0/k);
const double it_type4 = (8.0/16.0)*(1.0/k);
const double P0_coeff2_rho1 = k*rho1/hx2;
const double P_coeff2_rho1 = k*rho1/(4*hx2);
const double P0_coeff2_rho2 = k*rho2/hx2;
const double P_coeff2_rho2 = k*rho2/(4*hx2);
const double P0_coeff2_rho3 = k*rho3/hx2;
const double P_coeff2_rho3 = k*rho3/(4*hx2);
const double P0_coeff2_rho4 = k*rho4/hx2;
const double P_coeff2_rho4 = k*rho4/(4*hx2);
const double c_center_empty = 1.0 + (4.0*dc*k)/hx2 + k*phi;
const double c_center_cell  = 1.0 + (4.0*dcell*k)/hx2+k*phi;
const double c_periph_empty = dc*k/hx2;
const double c_periph_cell  = dcell*k/hx2;

const double tol=5e-5;
const int max_iter=10;

const double dTAF=0.01;
const double alpha_TAF=10;//production TAF
const double lambda_TAF=0.001;
const double p_angio=1e2;
const double a_center=1+(4.*dTAF*k)/hx2+k*lambda_TAF;
const double a_periph=dTAF*k/hx2;

int main(){
    // Model variables
    double *n, *f, *f_new, *m, *m_new, *c, *c_new, *omega_conso, *a, *a_new, *v;
    double P0, P1, P2, P3, P4;
    n = new double[nx*ny]();
    f = new double[nx*ny]();
    f_new = new double[nx*ny]();
    m = new double[nx*ny]();
    m_new = new double[nx*ny]();
    c = new double[nx*ny]();
    c_new = new double[nx*ny]();
    omega_conso = new double[nx*ny]();
    a = new double[nx*ny]();
    a_new = new double[nx*ny]();
    v = new double[nx*ny]();

    vector<Cell> cells;
    cells.reserve(100000);
    vector<Cell> new_cells;
    new_cells.reserve(100000);
    vector<Cell> dead_cells;
    dead_cells.reserve(100000);
    vector<Vessel> vessels;
    vessels.reserve(10000);

    // For random numbers
    random_device rd;
    mt19937 gen(rd());
    std::uniform_int_distribution<int> ic_x(315, 325);
    std::uniform_int_distribution<int> ic_y(195, 205);
    std::uniform_real_distribution<double> loop(0, 1);
    std::uniform_int_distribution<int> proliferate(1,4);
    std::uniform_int_distribution<int> age_rand(0,2000);

    // Initial conditions
    int id_nb = 1;
    while (cells.size() < 50){
        int i_rand = ic_x(gen);
        int j_rand = ic_y(gen);
        int ij= i_rand + j_rand*nx;

        double x_rand =(i_rand+0.5)*hx;
        double y_rand =(j_rand+0.5)*hy;
        double r2=(x_rand-0.8)*(x_rand-0.8)+(y_rand-0.5)*(y_rand-0.5);

        if(r2<0.0003 && n[ij]==0){
            int new_age=age_rand(gen);
            cells.push_back({i_rand, j_rand, new_age, id_nb, 1});
            n[ij] += 1;
            id_nb += 1;
        }
    }
    for(int j=0; j<ny; j++){
        for(int i=0; i<nx; i++){
            int ij= i + j*nx;
            f[ij] = 1;
            a[ij]=0;
            if(i==1){v[ij]=1;c[ij]=1;}
            else{v[ij]=0;c[ij]=1.- i/400.;}
        }
    }

    for(int j=0; j<nx; j++){
        if(j==68 || j==120 || j==200 || j==260 || j==336){
            vessels.push_back({2,j,0,true,1,j});
            int ij=2+j*nx;
            v[ij]=1;
            c[ij]=1e-3;
        }
    }

    // Prepare binary file
    ofstream continuous("continuous_data.bin", ios::binary);
    ofstream discrete("discrete_data.txt");

    // Simulation loop
    cout << "Debut de la simulation" << endl;
    double t1=clock();
    for(int t_step=0; t_step<t+1; t_step++){

        //Update f and m and prepare c
        for(int j=1; j<ny-1; j++){
            for(int i=1; i<nx-1; i++){
                int ij = i + j*nx;
                f_new[ij]=f[ij]+k*(0.05*v[ij] -0.1*v[ij]*f[ij] -eta*m[ij]*f[ij]);
                if(f[ij]<0)cout<<f[ij]<<endl;
                m_new[ij] = m[ij]*m_center_coeff + m_periph_coeff*(m[ij+1]+m[ij-1]+m[ij+nx]+m[ij-nx]);
                c_new[ij]=c[ij];
                omega_conso[ij]=0;
                // Boundary
                if(i==1){m_new[ij-1] = m_new[ij]; f[ij-1] = f[ij]; c_new[ij-1]=c_new[ij]; omega_conso[ij-1]= omega_conso[ij];}
                if(i==nx-2){m_new[ij+1] = m_new[ij]; f[ij+1] = f[ij]; c_new[ij+1]=c_new[ij]; omega_conso[ij+1]= omega_conso[ij];}
                if(j==1){m_new[ij-nx] = m_new[ij]; f[ij-nx] = f[ij]; c_new[ij-nx]=c_new[ij]; omega_conso[ij-nx]= omega_conso[ij];}
                if(j==ny-2){m_new[ij+nx] = m_new[ij]; f[ij+nx] = f[ij]; c_new[ij+nx]=c_new[ij]; omega_conso[ij+nx]= omega_conso[ij];}
            }
        }
        swap(f,f_new);
        for(auto it = cells.begin(); it != cells.end(); ++it){
            int i = it->i;
            int j = it->j;
            int ij = i+j*nx;

            double kappa_type=0, omega_type=0;
            if(it->type==1){kappa_type = kappa1; omega_type = omega1;}
            else if(it->type==2){kappa_type = kappa2; omega_type = omega2;}
            else if(it->type==3){kappa_type = kappa3; omega_type = omega3;}
            else if(it->type==4){kappa_type = kappa4; omega_type = omega4;}

            int Ae = n[ij-nx]+ n[ij+1]+n[ij+nx] +n[ij-1];
            if(Ae == 4){omega_type = omega5;}
            m_new[ij] += kappa_type*k*n[ij];
            omega_conso[ij] = omega_type;

            // Boundary
            if(i==1){m_new[ij-1] = m_new[ij]; omega_conso[ij-1]= omega_conso[ij];}
            if(i==nx-2){m_new[ij+1] = m_new[ij]; omega_conso[ij+1]= omega_conso[ij];}
            if(j==1){m_new[ij-nx] = m_new[ij]; omega_conso[ij-nx]= omega_conso[ij];}
            if(j==ny-2){m_new[ij+nx] = m_new[ij]; omega_conso[ij+nx]= omega_conso[ij];}
        }
        swap(m, m_new);

        // Gauss-Seidel
        for(int iter=0; iter<max_iter; iter++){
            double diff_max=0.0;

            // Red
            #pragma omp parallel for reduction(max:diff_max)
            for(int j=1; j<ny-1; j++){
                for(int i=1; i<nx-1; i++){
                    if((i+j)%2==0){
                        int ij = i + j*nx;
                        if(v[ij]==1){c_new[ij]=1;}
                        else{
                            double old=c_new[ij];
                            double c_periph, c_center;
                            c_periph = (n[ij] == 0) ? c_periph_empty : c_periph_cell;
                            c_center = (n[ij] == 0) ? c_center_empty : c_center_cell;
                            c_new[ij]= (c[ij] + c_periph*(c_new[ij-1]+c_new[ij+1]+c_new[ij-nx]+c_new[ij+nx])
                                +k*nu*f[ij]- k*omega_conso[ij]*n[ij])/c_center;
                            double diff=fabs(c_new[ij]-old);
                            if(diff>diff_max)diff_max=diff;
                        }
                        double old_a=a_new[ij];
                        double TAF_prod=0;
                        if(n[ij]>0 && c_new[ij]<0.1){
                            TAF_prod=k*alpha_TAF;
                        }
                        double a_center2=a_center+ k*0.1*v[ij];
                        a_new[ij]=(a[ij]+a_periph*(a_new[ij-1]+a_new[ij+1]+a_new[ij-nx]+a_new[ij+nx])+TAF_prod)/a_center2;
                    }

                }
            }

            // Black
            #pragma omp parallel for reduction(max:diff_max)
            for(int j=1; j<ny-1; j++){
                for(int i=1; i<nx-1; i++){
                    if((i+j)%2!=0){
                        int ij = i + j*nx;
                        if(v[ij]==1){c_new[ij]=1;}
                        else{
                            double old=c_new[ij];
                            double c_periph, c_center;
                            c_periph = (n[ij] == 0) ? c_periph_empty : c_periph_cell;
                            c_center = (n[ij] == 0) ? c_center_empty : c_center_cell;
                            c_new[ij]= (c[ij] + c_periph*(c_new[ij-1]+c_new[ij+1]+c_new[ij-nx]+c_new[ij+nx])
                                +k*nu*f[ij]- k*omega_conso[ij]*n[ij])/c_center;
                            double diff=fabs(c_new[ij]-old);
                            if(diff>diff_max)diff_max=diff;
                        }
                        double old_a=a_new[ij];
                        double TAF_prod=0;
                        if(n[ij]>0 && c_new[ij]<0.1){
                            TAF_prod=k*alpha_TAF;
                        }
                        double a_center2=a_center+ k*0.1*v[ij];
                        a_new[ij]=(a[ij]+a_periph*(a_new[ij-1]+a_new[ij+1]+a_new[ij-nx]+a_new[ij+nx])+TAF_prod)/a_center2;
                    }

                }
            }
            for(int j=1; j<ny-1; j++){
                for(int i=1; i<nx-1; i++){
                    // Boundary
                    int ij=i+j*nx;
                    if(i==1) c_new[ij-1] = c_new[ij];
                    if(i==nx-2) c_new[ij+1] = c_new[ij];
                    if(j==1) c_new[ij-nx] = c_new[ij];
                    if(j==ny-2) c_new[ij+nx] = c_new[ij];
                }
            }
            if(diff_max<tol)break;
        }
        swap(c,c_new);
        swap(a,a_new);

        // Blood vessels update
        if(t_step>=10*2000){ // We wait for TAF to be diffused in the cell
            vector<Vessel> new_vessels;
            for(auto it=vessels.begin(); it!= vessels.end(); it++){
                if(!it->active)continue;

                int i=it->i;
                int j=it->j;
                int ij=i+j*nx;

                double chi_a=chi/(1+alpha*a[ij]);

                double P1 =k*dv/hx2 -(k/(4*hx2))*(chi_a*(a[ij+1]-a[ij-1]) +rho5*(f[ij+1]-f[ij-1]));
                double P2 =k*dv/hx2 +(k/(4*hx2))*(chi_a*(a[ij+1]-a[ij-1]) +rho5*(f[ij+1]-f[ij-1]));
                double P3 =k*dv/hx2 -(k/(4*hy2))*(chi_a*(a[ij+nx]-a[ij-nx]) +rho5*(f[ij+nx]-f[ij-nx]));
                double P4 =k*dv/hx2 +(k/(4*hy2))*(chi_a*(a[ij+nx]-a[ij-nx]) +rho5*(f[ij+nx]-f[ij-nx]));
                // double P0=1 -4*k*dv/hx2 +k*alpha*chi_a*((a[ij+1]-a[ij-1])*(a[ij+1]-a[ij-1])
                //         +(a[ij+nx]-a[ij-nx])*(a[ij+nx]-a[ij-nx]))/(4*hx2*(1+alpha*a[ij]))
                //         -k*chi_a*(a[ij+1]+a[ij-1]-4*a[ij]+a[ij+nx]+a[ij-nx])/hx2
                //         -k*rho5*(f[ij+1]+f[ij-1]-4*f[ij]+f[ij+nx]+f[ij-nx])/hx2;

                // if(t_step%2000==0)cout<<"|P0:"<<P0<<"|P1:"<<P1<<"|P2:"<<P2<<"|P3:"<<P3<<"|P4:"<<P4<<endl;

                if(P1<0) P1=0; if(P2<0) P2=0; if(P3<0) P3=0; if(P4<0) P4=0;
                double P0=1-(P1+P2+P3+P4);
                if(P0<0) P0=0;

                double sum=P0+P1+P2+P3+P4;
                double nb=loop(gen)*sum;
                int new_i=i, new_j=j;

                if(nb<P0){}
                else if(nb<P0+P1){new_i=i-1;}
                else if(nb<P0+P1+P2){new_i=i+1;}
                else if(nb<P0+P1+P2+P3){new_j=j-1;}
                else{new_j=j+1;}

                int new_ij=new_i+new_j*nx;
                if(new_i>0 && new_i<nx-1 && new_j>0 && new_j<ny-1 && n[new_ij]==0){ // Not boundary

                    // Anastomosis
                    if(v[new_ij]==1 && new_ij!=ij && new_ij!=1){
                        if(new_i==it->old_i && new_j==it->old_j){
                            new_i=i;
                            new_j=j;
                            new_ij=ij;
                        }
                        else{ // Anastomosis
                            it->active=false;
                            v[new_ij] =1;
                            c[new_ij] +=1e-3;
                        }
                    }

                    if(new_ij!=ij && it->active){
                        it->old_i=i;
                        it->old_j=j;
                        it->i=new_i;
                        it->j=new_j;
                        v[new_ij] =1;
                        c[new_ij] += 0.1;
                    }
                }

                // Branching
                it->age+=k;
                if(it->age>=0.2){
                    double prob_branch=a[ij]*a[ij]*a[ij]*p_angio; // Depend on TAF concentration
                    if(loop(gen)<prob_branch){
                        int new_i2=-1, new_j2=-1, new_ij2=-1;
                        if(v[ij+nx]==0){
                            new_i2=i;
                            new_j2=j+1;
                            new_ij2=ij+nx;
                        }
                        else if(v[ij-nx]==0){
                            new_i2=i;
                            new_j2=j-1;
                            new_ij2=ij-nx;
                        }
                        if(new_i2!=-1){
                            new_vessels.push_back({new_i2, new_j2, 0, true, i, j});
                            v[new_ij2]=1;
                            c[new_ij2]+=1e-3;
                            it->age=0;
                        }
                    }
                }
            }
            vessels.insert(vessels.end(), new_vessels.begin(), new_vessels.end());
        }

        new_cells.clear();
        // Moving cells
        for(int idx=0; idx < cells.size(); idx++){
            int i = cells[idx].i;
            int j = cells[idx].j;
            int ij = i+j*nx;

            // Necrosis?
            if(c[ij]<0.05){
                cells[idx].type = 5;
                n[ij]=0;

                dead_cells.push_back(cells[idx]);
                cells[idx] = cells.back();
                cells.pop_back();
                idx--;
                continue;
            }

            int Ae = n[ij-nx]+ n[ij+1]+n[ij+nx] +n[ij-1];
            int Ai=5; // Ignore quiescent cells
            if(cells[idx].type==1) Ai=3;
            if(cells[idx].type==2) Ai=2;
            if(cells[idx].type==3) Ai=1;
            if(cells[idx].type==4) Ai=1;

            if(Ae >= Ai){

                // Calculate probabilities of movement
                double P0_type, P_type;
                if(cells[idx].type==1){
                    P0_type = P0_coeff2_rho1;
                    P_type = P_coeff2_rho1;
                }
                else if(cells[idx].type==2){
                    P0_type = P0_coeff2_rho2;
                    P_type = P_coeff2_rho2;
                }
                else if(cells[idx].type==3){
                    P0_type = P0_coeff2_rho3;
                    P_type = P_coeff2_rho3;
                }
                else if(cells[idx].type==4){
                    P0_type = P0_coeff2_rho4;
                    P_type = P_coeff2_rho4;
                }
                P0 = P0_coeff1 - P0_type*(f[ij-1]+f[ij+1]-4*f[ij]+f[ij-nx]+f[ij+nx]);
                P1 = P_coeff1 -P_type*(f[ij+1]-f[ij-1]);
                P2 = P_coeff1 +P_type*(f[ij+1]-f[ij-1]);
                P3 = P_coeff1 -P_type*(f[ij+nx]-f[ij-nx]);
                P4 = P_coeff1 +P_type*(f[ij+nx]-f[ij-nx]);
                // if(P0<0) P0=0; if(P1<0) P1=0; if(P2<0) P2=0; if(P3<0) P3=0; if(P4<0) P4=0;

                // Moving cells
                double sum=P0+P1+P2+P3+P4;
                double nb = loop(gen)*sum;
                int new_i=i;
                int new_j=j;
                if(nb < P0){}
                else if(nb < P0+P1){new_i=i-1;}
                else if(nb < P0+P1+P2){new_i=i+1;}
                else if(nb < P0+P1+P2+P3){new_j=j-1;}
                else{new_j=j+1;}
                int new_ij=new_i+new_j*nx;
                if(new_i>0 && new_i<nx-1 && new_j>0 && new_j<ny-1){ // Boundary
                    if(new_ij != ij && n[new_ij] == 0){ // Empty cell
                        cells[idx].i=new_i;
                        cells[idx].j=new_j;
                        n[new_ij]=1;
                        n[ij]=0;
                        // Update for prolif
                        i=new_i;
                        j=new_j;
                        ij=new_ij;
                    }
                }
            }
            cells[idx].age++;

            int it_age = 0;
            if(cells[idx].type==1) it_age=it_type1;
            if(cells[idx].type==2) it_age=it_type2;
            if(cells[idx].type==3) it_age=it_type3;
            if(cells[idx].type==4) it_age=it_type4;

            // If the cell is mature
            if(it_age!=0 && cells[idx].age % it_age==0){
                if((n[ij-1]==0 && i>1) || (n[ij+1]==0 && i<nx-2) || (n[ij-nx]==0 && j>1) || (n[ij+nx]==0&&j<ny-2)){ // If empty neighbour cell

                    // Prolifération
                    int size_t = new_cells.size();
                    while(new_cells.size() == size_t){
                        int invasion = proliferate(gen);

                        if(invasion == 1 && n[ij-1] == 0 && i>1){
                            new_cells.push_back({i-1, j, 0, id_nb, cells[idx].type});
                            n[ij-1] += 1;
                        }
                        else if(invasion == 2 && n[ij+1] == 0 && i<nx-2){
                            new_cells.push_back({i+1, j, 0, id_nb, cells[idx].type});
                            n[ij+1] += 1;
                        }
                        else if(invasion == 3 && n[ij-nx] == 0 && j>1){
                            new_cells.push_back({i, j-1, 0, id_nb, cells[idx].type});
                            n[ij-nx] += 1;
                        }
                        else if(invasion == 4 && n[ij+nx] == 0 && j<ny-2){
                            new_cells.push_back({i, j+1, 0, id_nb, cells[idx].type});
                            n[ij+nx] += 1;
                        }
                    }
                    id_nb++;

                    if(cells[idx].type>=1 && cells[idx].type<4){
                        double nb = loop(gen);
                        if(nb<pmut){
                            cells[idx].type++;
                        }
                    }
                }
            }
        }

        cells.insert(cells.end(), new_cells.begin(), new_cells.end());

        // Saving
        if(t_step%saving_step==0){
            continuous.write(reinterpret_cast<char*>(f), nx*ny*sizeof(double));
            continuous.write(reinterpret_cast<char*>(m), nx*ny*sizeof(double));
            continuous.write(reinterpret_cast<char*>(c), nx*ny*sizeof(double));
            continuous.write(reinterpret_cast<char*>(a), nx*ny*sizeof(double));
            continuous.write(reinterpret_cast<char*>(v), nx*ny*sizeof(double));
            for(auto it = cells.begin(); it != cells.end(); ++it){
                discrete << t_step << " " << it->i << " " << it->j << " " << it->type << "\n";
            }
            for(auto it = dead_cells.begin(); it != dead_cells.end(); ++it){
                discrete << t_step << " " << it->i << " " << it->j << " " << it->type << "\n";
            }
            cout << t_step << endl;
        }
    }
    continuous.close();
    discrete.close();

    double t2=clock();
    double temps=(t2-t1)/CLOCKS_PER_SEC;
    std::cout << "Fin de simulation | Temps : " << temps << " secondes." << std::endl;

    delete [] n;
    delete [] f;
    delete [] m;
    delete [] m_new;
    delete [] c;
    delete [] c_new;
    delete [] omega_conso;
    delete [] a;
    delete [] a_new;
    delete [] v;

    cout << "Creation des animations" << endl;
    std::system("python -m plotter2");
    cout << "Fin" << endl;

    return 0;
}