#include <iostream>
#include <fstream>
using namespace std;

const double Et=100;
const double Pt=0.33514;

int main(){
    double t_max=1;
    double dt=0.0001;
    double Ec, Em, Eb, b, C, P;
    double dt_Ec, dt_Em, dt_Eb, dt_Et, dt_b, dt_C;
    double nu=100, kplus=100, kminus=19, km=0.01, k2=0.03;
    double c=0, d=0;
    double cT=50,a=2;

    Ec=0;
    Eb=100;
    C=km/k2;
    b=(km/k2*(kminus+k2))/(kplus*(Pt-km/k2));

    ofstream cadherin("cadherin.bin", ios::binary);

    for(double t=0; t<t_max/dt; t++){
        if(t>=0.4/dt && t<=0.8/dt)d=0.1; //10
        else d=0;

        if(b>cT){
            dt_Eb=-(a+d)*Eb;
            dt_b=(a+d)*Eb-kplus*b*(Pt-C)+kminus*C+km;
        }
        else{
            dt_Eb=nu*(Et-Ec-Eb)*b-d*Eb;
            dt_b=-nu*(Et-Ec-Eb)*b+d*Eb-kplus*b*(Pt-C)+kminus*C+km;

        }
        dt_Ec=-c*Ec+d*Eb;
        dt_C=kplus*b*(Pt-C)-kminus*C-k2*C;

        Ec+=dt*dt_Ec;
        Eb+=dt*dt_Eb;
        b+=dt*dt_b;
        C+=dt*dt_C;

        if(int(t)%100==0){
            double t_actuel=t*dt;
            double Eb_adim=Eb/Et, C_adim=C/Et, Ec_adim=Ec/Et, b_adim=b/Et;
            cadherin.write(reinterpret_cast<const char*>(&t_actuel), sizeof(double));
            cadherin.write(reinterpret_cast<const char*>(&Eb_adim), sizeof(double));
            cadherin.write(reinterpret_cast<const char*>(&C_adim), sizeof(double));
            cadherin.write(reinterpret_cast<const char*>(&Ec_adim), sizeof(double));
            cadherin.write(reinterpret_cast<const char*>(&b_adim), sizeof(double));
        }
    }

    cadherin.close();

    system("python -m plotter");

    return 0;
}