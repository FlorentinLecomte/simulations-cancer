function [c,f,s] = modelpde(x,t,u,dudx)
dn=0.0005;
dm=0.0005;
dc=0.5;
rho=0.01;
eta=50;
kappa=1;
sigma=0;
nu=0.5;
omega=0.57;
phi=0.025;

c = [1; 1; 1; 1];
f = [dn*dudx(1) - rho*u(1)*dudx(2); 0; dm*dudx(3); dc*dudx(4)];
s = [0;-eta*u(3)*u(2); kappa*u(1)- sigma*u(3); nu*u(2)-omega*u(1)-phi*u(4)];
end