function [c,f,s] = modelpde(x,t,u,dudx)
dn=0.001;
dm=0.001;
gamma=0.005;
eta=10;
alpha=0.1;
beta=0;

c = [1; 1; 1];
f = [dn*dudx(1) - gamma*u(1)*dudx(2); 0; dm*dudx(3)];
s = [0;-eta*u(3)*u(2); alpha*u(1)- beta*u(3)];
end