% File for the space domain and initial condition.

spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -5;		b = 5;		% x-domain
c = -5;		d = 5;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC


for i=1:nx
    for j =1:ny
        DCC_IC(i,j) = yc(j)>=0.7*sin(0.9*xc(i))+0.008*xc(i)^3+4;
        CSC_IC(i,j) = 0;
        FIB_IC(i,j) = 0;

    end
end

ECM_IC = 1 - DCC_IC;
MMP_IC = 0.05 * DCC_IC;
Fib_IC = 0 * PCC_IC;