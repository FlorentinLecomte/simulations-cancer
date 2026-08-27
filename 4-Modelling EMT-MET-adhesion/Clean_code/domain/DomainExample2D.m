% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -2;		b = 2;		% x-domain
c = -2;		d = 2;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC
                
                

% new flower/ball mix
epsilon = 0.3;
for i=1:nx
    for j =1:ny
        DCC_IC(i,j) =  exp(-xc(i)^2/epsilon-yc(j)^2/epsilon);
    end
end

ECM_IC = 1 - DCC_IC;
MMP_IC = 0.05 * DCC_IC;
Fib_IC = 0 * PCC_IC;