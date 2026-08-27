% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -2;		b = 2;		% x-domain
c = -2;		d = 2;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC

% Proliferating CCs
PCC_IC = zeros(nx,ny);

epsilon = 0.3;
for i=1:nx
    for j =1:ny
        PCC_IC(i,j) =  exp(-(xc(i)^2+yc(j)^2)/epsilon);
            if (yc(j) >= xc(i)-0.1 && yc(j) <= xc(i)+0.1) || (xc(i) >= -0.05 && xc(i) <= 0.05)
%                 ECM_IC(i,j) = 1 ;
                PCC_IC(i,j) = 0 ;
%             else
%                 ECM_IC(i,j) = 1 - PCC_IC(i,j);
%                 ECM_IC(i,j) = 1 - PCC_IC(i,j);
            end
    end
end

% Migrating CCs
MCC_IC = 0.05 * PCC_IC;

% ECM (with positivity constraint)
ECM_IC = max(1 - PCC_IC - MCC_IC, 0) ;

% Fibroblast cells
Fib_IC = 0 * PCC_IC;

% Acidity
AC_IC = PCC_IC + 0.2;
