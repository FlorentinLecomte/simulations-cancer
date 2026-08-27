% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -5;		b = 5;		% x-domain
c = -5;		d = 5;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC


epsilon = 1;
for i=1:nx
    for j =1:ny
        PCC_IC(i,j) =  exp(-xc(i)^2/epsilon-yc(j)^2/epsilon) ;
            %( xc(i)^2 + yc(j)^2 < 2 + (yc(j)<=0)*sin(5*atan2(yc(j),xc(i))) );%exp(-xc(i)^2/epsilon-yc(j)^2/epsilon);
        ECM_IC(i,j) = 0.5*( abs(xc(i)-1)<0.5 || abs(xc(i)+1)<0.5) ;% + 0.5*( abs(yc(j) )<0.5) ;%1 - PCC_IC(i,j);% + rand() * 0.5;
    end
end


MCC_IC = 0.05 * PCC_IC;
AC_IC = 0.3 * PCC_IC;
Fib_IC = 0 * PCC_IC;