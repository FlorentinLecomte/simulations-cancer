% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -5e+3;		b = 5e+3;		% x-domain
c = -5e+3;		d = 5e+3;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                        % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC


epsilon = 0.1;
g = (1.5e+3)*sin( (8e-4)*xc(:)  )  + (2.4e-8)*xc(:).^3 + (2e+3);
for i=1:nx
    for j =1:ny
        PCC_IC(i,j) =... % exp(-xc(i)^2/epsilon-yc(j)^2/epsilon) ;
            ...%( xc(i)^2 + yc(j)^2 < 2 + (yc(j)<=0)*sin(5*atan2(yc(j),xc(i))) );%exp(-xc(i)^2/epsilon-yc(j)^2/epsilon);
            0.7*exp( -5*(yc(j)<g(i))*(yc(j)-g(i))^2 );
%             (yc(j)<g(i));
%        ECM_IC(i,j) = 0.5*( abs(xc(i)-1)<0.5 || abs(xc(i)+1)<0.5) ;% + 0.5*( abs(yc(j) )<0.5) ;%1 - PCC_IC(i,j);% + rand() * 0.5;

    end
end

% close, imagesc( PCC_IC )
% pause

% rescale the initial conditions to more "proper" values
% PCC_IC = 0.5*PCC_IC/max(PCC_IC(:));

MCC_IC = 0*0.05 * PCC_IC;

[zm,xm,ym] = TerrainNew(8,2);

zm = ( zm - min(zm(:)) )/( max(zm(:)) - min(zm(:)) ) ; 
zm = 0.4*zm + 0.05;
ECM_IC = zm;
ECM_IC(:) = 0.2;%1 - (PCC_IC(:)-MCC_IC(:));


AC_IC = 0*0.3 * PCC_IC;
Fib_IC = 0 * PCC_IC;
