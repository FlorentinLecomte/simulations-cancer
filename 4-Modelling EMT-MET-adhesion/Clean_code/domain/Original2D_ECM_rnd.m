% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -2;		b = 2;		% x-domain
c = -2;		d = 2;		% y-domain (needed if spDim = 2)

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC


epsilon = 0.3;
for i=1:nx
    for j =1:ny
        PCC_IC(i,j) =  exp(-(xc(i)^2+yc(j)^2)/epsilon);
        ECM_IC(i,j) = 1 - PCC_IC(i,j);
    end
end

% extarct the extremes of the above IC ECM (only for consistency)
minECM = min(ECM_IC(:));  maxECM = max(ECM_IC(:));
% create new ECM by tintroducing randomness
[ECM_IC, ~, ~] = TerrainNew( log2(size(PCC_IC,1)), 3);
% rescale it between the min and max of the ECM that you previously
% extracted
ECM_IC = (ECM_IC - min(ECM_IC(:)))/( max(ECM_IC(:)) - min(ECM_IC(:)) ); % scale between 0 and 1
ECM_IC = minECM + ECM_IC*(maxECM-minECM); % scale between minECM and maxECM

    


MCC_IC = 0.05 * PCC_IC;
AC_IC = PCC_IC + 0.2;
Fib_IC = 0 * PCC_IC;