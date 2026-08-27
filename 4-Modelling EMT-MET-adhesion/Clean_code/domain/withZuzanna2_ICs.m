% File for the space domain and initial condition.

% domain parameters
spDim = 2;              % space dimension (1 or 2)
% domain will be [a,b] in 1D or [a,b]x[c,d] in 2D
a = -5;		b = 5;		% x-domain
c = -5;		d = 5;		% y-domain (needed if spDim = 2)

%a-variable 
alpha_min = 0;   alpha_max = 1;       % alpha-level bounds
alpha_num_levels = 11;           % number of a-levels to discretise

alevels = linspace(alpha_min, alpha_max, alpha_num_levels);
delta_a = alevels(2)-alevels(1);


modP.comps = 2;
modP.compNames = char('GFs', 'ECM');%, 'FIB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CellCenters;    % this script creates xc (and yc) and zero vecors/matrices
                % DCC_IC, CSC_IC, FIB_IC, ECM_IC, and MMP_IC

epsilon = 0.3;
for i=1:nx
    for j =1:ny
        CC_IC(i, j) =  100*exp(-(xc(i)^2+yc(j)^2)/epsilon);
        if xc(j)^2+yc(i)^2 > 0.1
            CC_IC(i, j) = 0;
        end
        ECM_IC(i, j) = 1 - CC_IC(i, j);
    end
end

% extarct the extremes of the above IC ECM (only for consistency)
% minECM = min(ECM_IC(:));  maxECM = max(ECM_IC(:));
% create new ECM by tintroducing randomness
[ECM_IC, ~, ~] = TerrainNew( log2(size(CC_IC,1)), 3);
% rescale it between the min and max of the ECM that you previously extracted
% for i=1:nx
%     for j =1:ny
%         ECM_IC(i,j) = xc(j)+yc(i);
%     end
% end

ECM_IC = (ECM_IC - min(ECM_IC(:)))/( max(ECM_IC(:)) - min(ECM_IC(:)) ); % scale between 0 and 1
% ECM_IC = minECM + ECM_IC*(maxECM-minECM); % scale between minECM and maxECM


CC_IC = 0.02 * CC_IC;   % 1st component of S matrix
GF_IC = 8e-1*CC_IC/max(CC_IC(:));    % 2nd component of S matrix
%the FIB and the ECM share the same tissue 
% FIB_IC = 0.8*ECM_IC;    % 3nd component of S matrix
ECM_IC = 0.2*ECM_IC;

%setup the IC matrix
CC_IC_all = zeros(nx, ny, alpha_num_levels);
% for level = 1:numlevels
%     CC_IC_all(: , : , level) = CC_IC; 
% end
CC_IC_all(: , : , 1) = CC_IC; 

clear CC_IC
