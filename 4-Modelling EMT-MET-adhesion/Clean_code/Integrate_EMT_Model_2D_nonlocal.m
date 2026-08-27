function [S, C] = Integrate_EMT_Model_2D_nonlocal2(S_start, C_start, modP, numP, domP, t_start, T, FNonLocal)
% Function to integrate the system from t_start to T with initial data
% S_start in 2D using an IMEX-RK scheme.

% in S we keep the components of the solutions, i.e.
% gf = S(:,:,1);
% ecm = S(:,:,2);
% fib = S(:,:,3);


% unwrap numerical parameters
delta_x = numP.dx;
delta_y = numP.dy;
max_delta_t = numP.mdt;
CFL = numP.CFL;
nx = numP.Nx;
ny = numP.Ny;
% aa_levels = nP.aa_levels;
alpha_num_levels = numP.alpha_num_levels;

% initialize variables and load scheme coefficients
t = t_start;
S = S_start;
C = C_start;
d = max(delta_x,delta_y);

[Ae,be,we,A,~,w] = GetScheme(numP.scIndex);


while t<T
    delta_S = zeros(nx, ny, modP.comps);
    delta_C = zeros(nx, ny, alpha_num_levels);

    RKstages = size(A,1);			% number of RK stages
    K = zeros( nx, ny, modP.comps, RKstages);
    f = zeros( nx, ny, modP.comps, RKstages);
    g = zeros( nx, ny, modP.comps, RKstages);

    KC = zeros( nx, ny, alpha_num_levels, RKstages);
    fC = zeros( nx, ny, alpha_num_levels, RKstages);
    gC = zeros( nx, ny, alpha_num_levels, RKstages);

    for stage=1:RKstages				% iterate through the stages (IMEX-RK)
        temp = S;
        tempC = C;
        % computation of the preliminary (purely explicit) RK stage temp
        for j=1:stage-1
            temp = temp + delta_t*( Ae(stage,j)*f(:,:,:,j) + A(stage,j)*g(:,:,:,j) );
            tempC = tempC + delta_t*( Ae(stage,j)*fC(:,:,:,j) + A(stage,j)*gC(:,:,:,j) );
        end
        % computation of the final (partly implicit) RK stage K(:, :, i)
        if stage==1
            K(:,:,:,1) = S;
            KC(:,:,:,1) = C;
            [g(:,:,:,1), gC(:,:,:,1)] = EXPLICIT_DIFFUSION(S, C, modP, numP, domP);
            time = t;
            [f(:,:,:,1), fC(:,:,:,1), a] = EXPLICIT_ADVECTION_REACTION(K(:,:,:,1), KC(:,:,:,1), time, modP, numP, domP, FNonLocal);
            delta_t = min([CFL*d/(2*a), T-t, max_delta_t]);
        else
            time = t + be(stage) * delta_t;
            [K(:,:,:,stage), g(:,:,:,stage),KC(:,:,:,stage), gC(:,:,:,stage)] = IMPLICITPART(temp, tempC, delta_t*A(stage,stage), modP, numP, domP);
            [f(:,:,:,stage), fC(:,:,:,stage),~] = EXPLICIT_ADVECTION_REACTION(K(:,:,:,stage), KC(:,:,:,stage), time, modP, numP, domP, FNonLocal);
        end
        delta_S = delta_S + delta_t*(we(stage)*f(:,:,:,stage)+w(stage)*g(:,:,:,stage));
        delta_C = delta_C + delta_t*(we(stage)*fC(:,:,:,stage)+w(stage)*gC(:,:,:,stage));
    end

    S = S + delta_S;
    C = C + delta_C;
    t = t + delta_t;

    if numP.checkPos
        if ( min(S(:))<0 ) || ( min(C(:))<0 )
            [ min(S(:)), min(C(:)) ]
            error('Solution became negative.');
        end
    end
end % while t<T


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Y,YC,a] = EXPLICIT_ADVECTION_REACTION(S, C, time, modP, numP, domP, FNonLocal)
% function to evaluate the explicit term i.e. advection & reaction

%% start with the reactions
[Y ,YC] = Reaction(S, C, time, modP, numP, domP, 2);

%% unwrap parameters
theta = numP.theta; delta_x = numP.dx; delta_y = numP.dy;
nx = numP.Nx; ny = numP.Ny;
hapto_L = modP.hapto_L;
hapto_sl = modP.hapto_sl;
hapto_a0 = modP.hapto_a0;
alpha_num_levels = numP.alpha_num_levels;

%% now the taxis terms
% CC = sum(C,3);
% GF = S(:,:,1);
ECM = S(:,:,2);
% FIB = S(:,:,3);

% adiff = numP.alevels - hapto_a0;
% hapto = hapto_L./( 1 + exp(-hapto_sl*adiff) ) ;
% clf;
% plot(adiff, hapto)
% pause
a = 1e-6; % impact on delta_t = min([CFL*d/(2*a), T-t, max_delta_t]);

% clf;
%the haptotaxis coefficient (only a dependent, no x or y)
adiff = numP.alevels - hapto_a0;
hapto = hapto_L./( 1 + exp(-hapto_sl*adiff) ) ;
% plot(adiff, hapto);
% pause

kappa2=1;

d_alpha=1/alpha_num_levels;
o=sum(C,3)*d_alpha+ECM; %sum(,3) we sum on the phenotype to get a 2D matrix

Scv1=modP.Scv1;
Scv2=1;
Scv3=3/(pi*FNonLocal.R^2);
Scv=Scv1*Scv2*Scv3;

G=Scv * (1-o).^kappa2 .* ECM;
G_translate=G.';

% (ii) Evaluate the nonlocal term for G on all grid cell interfaces
% (in the call, ruleID = 3 is used, the fast FFT-based evaluation)
%fprintf('\nStarting evaluation of nonlocal term...\n'); %tstart=tic;
[A11, A22] = evalIntegral2D(G_translate, FNonLocal.mask, 3);
A_cvx=A11.';
A_cvy=A22.';
%elapsed=toc(tstart);
%fprintf('Completed evaluation of nonlocal term in %d second.\n', elapsed)
%fprintf('Dimension of A11: (%d,%d)\n', size(A11));
%fprintf('Dimension of A22: (%d,%d)\n', size(A22));

for a_level = 1:alpha_num_levels
    % haptotaxis
        %the haptotaxis coefficient (only a dependent, no x or y)
        % adiff = numP.alevels(a_level) - hapto_a0;
        % hapto = hapto_L./( 1 + exp(-hapto_sl*adiff) ) ;

        %... on the interfaces
        hapto_xif = hapto(a_level);%(hapto(1:end-1, :) + hapto(2:end, :)) / 2;
        hapto_yif = hapto(a_level);%(hapto(:, 1:end-1) + hapto(:, 2:end)) / 2;
        % derivatives for characteristic velocities; with grad ECM (no variable of the local eq.)
        % DVx = [zeros(1,ny);	hapto_xif .* diff(ECM, 1, 1)/delta_x;	zeros(1,ny)];
        % DVy = [zeros(nx,1);	hapto_yif .* diff(ECM, 1, 2)/delta_y;	zeros(nx,1)];

        DVx =	hapto_xif .* A_cvx;	
        DVy = 	hapto_yif .* A_cvy 	;

        DFx = DVx; DFy =  DVy ;
        [nF_MCCx, nF_MCCy] = NUMFLUX(C(:,:,a_level), DFx, DFy, theta);
        % compute flux differences for the RHS i.e. the  - \nabla (F)
        Flux_MCC =	-1/delta_x * diff(nF_MCCx, 1, 1) ...
            -1/delta_y * diff(nF_MCCy, 1, 2);
        a = max( a, max(abs( [DVx(:); DVy(:)] )) );
        YC(:, :, a_level) = Flux_MCC + YC(:, :, a_level);
end % for a_level


    % cross-taxis for the GFs (BEWARE OF THE SIGN OF cross_tax.)
%     cross_tax = -1 * xi_2 * (1 + ECM) ./ (1 + (CC + GF) + ECM).^2;
%     cross_tax_xif = (cross_tax(1:end-1, :) + cross_tax(2:end, :)) / 2;
%     cross_tax_yif = (cross_tax(:, 1:end-1) + cross_tax(:, 2:end)) / 2;
%     DPx = [zeros(1,ny);	cross_tax_xif .* diff(GF, 1, 1)/delta_x;	zeros(1,ny)];
%     DPy = [zeros(nx,1)	cross_tax_yif .* diff(GF, 1, 2)/delta_y	zeros(nx,1)];

    % acidotaxis
%     ac_tax = -xi_3 * (1 + ECM) ./ (1 + (CC + GF) + ECM).^2;
%     ac_tax_xif = (ac_tax(1:end-1, :) + ac_tax(2:end, :)) / 2;
%     ac_tax_yif = (ac_tax(:, 1:end-1) + ac_tax(:, 2:end)) / 2;
%     DHx = [zeros(1,ny);	ac_tax_xif .* diff(FIB, 1, 1)/delta_x;	zeros(1,ny)];
%     DHy = [zeros(nx,1)	ac_tax_yif .* diff(FIB, 1, 2)/delta_y	zeros(nx,1)];

    % NUMFLUX on the summ of all fluxes (just haptotaxis in this case)
    % DFx = DVx + DPx + DHx; DFy =  DVy + DPy + DHy;
    % [nF_MCCx, nF_MCCy] = NUMFLUX(CC, DFx, DFy, theta);

% compute flux differences for the RHS i.e. the  - \nabla (F)
%     Flux_MCC =	-1/delta_x * diff(nF_MCCx, 1, 1) ...
%         -1/delta_y * diff(nF_MCCy, 1, 2);

% compute maximal speed
    %a = max( abs( [DVx(:); DVy(:)]  + [DPx(:); DPy(:)] ) );
%     a = max( abs( [DVx(:); DVy(:)]  ) );

% add reaction and advection terms
% Y(:,:,1) = Flux_MCC	+ Y(:,:,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [K,g, KC, gC]= IMPLICITPART(temp, tempC, d, modP, numP, domP)
% function for the implicit part of the IMEX that solves
% (Id - d * Impl) K = temp. The Impl term of the system is the diffusion
% part.

% unwrap some things
nx = numP.Nx;
ny = numP.Ny;
% L_diff = numP.Ld;
%
% Dc_L = modP.Dc_L;
% Dc_sl = modP.Dc_sl;
% Dc_a0 = modP.Dc_a0;
%
% Dgf = modP.Dgf;


pr = nx*ny;
Id = speye(pr);
K = temp;%zeros(nx,ny,mP.comps);
KC = tempC;

% cc = sum(tempC,3);
% gf = temp(:, :, 1);
% ecm = temp(:, :, 2);
% fib = temp(:, :, 3);


%first for K i.e. S i.e. no cancer cells
for equation = [1:modP.comps]
    % calculate the discrete diffusion for very system equation
    [discrete_diff, exists] = Diffusion_discrete(K, KC, modP, numP, domP, 'system', equation);
    if exists == 1 %if this equation has diffusion,
    % solve with Krylov or otherwise
    if numP.krylov
        [tt, ~] = bicgstab((Id-d*(discrete_diff)) ,reshape(temp(:,:,equation),pr,1));
    else
        tt = (Id-d*(discrete_diff))\reshape(temp(:,:,equation),pr,1);
    end

    K(:,:,equation) = reshape(tt,nx,ny);
    end %if exists
end % for equation loop


% now for KC i.e. C(:,:,alpha) alpha: a_level of mesenchymality
for a_level = 1:numP.alpha_num_levels
    % calculate the discrete diffusion for every a level
    [discrete_diff, exists] = Diffusion_discrete(K, KC, modP, numP, domP, 'cancer', a_level);
    if exists == 1
    if numP.krylov
        [tt, ~] = bicgstab((Id-d*(discrete_diff)) ,reshape(tempC(:,:,a_level),pr,1));
    else
        tt = (Id-d*(discrete_diff))\reshape(tempC(:,:,a_level),pr,1);
    end
    KC(:,:,a_level) = reshape(tt,nx,ny);
    end % if exists
end
% evaluate the Impl term (diffusion) for the new K
[g, gC] = EXPLICIT_DIFFUSION(K, KC, modP, numP, domP);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [g, gC] = EXPLICIT_DIFFUSION(S, C, modP, numP, domP)
% function to evaluate the implicit term (diffusion) for S

nx = numP.Nx;
ny = numP.Ny;
alpha_num_levels = numP.alpha_num_levels;

Dc_M = modP.Dc_M;
Dc_sl = modP.Dc_sl;
Dc_a0 = modP.Dc_a0;

Dgf = modP.Dgf;

pr = nx*ny;

L_diff = numP.Ld;

g = zeros(nx, ny, modP.comps);
gC = zeros(nx, ny, alpha_num_levels);


% cc = sum(C,3);
% gf = S(:, :, 1);
% ecm = S(:, :, 2);
% fib = S(:, :, 3);
        %if non-lienar diffusion
%         f = Dc * (~mP.degDiff + cc.*gf + cc.*ecm  + gf .* ecm ) ./ ...
%            ( 1+ cc .* (gf + ecm) );
%         f = Dgf*(0*gf+1);
%         discrete_diff = GiveDiffusionMatrixFV([nP.Nx,nP.Ny],...
%             [nP.dx,nP.dy],2, f); % non-linear diffusion
        %if linear diffusion
%         discrete_diff = Dgf*L_diff;
%for the non-cancerous componens
for equation = [1]  %tackle these equations
    if equation == 1
        % fib = S(:, :, 3);
        %if non-lienar diffusion
%         f = Dc * (~mP.degDiff + cc.*gf + cc.*ecm  + gf .* ecm ) ./ ...
%            ( 1+ cc .* (gf + ecm) );
%         f = Dgf*(0*gf+1);
%         discrete_diff = GiveDiffusionMatrixFV([nP.Nx,nP.Ny],...
%             [nP.dx,nP.dy],2, f); % non-linear diffusion
        %if linear diffusion
%        discrete_diff = Dgf*L_diff;
        discrete_diff = Diffusion_discrete(S, C, modP, numP, domP, 'system', equation);
    end
    curr_comp = reshape(S(:,:,equation),pr,1);
    g(:,:,equation) = reshape(discrete_diff*curr_comp,nx,ny);
end

%for the cancer components
% now for gC i.e. C(:,:,alpha) alpha: level of mesenchymality
for a_level = 1:numP.alpha_num_levels       %tackle all the alpha a_levels
    % nonlinear diffusion
        %f = Dc * (~mP.degDiff + cc.*gf + cc.*ecm  + gf .* ecm ) ./ ...
        %      ( 1+ cc .* (gf + ecm) );
        %discrete_diff = GiveDiffusionMatrixFV([nP.Nx,nP.Ny],...
        %    [nP.dx,nP.dy],2, f);
     %linear diffusion
%     discrete_diff = Dc_L./( 1 + exp(-Dc_sl*(nP.alevels(a_level)-Dc_a0)) )*L_diff;
    discrete_diff = Diffusion_discrete(S, C, modP, numP, domP, 'cancer', a_level);
    curr_comp = reshape(C(:,:,a_level),pr,1);
    gC(:, :, a_level) = reshape(discrete_diff*curr_comp,nx,ny);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nFx, nFy] = NUMFLUX(U, DFx, DFy, theta)
% function to compute the numerical flux for U wrt. the velocities
% DFx and DFy, using the theta-minmod limiter
[nx, ny] = size(U);

U_xtnd=[    nan		U(1,:)		nan;
            U(:,1)	U			U(:,ny);
            nan		U(nx,:)	    nan];
% The matrix U is in the form
%	  | y=1 | y=2 | y=3 | .....
%------------------------------------
% x=1 |		|	  |
% x=2 |		|	  |
% x=3 |		|	  |
% ... |		|	  |

% computation of minmod slopes DU
DU = zeros(size(U_xtnd,1)-2, size(U_xtnd,2)-2, 2);

% use two copies, one for the y- and one for the x-direction Ufx, Ufy and
% compute the minmod slopes in both directions
Ufx = U_xtnd(:,2:end-1);
tdU = theta * diff(Ufx, 1, 1);
fU = tdU(1:end-1,:); bU = tdU(2:end, :);
cU =(Ufx(3:end,:)-Ufx(1:end-2,:))/2;

sdU = sign(tdU);
sfU = sdU(1:end-1, :); sbU = sdU(2:end, :);
DU(:,:,1) = (sfU==sbU).*sfU.*min(min(abs(fU),abs(bU)),abs(cU));

% along the y-direction
Ufy = U_xtnd(2:end-1,:);
tdU = theta * diff(Ufy, 1, 2);
fU= tdU(:, 1:end-1); bU = tdU(:, 2:end);
cU=(Ufy(:,3:end)-Ufy(:,1:end-2))/2;

sdU = sign(tdU);
sfU = sdU(:, 1:end-1); sbU = sdU(:, 2:end);
DU(:,:,2) = (sfU==sbU).*sfU.*min(min(abs(fU),abs(bU)),abs(cU));

% compute interface values
U_xplus = U + DU(:,:,1)/2;
U_xminus= U - DU(:,:,1)/2;
% and the same along the y direction (ONLY FOR UNIFORM GRID)
U_yplus = U + DU(:,:,2)/2;
U_yminus= U - DU(:,:,2)/2;

% Upwinding
POSITIVEx = (DFx > 0);
nFx = DFx.*( POSITIVEx .* [U(1, :); U_xplus] ...
    + ~POSITIVEx.*[U_xminus; U( nx, :)] );

POSITIVEy = (DFy > 0);
nFy = DFy .* ( POSITIVEy .* [U(:,1)  U_yplus] ...
    + ~POSITIVEy .*[U_yminus  U(:, ny)] );
end

