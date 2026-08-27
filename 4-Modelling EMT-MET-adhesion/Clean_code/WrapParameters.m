% very boring script that wraps numerical parameters in numP, model
% parameters in modP, and domain parameters in domP

numP.Nx = nx;
numP.CFL = CFL;

numP.mdt = max_delta_t;
numP.nos = nos;

numP.dx = delta_x;
numP.theta = theta;

if spDim == 2
    numP.dy = delta_y;
    numP.Ny = ny;
end

% for the a-levels
numP.alpha_min = alpha_min;
numP.alpha_max = alpha_max;
numP.alpha_num_levels = alpha_num_levels;
numP.alevels = alevels;
numP.da = delta_a;

numP.scIndex = scIndex;

numP.krylov = krylov;
numP.checkPos = checkPos;

numP.Ld = L_diff;

domP.Omega = Omega;

domP.a = a; domP.b = b;
domP.spDim = spDim;

if spDim == 2
    domP.c = c; domP.d = d;
end


%%%%%%%% adjust model parameters %%%%%%%%%%%%%%%%%%%%%
%for the (sigmoid wrt a) haptotaxis function \lamdba = L/(1 + exp(-k(x-x)0)) )
modP.Dc_M = Dc_M;
modP.Dc_sl = Dc_sl;
modP.Dc_a0 = Dc_a0;

% modP.xi_1 = xi_1;
% modP.xi_2 = xi_2;

modP.Dgf =Dgf;

%modP.mu = mu;

modP.ecm_degr_rate = ecm_degr_rate;
modP.degr_factor = degr_factor;
modP.ecm_prod_rate = ecm_prod_rate;

modP.prol_rate = prol_rate;
modP.prol_factor = prol_factor;

modP.gf_prod_rate = gf_prod_rate;
modP.gf_degr_rate = gf_degr_rate;

%for the (sigmoid wrt a) haptotaxis function \lamdba = L/(1 + exp(-k(x-x)0)) )
modP.hapto_L = hapto_L;
modP.hapto_sl = hapto_sl;
modP.hapto_a0 = hapto_a0;

% for the EMT & MET switches (gf dependent) %1/(1+exp(slp*(gf-gf0)))
% see reaction
modP.EMT_slp = EMT_slp;             modP.EMT_gf0 = EMT_gf0;   %increasing with inflection at EMTgf0
modP.MET_slp = MET_slp;             modP.MET_gf0 = MET_gf0;   %decreasing with inflection at METgf0

% EMT and MET functions (a dependent) % L*exp(-lambda*(a0-a).^2)
% see reaction
modP.EMT_L = EMT_L;         modP.EMT_gamma = EMT_gamma;
modP.MET_L = MET_L;         modP.MET_gamma = MET_gamma;


% Non-local term parameters
modP.Scv1=1;
R = 0.3;
