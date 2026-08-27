
%% cancer cells params
    % Diffusion: sigmoid (wrt a) function i.e L/(1 + exp(-k(x-x_0)) )
    %D_c^M c^2 1/( 1 + exp(-Dc_sl*(a-Dc_a0)) )
    Dc_M = 0*1e-7;   %maximum diffusion coefficient
    Dc_sl = 20;    %maximum slope of the sigmoid wrt alpha
    Dc_a0 = 0.5;   %sigmoid midpoint of the sigmoid wrt alpha

    % Haptotaxis: sigmoid (wrt a) function i.e. L/(1 + exp(-k(x-x_0)) )
    hapto_L = 0.5;%4e-1;%0.1;      %maximum haptotaxis coefficient
    hapto_sl = 50;    %maximum slope of the logistic
    hapto_a0 = 0.5;   %sigmoid midpoint

    % EMT & MET (gf & a dependent)
    % EMT:: EMT_L/(1+exp(-EMT_slp*(gf-gf0))) \( \int_0^a k(x,a,a')c(t,x,a')da' - c(t,x,a)\int_a^1 k(x,a,a')da' \)
    % k(x,a,a')= exp(-EMT_gamma*(a-a').^2)
    EMT_L = 10*1e-1;    %maximum EMT coefficient
    EMT_slp = 1e+3;             
    EMT_gf0 = 1e-2;%1e-2;   %increasing with inflection at EMTgf0
    EMT_gamma = 20;  % exponent coeff. 
    
    % MET:: MET_L/(1+exp( MET_slp*(gf-gf0))) \( \int_0^a k(x,a,a')c(t,x,a')da' - c(t,x,a)\int_a^1 k(x,a,a')da' \)
    % k(x,a,a')= exp(-MET_gamma*(a-a').^2)
    MET_L = 1;%200*EMT_L;    %maximum MET coefficient
    MET_slp =  1e+4;             
    MET_gf0 = 1e-3;   %decreasing with inflection at METgf0
    MET_gamma = 20;   %exponent coeff. 
    
    % cell proliferation
    % prol_rate*exp(-prol_factor*(alevels.^2))*cca*( 1 - sum(XC,3)*da - ecm )
    prol_rate = 1e+0;%5;                
    prol_factor = 7;                    % (for the exponential; see reaction)

%% growth factors
    Dgf = 2e+1*1e-2;                    % diffusion
    gf_prod_rate = 3e-1;
    gf_degr_rate = 3e-1;           


%% ecm degradation
    ecm_degr_rate = 1e-2;%5;       % rate (check the paper for the proper coefficie)
    degr_factor = 1e+2;           
    ecm_prod_rate = 0.0001;


