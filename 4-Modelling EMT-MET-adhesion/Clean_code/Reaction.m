function [Y, YC] = Reaction(S, SC, t, modP, numP, domP, spDim)
% function to evaluate the reaction terms of the System at X
% works for 1D, 2D and 0D input data

Y = 0 * S;
YC = 0 * SC;

% unwrap parameters
%mu = modP.mu; %R_0 = mP.R_0; b = mP.b; y_ref = mP.y_ref;
alevels = numP.alevels; 
da = numP.da;

ECM_degr_rate = modP.ecm_degr_rate;       ECM_degr_factor =  modP.degr_factor;
ECM_prod_rate = modP.ecm_prod_rate;
prol_rate = modP.prol_rate;               prol_factor =  modP.prol_factor;

gf_prod_rate = modP.gf_prod_rate;
gf_degr_rate = modP.gf_degr_rate;

EMT_slp = modP.EMT_slp;         EMT_gf0 = modP.EMT_gf0;   %increasing with inflection at EMTgf0
MET_slp = modP.MET_slp;         MET_gf0 = modP.MET_gf0;   %decreasing with inflection at METgf0

EMT_L = modP.EMT_L;             EMT_gamma = modP.EMT_gamma;
MET_L = modP.MET_L;             MET_gamma = modP.MET_gamma;

%% extract components
switch spDim
    case 1
        gf = S(:, 1);
        ecm = S(:, 2);
%         fib = X(:, 3);
    case 2
        gf = S(:, :, 1);
        ecm = S(:, :, 2);
%         fib = X(:, :, 3);
    otherwise
        error('Space dimension has to be 0, 1, or 2.');
end

%% Cancer cell reaction terms
    % proliferation rate cancer cells (alpha dependent component) 
    % the coeff of prolif is added later
    cell_prol = prol_rate*exp(-prol_factor*(alevels.^2));  %size of alevels
%     cell_prol_2 = 2*prol_rate./(1+exp(0.4*prol_factor*(alevels)));
% clf;
% plot(alevels, cell_prol, alevels, cell_prol_2)
% pause
% cell_prol= prol_rate*(alevels<0.5);


%clf; plot(alevels, cell_prol), pause
    [~,~,cell_prol_alpha] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , cell_prol );

    % the gf dependence of the EMT/MET
    EMT_gf_dep = EMT_L./(1+exp(-EMT_slp*(gf-EMT_gf0)));
    MET_gf_dep = MET_L./(1+exp(MET_slp*(gf-MET_gf0))); %NOTE the "missing" "-" sign in exp()

    % EMT_gf_dep = EMT_L*(gf>=EMT_gf0);
    % MET_gf_dep = MET_L*(gf<MET_gf0);

% clf; 
% GF = linspace(0,0.5,1e+3); 
% EMT_GF_DEP = EMT_L./(1+exp(-EMT_slp*(GF-EMT_gf0)));
% EMT_NEW = EMT_L*(GF>EMT_gf0);
% plot(GF, EMT_GF_DEP); hold on; plot(GF, EMT_NEW);
% pause


    % now for XC i.e. C(:,:,alpha), alpha the level of mesenchynality
    for alpha = 1: numP.alpha_num_levels % go through the levels, one after the other
        switch spDim % consider the dimensions
            case 1
                cca = SC(:, alpha);
            case 2
                cca = SC(:, :, alpha);
            otherwise
                error('Space dimension has to be 0, 1, or 2.');
        end

        % create the 3 dimensional matrix with 
        adiff = numP.alevels(alpha) - numP.alevels;

%         % the gf dependence of the EMT/MET
%         EMT_gf_dep = EMT_L./(1+exp(-EMT_slp*(gf-EMT_gf0)));
%         MET_gf_dep = MET_L./(1+exp(MET_slp*(gf-MET_gf0))); %NOTE the "missing" - in exp()
% % EMT_gf_dependence = EMT_L*(gf>=EMT_gf0);
% % MET_gf_dependence = MET_L*(gf<MET_gf0);

        % the EMT/MET kernel dependence on alpha...
        EMT_alpha_dep = exp(-EMT_gamma*adiff.^2);
        EMT_alpha_dep_on= 3./(1+exp(0.5*EMT_gamma*(adiff)));
        EMT_alpha_dep_off = 2./(1+exp(-0.5*EMT_gamma*(adiff)));

        MET_alpha_dep = exp(-MET_gamma*adiff.^2);
        MET_alpha_dep_on= 3./(1+exp(-0.5*MET_gamma*(adiff)));
        MET_alpha_dep_off = 2./(1+exp(0.5*MET_gamma*(adiff))); 
% clf;
% % plot(adiff, MET_alpha_dep.*(adiff<0))
% % pause
% nexttile;
% plot(adiff,MET_alpha_dep.*(adiff<0), adiff, MET_alpha_dep_on.*(adiff<0));
% nexttile;
% plot(adiff,MET_alpha_dep.*(adiff>=0), adiff, MET_alpha_dep_off.*(adiff>=0));
% pause
 
        % which is directional to account for the EMT 
        [~,~,EMT_on] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , EMT_alpha_dep_on.*(adiff>0) );
        [~,~,EMT_off] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , EMT_alpha_dep_off.*(adiff<=0) );
        % directional MET kernel
        [~,~,MET_on] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , MET_alpha_dep_on.*(adiff<0) );
        [~,~,MET_off] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , MET_alpha_dep_off.*(adiff>=0) );

    % combine EMT, MET, and proliferation
        cca_new = EMT_gf_dep.*( sum(SC.*EMT_on,3)*da - cca.*sum(EMT_off,3)*da ) ...
              + MET_gf_dep.*( sum(SC.*MET_on,3)*da - cca.*sum(MET_off,3)*da ) ...
              + cell_prol_alpha(:,:,alpha).*cca.*( 1 - sum(SC,3)*da - ecm );
        switch spDim
            case 1
                YC(:, alpha) = cca_new;
            case 2
                YC(:, :, alpha) = cca_new;
            otherwise
                error('Space dimension has to be 0, 1, or 2.');
        end
    end %for alpha=1:numlevels


%


%% reaction terms of the non-cancerous part of the system
    % GF
    gf_react = gf_prod_rate.* ecm .* sum(SC,3)*da - gf_degr_rate*gf;
    % ECM
    % ecm_degr_kernel = exp(-ECM_degr_factor*((alevels-alevels(end)).^2)); 
    ecm_degr_kernel = 1./(1+exp(-ECM_degr_factor*(alevels-0.5))); 

% clf; 
% plot(alevels, ecm_degr_kernel);%, alevels, ecm_degr_kernel_new)
% pause

    [~,~,ecm_eta] =  ndgrid( zeros(numP.Ny,1) , zeros(1,numP.Nx) , ecm_degr_kernel );
    ecm_react = - ECM_degr_rate*ecm.*sum( ecm_eta.*SC, 3)*da + ECM_prod_rate*ecm .*(1 - sum(SC,3)*da - ecm);
    
%% transfer to appropriate output
switch spDim
    case 1
        Y(:, 1) = gf_react;
        Y(:, 2) = ecm_react;
%         Y(:, 3) = fib_new;
    case 2
        Y(:, :, 1) = gf_react;
        Y(:, :, 2) = ecm_react;
%         Y(:, :, 3) = fib_new;
    otherwise
        error('Space dimension has to be 0, 1, or 2.');
end
end