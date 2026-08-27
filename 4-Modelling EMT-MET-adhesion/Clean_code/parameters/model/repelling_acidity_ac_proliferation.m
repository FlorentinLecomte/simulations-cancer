run('fullModel.m');

xi_2 = 0;                      % cross -diffusion 
xi_3 = 0.1;%0.5;                      % new: accidity taxis coefficient

% acid
Dac = 0.07;                    % diffusion
rho = 0.55;                        % production
eta = 0.05;                         % decay rare

delta = 0.2;          

acidityProlif = 1;
mu(2) = 0.1; % corresponds to mu0 in the paper