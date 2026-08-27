% File for model parameters
Dc = 1e-3;%1e-2;                  % diffusion
xi_1 = 0.4;%1;                      % haptotaxis
xi_2 = 0.1;%1;                      % cross -diffusion 
xi_3 = 0;                      % new: accidity taxis coefficient
kD = 1;

mu = [0, 0.1, 0.15];%[ 0, 1, 0.5 ];           % proliferation (cell, cell, matrix)
% the first proliferation is neglected

constEMT = 0; 
constMET = 0;

delta = 0.3;%5;                    % matrix degradation

gamma0 = 0.1;

degDiff = 0;
% for nonconstant EMT and MET 
R_0 = 1;%2;                          % total amount of receptors
b = 2;                              % main coeff in MET 
p = 2;                              % plays the role of q in the text
mu_y = 1.1;
mu_zeta = 1.1;
sig_y = 0.5;
sig_zeta = 0.3;
y_ref = 0.6;

% acid
Dac = 0*1e-1;                    % diffusion
rho = 0*0.4;                        % production
eta = 0*0.2;                         % decay rare


