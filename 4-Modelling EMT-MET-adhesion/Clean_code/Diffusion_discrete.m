function [discrete_diff,exists] = Diffusion_discrete(S, C, modP, numP, domP, submodel, syseq_or_alevel)

% function for the implicit part of the IMEX that solves
% (Id - d * Impl) K = temp. The Impl term of the system is the diffusion
% part.

% cc = sum(C,3);
% gf = S(:, :, 1);
% ecm = S(:, :, 2);



% unwrap some things
nx = numP.Nx;           ny = numP.Ny;
dx = numP.dx;           dy = numP.dy;
% pr = nx*ny;
L_diff = numP.Ld;       %Id = speye(pr);
alevels = numP.alevels;

Dc_M = modP.Dc_M;
Dc_sl = modP.Dc_sl;
Dc_a0 = modP.Dc_a0;

Dgf = modP.Dgf;


% do not return any diffusion (as the equation might be dummy)
% unless otherwise stated
exists = 0; 
discrete_diff = NaN; 

switch submodel
    case 'cancer' %cancer cell subsystem (linear wrt space; nonlinear wrt alpha)
        cc = C(:,:,syseq_or_alevel);
        level = alevels(syseq_or_alevel);
        %f = Dc_L./( 1 + exp(-Dc_sl*(numP.alevels(syseq_or_alevel)-Dc_a0)) );
        %non-linear diffusivity depending on the level of a
%         Dc = cc.^2;%cc.^2./(0.5 + cc.^2); % non linearity wrt c
%         Dc = cc.^2.*Dc_M./( 1 + exp(-Dc_sl*(numP.alevels(syseq_or_alevel)-Dc_a0)) ); %dependence on a
%         Dc = Dc_M*cc.^2.*exp(-5*(alevels(syseq_or_alevel)-1).^2 );
%         Dc = Dc_M*cc.^2./( 1+exp(-Dc_sl*(level-Dc_a0)) );
        Dc = Dc_M*cc.^2./( 1+exp(-Dc_sl*(level-Dc_a0)) );
%         1./( 1+exp(-Dc_sl*(alevels-Dc_a0)) )
%         pause
%         Dc = Dc_M;
        if isscalar(Dc)
            discrete_diff = Dc*L_diff;
        else
            discrete_diff = GiveDiffusionMatrixFV([nx,ny],[dx,dy],2, Dc);
        end
        exists = 1;
    case 'system' %non-cancerous subsystem
        if syseq_or_alevel == 1             % GF %linear diffusion
            discrete_diff = Dgf*L_diff;  
            exists = 1;
        %     elseif syseq_or_alevel == 2       % ECM % no diffusion
        %         discrete_diff = 0; 
        end
end % for switch


end
