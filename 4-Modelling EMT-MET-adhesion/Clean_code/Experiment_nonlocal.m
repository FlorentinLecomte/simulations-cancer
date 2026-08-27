% use this script to start/continue an Experiment
% the settings are loaded from the file with the name written
% in the setup variable in the settings folder. If the setup
% variable is not specified 'DefaultSetup' will be loaded

if ~exist('steps','var') && ~exist('finished','var')
%    close all;
    clearvars -except setup sim_cue cue_counter err_cue par_var ...
        par_study_res par_range par_value;
    steps = 1;
    finished = 0;
end


if finished
    error('Simulation finished already.');
end

dispstat('', 'init');

if steps == 1
    finished = 0;
    % load setup parameters
    if ~exist('setup', 'var')
        setup = 'DefaultSetup';
    end
    run(['settings', filesep, setup]);

    if (~exist(['parameters/model/', model_parameters,'.m'], 'file') || ...
            ~exist(['parameters/simulation/', simulation_parameters,'.m'], 'file') || ...
            ~exist(['domain/', initial_condition,'.m'], 'file'))
        error('Something is wrong with your setup file.');
    end

    % load, process, and wrap parameters
    run(['parameters/model/', model_parameters]);
    run(['parameters/simulation/', simulation_parameters]);
    if ~exist('fixedSeed', 'var')
        fixedSeed = 1;
    end

    if fixedSeed
        rng(44760, 'twister');
    end

    % INVOKE ICs from file
    run(['domain/', initial_condition]);


    dt = T/nos;

    delta_x = (b-a)/nx;

    if spDim == 1
        Omega = b -a;
        L_diff = GiveDiffusionMatrixFV(nx,delta_x,1);
    elseif spDim == 2
        Omega = (b-a)*(d-c);							% domain size
        delta_y = (d-c)/ny;
        L_diff = GiveDiffusionMatrixFV([nx,ny],[delta_x,delta_y],2);
    else
        error('invalid choice of space dimension')
    end

    t = 0:dt:T;

    WrapParameters;



    if ~exist('par_var', 'var')
        par_var.active = 0;
    end

    if par_var.active
        pvalue = getfield(modP, par_var.parameter);
        pvalue(par_var.index) = par_value;
        modP = setfield(modP, par_var.parameter, pvalue);
    end

    % in S we keep the components of the solutions, i.e.
    % data matrix:
    % (nx lines) x (ny columns), so, x is along the vertical direction
    % i.e the first columns corr. to y=min and the last to y=max
    % and the first row corr. to x=min and the last to x=max
    % it has also 4 components of the system
    % and nos saved positions in time
    if spDim == 1
        %%%%%%%%%%% adjust to your system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S = zeros(nx, modP.comps, nos+1);
%         S(:, 1, 1) = CC_IC;
        S(:, 1, 1) = GF_IC;
        S(:, 2, 1) = ECM_IC;
        if modP.comps > 2
            S(:, 3, 1) = FIB_IC;
        end
        C = zeros(nx, alpha_num_levels, nos+1);
        C(:, :, 1) = CC_IC_all;
    elseif spDim == 2
        %%%%%%%%%%% adjust to your system %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        S = zeros(nx, ny, modP.comps, nos+1);
%         S(:, :, 1, 1) = CC_IC;
        S(:, :, 1, 1) = GF_IC;
        S(:, :, 2, 1) = ECM_IC;
        if modP.comps> 2
            S(:, :, 3, 1) = FIB_IC;
        end
        C = zeros(nx, ny, alpha_num_levels, nos+1);
        C(:, :, :, 1) = CC_IC_all;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

    % to make the method start from initial time
    steps = 2;
    dispstat('Simulation started', 'keepthis', 'timestamp');
else
    dispstat('Simulation resumed', 'keepthis', 'timestamp');
end

fileName = ['results', filesep, setup, 'Data.mat'];

if exist(fileName, 'file')
    prev = load(fileName, 'modP', 'numP', 'domP', 'S', 't');
    if strcmp(DataHash(modP),DataHash(prev.modP)) && strcmp(DataHash(numP), DataHash(prev.numP)) && ...
            strcmp(DataHash(domP), DataHash(prev.domP)) && strcmp(DataHash(S(:, :, :, 1)), DataHash(prev.S(:, :, :, 1))) && ...
            strcmp(DataHash(t), DataHash(prev.t))
        fprintf('Results have been computed before. Loading them now...\n');
        load(fileName);
        clearvars -except S t modP numP domP finished sim_cue cue_counter err_cue ...
            setup par_study_res par_var par_range therapy;
        return;
    else
        fprintf('Results of a modified version of this experiment exist.\n')
        prevFile = [fileName(1:end-4), 'Moved', datestr(now, 'dd-HHMM'), '.mat'];
        fprintf('Moving them to %s...\n', prevFile)
        movefile(fileName, prevFile);
    end
end

%%Non local mask%%
FNonLocal.N1 = numP.Nx;       % N1 value
FNonLocal.N2 = numP.Ny;       % N2 value
FNonLocal.h  = numP.dx%*numP.dy;      % grid spacing in both directions
FNonLocal.RuleId = 1;     % select integration rule
                           % 1 = trapezoidal rule (good choice)
FNonLocal.NR = 100;       % initial number of subdivisions of integral
FNonLocal.weightTol =5e-4;% desired tolerance for integration weights
FNonLocal.R  = 0.3;       % sensing radius
% function handle for Omega(r)
FNonLocal.Omega = @(r)((1-r/FNonLocal.R).*(r<=FNonLocal.R));
FNonLocal.BCs = 'zzzz';   % for zero-flux BCs

%FNonLocal

% now we can setup the integration weights (mask) for the nonlocal term
% on the vertical grid cell faces (x1-direction)
% (i) first for periodic boundary conditions
clear FNonLocal.mask % just in case there is something
FNonLocal.mask.BCs = FNonLocal.BCs;
nonlocal = setupIntegralRule2D_weights('x1', ...
    FNonLocal.R, ...
    FNonLocal.h, ...
    FNonLocal.Omega, ...
    FNonLocal.RuleId, ...
    FNonLocal.NR, ...
    FNonLocal.weightTol);
% (ii) and second correct for zero flux BCs and store final mask
FNonLocal.mask.x1Dir = setupIntegralRule2D_BCs(nonlocal, ...
    FNonLocal.N1, ...
    FNonLocal.N2, ...
    FNonLocal.BCs);

% and now we do the same for the nonlocal term on the horizontal
% grid cell faces (x2-direction); in this direction we save some
% time and do not iterate to achieve weightTol for the weights
% but just use the final NR from the x1 direction above
nonlocal = setupIntegralRule2D_weights('x2', ...
    FNonLocal.R, ...
    FNonLocal.h, ...
    FNonLocal.Omega, ...
    FNonLocal.RuleId, ...
    FNonLocal.mask.x1Dir.NR, ...
    inf);
FNonLocal.mask.x2Dir = setupIntegralRule2D_BCs(nonlocal, ...
    FNonLocal.N1, ...
    FNonLocal.N2, ...
    FNonLocal.BCs);
fprintf('Precomputation of weights for nonlocal term evaluation completed.\n')


%time evolution
tic
while (steps-1)*dt <= T
    if spDim == 1
        S(:,:,steps) = Integrate_EMT_Model_1D(S(:,:,steps-1), ...
            modP, numP, domP, t(steps-1),(steps-1)*dt);

    elseif spDim == 2
        [S(:,:,:,steps), C(:,:,:,steps)] = Integrate_EMT_Model_2D_nonlocal(...
            S(:,:,:,steps-1), C(:,:,:,steps-1), ...
            modP, numP, domP, t(steps-1),(steps-1)*dt, FNonLocal);
    end
    dispstat(sprintf('##### Progress of computation: %5.2f %% #####', (steps-1)/nos * 100));
    steps = steps + 1;
end
dispstat('Computation finished ', 'keepthis', 'timestamp');
toc

%save the results
finished = 1;
%numP = rmfield(numP, 'Ld');
if ~ par_var.active
    save(fileName, 'S', 'C' , 't', 'modP', 'numP', 'domP', 'finished', ...
    'setup');
end
clearvars -except S C t modP numP domP finished sim_cue cue_counter err_cue ...
    setup par_study_res par_var par_range therapy;
