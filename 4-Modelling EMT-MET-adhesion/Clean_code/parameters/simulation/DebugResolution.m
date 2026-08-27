% File for numerical parameters
nx = 2000;				ny = 50;    % space resolution
CFL = 0.49;                         % Courant number
max_delta_t = 0.1;
scIndex = 1;                        % IMEX scheme (see GetScheme)

theta = 2;                          % theta for the minmod limiter

krylov = 0;                         % 1 for bicgstab instead of LR decomp. for linear systems
checkPos = 1;                       % 1 to cancel if negative values occur
nos = fix(100);                     % number of time instances that will be saved (must be integer)
fixedSeed = 0;