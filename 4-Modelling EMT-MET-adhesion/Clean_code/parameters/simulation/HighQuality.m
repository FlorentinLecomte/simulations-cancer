% File for numerical parameters
nx = 256;				ny = 256;   % space resolution
CFL = 0.49;                         % Courant number
max_delta_t = 0.1;
scIndex = 3;                        % IMEX scheme (see GetScheme)

theta = 2;                          % theta for the minmod limiter

krylov = 1;                         % 1 for bicgstab instead of LR decomp. for linear systems
checkPos = 0;                       % 1 to cancel if negative values occur
nos = 100;                            % number of time instances that will be saved (must be integer)