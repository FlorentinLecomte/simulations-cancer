% File for numerical parameters
nx = 128;				ny = 128;   % space resolution
CFL = 0.49;%0.39;                         % Courant number
max_delta_t = 0.1;
scIndex = 1;                        % IMEX scheme (see GetScheme)

theta = 2;                          % theta for the minmod limiter

krylov = 0;                         % 1 for bicgstab instead of LR decomp. for linear systems
checkPos = 0;                       % 1 to cancel if negative values occur
nos = 100;                            % number of time instances that will be saved (must be integer)