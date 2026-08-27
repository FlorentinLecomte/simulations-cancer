% File for numerical parameters
nx = 250;				ny = 250;   % space resolution
CFL = 0.49;                         % Courant number
max_delta_t = 0.1;
scIndex = 3;                        % IMEX scheme (see GetScheme)

theta = 2;                          % theta for the minmod limiter

krylov = 1;                         % 1 for bicgstab instead of LR decomp. for linear systems
checkPos = 0;                       % 1 to cancel if negative values occur
nos = fix(T*3);                            % number of time instances that will be saved (must be integer)