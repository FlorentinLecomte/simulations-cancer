% File for numerical parameters
nx = 5000;				ny = 1;   % space resolution
CFL = 0.49;                         % Courant number
max_delta_t = 0.1;
scIndex = 3;                        % IMEX scheme (see GetScheme)

theta = 2;                          % theta for the minmod limiter

krylov = 0;                         % 1 for bicgstab instead of LR decomp. for linear systems
checkPos = 0;                       % 1 to cancel if negative values occur
nos = fix(T*4.5);                            % number of time instances that will be saved (must be integer)