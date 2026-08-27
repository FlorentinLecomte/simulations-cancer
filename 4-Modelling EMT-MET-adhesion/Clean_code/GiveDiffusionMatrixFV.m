function L = GiveDiffusionMatrixFV(N, h, sDim, f)
% provide Diffusion matrix for uniform finite volumes in spDim = 1 or 2
% its a discrete version of the operator div(f div( . ))
% in the 2D case we assume h = [dx, dy] and N = [Nx, Ny].
% Homogeneous Neumann boundary conditions are assumed.

nx = N(1);
if nargin == 3 %if you don't have f, i.e. if you only have linear diffusion
    delta_x = h(1);
    
    
    ex = ones(nx,1);
    
    mx = [-1; -2*ones(nx-2,1); -1];
    Dxx = spdiags([ex mx ex], [-1 0 1], nx, nx); %1D discrete Laplacian in the x-direction ;
    
    switch sDim
        case 2
            delta_y = h(2);
            ny = N(2);
            ey = ones(ny,1);
            my = [-1; -2*ones(ny-2,1); -1];
            Dyy = spdiags([ey, my, ey], [-1 0 1], ny, ny); %1D discrete Laplacian in the y-direction ;
            L = kron(Dyy, speye(nx))*1/(delta_y)^2 + kron(speye(ny), Dxx)*1/(delta_x)^2 ;
        case 1
            L = Dxx * 1 / (delta_x)^2;
        otherwise
            error('Only 1D and 2D supported.')
    end
    
elseif nargin == 4
    switch sDim
        case 1
            fx_plus = (f(3:end)+f(2:end-1))/2;
            fx_minus = (f(2:end-1)+f(1:end-2))/2;
            mx = [-fx_minus(1); -(fx_minus+fx_plus); -fx_plus(end)];
            fx_plus2 = [0; fx_minus(1); fx_plus];
            fx_minus2 = [fx_minus; fx_plus(end); 0];
            L =  spdiags([fx_minus2 mx fx_plus2], [-1 0 1], N(1), N(1)) / h(1)^2;
        case 2
            %new and faster code
            ny = N(2);
            Lx = sparse(nx*ny,nx*ny);
            
            FXp = .5 * (f(3:end, :) + f(2:end-1, :));
            FXm = .5 * (f(2:end-1, :) + f(1:end-2, :));
            
            MX = [-FXm(1, :); -(FXp+FXm); -FXp(end, :)];
            
            FXp2 = [zeros(1, ny); FXm(1, :); FXp];
            FXm2 = [FXm; FXp(end, :); zeros(1, ny)];
            
            %LXX = spdiags([FXm2(:) MX(:) FXp2(:)], [-1 0 1], nx*ny, nx*ny);
            
            FYp = .5 * (f(:, 3:end) + f(:, 2:end-1));
            FYm = .5 * (f(:, 2:end-1) + f(:, 1:end-2));
            
            MY = [-FYm(:, 1), -(FYm+FYp), -FYp(:, end)];
            
            FYp2 = [zeros(nx, 1), FYm(:, 1), FYp];
            FYm2 = [FYm, FYp(:, end), zeros(nx, 1)];
            
            %LYY = spdiags([FYm2(:) MY(:) FYp2(:)], [-nx 0 nx], nx*ny, nx*ny);
            
            %L = LYY*1/h(2)^2 + LXX*1/h(1)^2;
            L = spdiags([FYm2(:)*h(2)^(-2) FXm2(:)*h(1)^(-2) ...
                MX(:)*h(1)^(-2)+MY(:)*h(2)^(-2) FXp2(:)*h(1)^(-2) ...
                FYp2(:)*h(2)^(-2)], [-nx -1 0 1 nx], nx*ny, nx*ny);
            
            % old and slow code:
%             for i=1:N(2)
%                 fx_plus = (f(3:end,i)+f(2:end-1,i))/2;
%                 fx_minus = (f(2:end-1,i)+f(1:end-2,i))/2;
%                 
%                 mx = [-fx_minus(1); -(fx_minus+fx_plus); -fx_plus(end)];
%                 
%                 fx_plus2 = [0;fx_minus(1);fx_plus];
%                 fx_minus2 = [fx_minus;fx_plus(end);0];
%                 
%                 Dxx = spdiags([fx_minus2 mx fx_plus2], [-1 0 1], nx, nx); %1D discrete Laplacian in the x-direction ;
%                 Lx((i-1)*nx+1:i*nx,(i-1)*nx+1:i*nx) = Dxx;
%             end
%             
%             Ly = sparse(nx*ny,nx*ny);
%             for j=1:N(1)
%                 fy_plus = (f(j,3:end)+f(j,2:end-1))/2;
%                 fy_minus = (f(j,2:end-1)+f(j,1:end-2))/2;
%                 my = [-fy_minus(1); -(fy_minus+fy_plus)'; -fy_plus(end)];
%                 fy_plus2 = [0;fy_minus(1);fy_plus'];
%                 fy_minus2 = [fy_minus';fy_plus(end);0];
%                 Dyy = spdiags([fy_minus2, my, fy_plus2], [-1 0 1], ny, ny); %1D discrete Laplacian in the y-direction ;
%                 e = zeros(nx,1);
%                 e(j) = 1;
%                 eMatrix = spdiags(e,0,nx,nx);
%                 Ly = Ly + kron(Dyy, eMatrix);
%             end
%             
%             L = Ly*1/h(2)^2 + Lx*1/h(1)^2;
        otherwise
            error('Only 1D and 2D supported.')
    end
%    fprintf('faster code deviation %f\n', norm(L(:)-Ln(:)))
end