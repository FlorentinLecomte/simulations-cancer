% small helper script to create an xc, yc grid and empty vectors to make
% the definition of the initial condition more comfortable

delta_x = (b-a)/nx;
x = a:delta_x:b;
xc = (x(1:end-1)+x(2:end))/2;
if spDim == 2
    delta_y = (d-c)/ny;
    y = c:delta_y:d;
    yc = (y(1:end-1)+y(2:end))/2;
elseif spDim == 1
    ny = 1;
end
