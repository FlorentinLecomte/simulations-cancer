function nonlocal = setupIntegralRule2D_weights(dir, R, h, Omega, ...
					ruleId, NR, varargin)
% nonlocal = setupIntegralRule2D_weights(dir, R, h, Omega, ...
% 				 ruleId, NR, varargin)
%
% Assume that a N2 times N1 matrix G has entries G(j,i) \approx
% g(x1_i,x2_j) and can be extended periodically in both dimensions 
% whenever required. The grid points (x1_i,x2_j) are on a uniform
% spatial grid with grid spacing h in both dimensions. We want to
% approximate the non-local term (two components) 
% A\{g()\}(x1,x2):= 1/R\int_0^R r \int_0^{2\pi} 
%                   \eta(\theta) \Omega(r) g((x1,x2)+r\eta(\theta)) 
%                   d\theta dr 
% (where \eta(\theta):=(\cos(\theta);\sin(\theta))^T) 
%      for (x1,x2) = (x1_i+h/2, x2_j) in the case dir == 'x1' 
%       or (x1,x2) = (x1_i, x2_j+h/2) in the case dir == 'x2' 
% by evaluating, for the first component of the integral:
%   tmp = nonlocal.Wx1 .* G(j+[-nonlocal.km:nonlocal.kp],...
%                           i+[-nonlocal.lm:nonlocal.lp]); 
%   I1 = sum(tmp(:));
% and replacing nonlocal.Wx1 by nonlocal.Wx2 for the second
% component of the integral. Here, data stored in the structure
% nonlocal computed and returned by calling this function is
% employed. The elements of that structure are detailed below.
% In the above formula we assume that (j,i) is sufficiently in the
% interior of G; otherwise one must make use of the periodicity of
% G. A faster method of evaluating the integral for all (x1,x2) =
% (x1_i+h/2, x2_j) or (x1,x2) = (x1_i, x2_j+h/2) is via FFT, see
% evalIntegral2D.m for two FFT options. The data required there is
% also prepared here and returned in the structure nonlocal.
%
% The function \Omega(r) is prescribed by the argument Omega of this
% Matlab function. Omega either be an integer (1 or 2) or a
% function handle. In the case that Omega is an integer, then 
%  Omega==1  ... constant \Omega(r),
%  Omega==2  ... linearly decaying \Omega(r).
% In the case that Omega is a function handle, then Omega must be
% callable with arguments r\in(0,+\infty) and returns \Omega(r); for
% r>R the return value of Omega is always corrected to zero (make
% sure Omega does not return NaN in that case). 
%
% A piecewise bilinear reconstruction of the function g(x1,x2),
% based on the grid data stored in matrix G is used to derive the
% weight matrices Wx1 and Wx2 for the approximation of the
% non-local term. Note that the matrix G is not required for the
% computation of these weights. The non-local term is discretised
% by a composite integration rule (selected with argument ruleId).   
% An iterative scheme is used to compute the weights Wx1 and Wx2 to
% a prescribed accuracy varargin{1}. Starting with NR and 2*NR
% subdivisions in the radial direction in the composite integration
% rule, the resulting weights are compared elementwise and the 2*NR
% weights are accepted if the maximum absolute difference is less
% than the prescribed accuracy varargin{1}. If not, the weights
% corresponding to 2^2*NR subdivisions in the radial direction are
% computed and so on. A maximum of 10 doublings of the initial NR
% is executed and an error is issued if the required accuracy is
% not met. Typically a quadratic decay of the error is observed. To
% use a specific value NR for the computation of the weights
% without any iteration, set varargin{1} = +inf. This can be used
% if for instance the weights with dir = 'x1' have been computed to
% a certain accuracy and the weights for dir ='x2' are supposed to
% be computed with the same resolution. 
%
% Input arguments
% ---------------
% dir       ... 'x1' or 'x2'
% R         ... sensitivity radius
% h         ... grid width in spatial grid
% Omega     ... function handle returning Omega(r) for r>=0 or integer.
% ruleId    ... = 1 --  trapezoidal rule, piecewise bilinear reconstruction
%               = 2 --  midpoint rule,    piecewise bilinear reconstruction
% NR        ... number of radial discretisation points (initial value)
% varargin{1} ... tolerance for computation of weights, default = h;
%                 if == +inf then do not iterate but just use NR
% varargin{2} ... use two different methods to compute matrices
%                 (true/false) default=false.
%
% Output data structure
% ---------------------
% The structure returned has the following fields:
% dir, R, h, ruleId ... copies of the respective input
% OmegaId           ... value of Omega if it is an integer
%                       and not a function handle
% NR0                       ... copy of the input NR
% weighttol                 ... value used as prescribed accuracy
%                               of weights
% NR                        ... NR used for the final computation
%                               of the weights 
% M                         ... number of grid cells of size h to
%                               cover sensitivity radius R
% lm, lp, km, kp            ... the weight matrices Wx1, Wx2 are 
%                               of size km+kp+1 times lm+lp+1
% Wx1, Wx2                  ... weight matrices for the first and
%                               second component of the non-local
%                               term
%
%
% Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
% Version of 1.1 of Sep 30, 2008
%   * changed input argument Omega which can now be a
%   function handle or an integer (1 or 2) to make it consistent
%   with the 1D version of the setup function. Details are
%   documented above.
% Version of 1.0 of Jan 15, 2008 
%

% List of subfunctions in this file
% ==========================================
% 
% function [l,k,delta1,delta2] = getBilinearInterpParams(evalAt, h, ci, cj, lp, kp)
% function w = getBilinearInterpWeights(evalAt, lm, lp, km, kp, h)
%

format short e

% deal with varargin values
weighttol = h;
if (nargin > 8)
  if ~isempty(varargin{1})
    weighttol = varargin{1};
  else
    disp('setupIntegralRule2D::using default for weighttol.')
  end
end

comparison = false;
if (nargin > 9)
  if ~isempty(varargin{2})
    comparison = varargin{2};
  else
    disp('setupIntegralRule2D::using default for comparison.')
  end
end

switch ruleId
 case 1
  disp(['setupIntegralRule2D::trapezoidal rule, ' ...
	'piecewise bilinear reconstruction']);  
 case 2
  disp(['setupIntegralRule2D::midpoint rule, ' ...
	'piecewise bilinear reconstruction']);  
 otherwise
  error('setupIntegralRule2D::unknown ruleId.');
end

nonlocal.dir = dir;
nonlocal.R = R;
nonlocal.h = h;
nonlocal.ruleId  = ruleId;
nonlocal.NR0 = NR;
nonlocal.weighttol = weighttol;
nonlocal.NR = NR;
nonlocal.M = ceil(R/h);

% check Omega
if isnumeric(Omega)
  switch Omega
   case 1
    OmegaId = Omega;
    nonlocal.OmegaId = Omega;
    OmegaH = @(r)((1/(pi*(nonlocal.R)^2))*(abs(r)<=nonlocal.R));
   case 2
    OmegaId = Omega;
    nonlocal.OmegaId = Omega;
    OmegaH = @(r)((3/(pi*(nonlocal.R)^2))*(1-r/nonlocal.R) ...
		  *(abs(r)<=nonlocal.R));
   otherwise
    error('setupIntegralRule2D::Unsupported value of OmegaId.');
  end
else
  if isa(Omega, 'function_handle') 
    OmegaH = @(r)(Omega(r)*(abs(r)<=nonlocal.R));
  else
    error(['setupIntegralRule2D::supplied Omega is of unsupported' ...
	   ' type.']); 
  end
end
  
% check value of M and define alpha
alpha = nonlocal.R-(nonlocal.M-1)*nonlocal.h;
if (alpha <= 0)
  if (alpha < -eps)
    disp(alpha)
    error('setupIntegralRule2D::computed value of alpha is <= 0.');
  else
    disp(['setupIntegralRule2D::adjusted value of alpha from slightly' ...
	  ' negative to 1e-20.']);
    alpha=1e-20;
  end
end
if (alpha > nonlocal.h)
  if (alpha-eps > nonlocal.h)
    disp(alpha-nonlocal.h);
    error('setupIntegralRule2D::computed value of alpha is > h.');
  else
    disp(['setupIntegralRule2D::adjusted value of alpha from just above' ... 
	  ' h to h.']);
    alpha = nonlocal.h;
  end
end

nonlocal.lm= nonlocal.M;
nonlocal.lp= nonlocal.M;
nonlocal.km= nonlocal.M;
nonlocal.kp= nonlocal.M;
switch nonlocal.dir
 case 'x1'
  nonlocal.lp = nonlocal.M+1;
 case 'x2'
  nonlocal.kp = nonlocal.M+1;
 otherwise
  error('setupIntegralRule2D::dir must be x1 or x2')
end

% the index corresponding to (x1_i,x2_j) in matrices Wx1 and Wx2 is (cj,ci)
cj = nonlocal.km+1;
ci = nonlocal.lm+1;

maxiter = 10;
itercount = 0;
converged = false;

while (~converged)
  itercount = itercount + 1;
  % initialise weight matrices
  Wx1 = zeros(nonlocal.km+nonlocal.kp+1, nonlocal.lm+nonlocal.lp+1);
  Wx2 = Wx1;
  if (comparison)
    Wx1b=Wx1;
    Wx2b=Wx2;
  end

  % the grid spacing along the radial direction
  DR = nonlocal.R/nonlocal.NR;
  for m=1:nonlocal.NR % loop over nR circle
    switch ruleId           % the radius of the current circle
     case 1
      rm = (m*nonlocal.R)/nonlocal.NR;       % for trapezoidal rule
     case 2
      rm = ((m-0.5)*nonlocal.R)/nonlocal.NR; % for midpoint rule
     otherwise
      error('never get here!')
    end
    % facRm = rm*Omega(rm)/R; % do not use DR
    facRm = rm*OmegaH(rm)/nonlocal.NR; % the radial factor in the integration formula 
    if (m==nonlocal.NR)        % adjusted if on the outermost radius 
      if (ruleId == 2)         % (trapezoidal rule) 
	facRm = facRm * 0.5;
      end
    end
    % compute the integral over the circle with radius rm
    curNTheta=max(10, round(2*pi*rm/DR));
    curDTheta=2*pi/curNTheta;
    for n = 1:curNTheta
      switch ruleId              % the current theta value
       case 1
	thetan = curDTheta*n;       % for trapezoidal rule
       case 2
	thetan = curDTheta*(n-0.5); % for midpoint rule
       otherwise
	error('never get here!')
      end
      etan   = [cos(thetan); sin(thetan)]; % corresponding unit outer normal
      
      switch dir
       case 'x1'
	% g needs to be evaluated at (x1_i,x2_j)+(h/2,0)+rm*etan;
	evalAt = rm*etan+[h/2;0];          % relative to (x1_i,x2_j)
       case 'x2'
	% g needs to be evaluated at (x1_i,x2_j)+(0,h/2)+rm*etan;
	evalAt = rm*etan+[0;h/2];          % relative to (x1_i,x2_j)
       otherwise
	error('setupIntegralRule2D::dir must be x1 or x2')
      end
      
      % evaluate integrand in that point by using bilinear interpolation
      
      % (i) get indices of lower left corner of 'dual' grid cell (relative to
      % (cj,ci) as (k,l)) containing evalAt and corresponding fractions
      % delta1 (x1-dir.) and delta2 (x2-dir.)  
      [l,k,delta1,delta2] = getBilinearInterpParams(evalAt, h, ci, cj, ...
				nonlocal.lp, nonlocal.kp); 
      % now evalAt(1) = h*(l+delta1))
      % and evalAt(2) = h*(k+delta2))
      %disp([l,k,delta1,delta2])
      
      gamma = facRm*etan*curDTheta;
      % (ii) update matrix entries accordingly
      % lower left corner
      Wx1(cj+k,ci+l) = Wx1(cj+k,ci+l) + gamma(1)*(1-delta1)*(1-delta2);
      Wx2(cj+k,ci+l) = Wx2(cj+k,ci+l) + gamma(2)*(1-delta1)*(1-delta2);
      % lower right corner
      Wx1(cj+k,ci+l+1) = Wx1(cj+k,ci+l+1) + gamma(1)*delta1*(1-delta2);
      Wx2(cj+k,ci+l+1) = Wx2(cj+k,ci+l+1) + gamma(2)*delta1*(1-delta2);
      % upper left corner
      Wx1(cj+k+1,ci+l) = Wx1(cj+k+1,ci+l) + gamma(1)*(1-delta1)*delta2;
      Wx2(cj+k+1,ci+l) = Wx2(cj+k+1,ci+l) + gamma(2)*(1-delta1)*delta2;
      % upper right corner
      Wx1(cj+k+1,ci+l+1) = Wx1(cj+k+1,ci+l+1) + gamma(1)*delta1*delta2;
      Wx2(cj+k+1,ci+l+1) = Wx2(cj+k+1,ci+l+1) + gamma(2)*delta1*delta2;
      
      if (comparison)
	wm =getBilinearInterpWeights(evalAt, nonlocal.lm, nonlocal.lp, ...
					     nonlocal.km, nonlocal.kp, h);
	Wx1b = Wx1b + gamma(1)*wm;
	Wx2b = Wx2b + gamma(2)*wm;
      end
    end % of theta loop
  end % of r loop
  
  if (comparison)
    wdiff1 = norm(Wx1(:)-Wx1b(:), inf);
    wdiff2 = norm(Wx2(:)-Wx2b(:), inf);
    if (wdiff1>1e-14)
      wdiff1
      error('error computing first weight matrix')
    else
      disp(['Wx1 computed accurately (max diff = ' num2str(wdiff1) ').']);
    end
    if (wdiff2>1e-14)
      wdiff2
      error('error computing second weight matrix')
    else
      disp(['Wx2 computed accurately (max diff = ' num2str(wdiff2) ').']);
    end
  end

  if (itercount  == 1) % run another time
    if (nonlocal.weighttol == +inf) % no iteration!
      converged = true;
    else
      disp(['iteration     NR      maxdiff']);
      disp(['---------------------------------']);
      fprintf(1,'%9i%7i           ---\n', itercount, nonlocal.NR);
      Wx1old = Wx1;
      Wx2old = Wx2;
      nonlocal.NR = 2*nonlocal.NR;
      if (comparison)
	Wx1bold = Wx1b;
	Wx2bold = Wx2b;
      end
    end
  else % make convergence test
    wdiff1 =  Wx1old - Wx1;
    wdiff2 =  Wx2old - Wx2;
    wdiff = norm([wdiff1(:); wdiff2(:)], inf);
      fprintf(1,'%9i%7i      %10.2e\n', itercount, nonlocal.NR, wdiff);
    if (wdiff < nonlocal.weighttol)
      disp('...converged')
      converged = true;
    else
      if (itercount >= maxiter) 
	disp('...Maximum number of iterations reached without convergence.')
	error('reached maxiter, no convergence');
      else
        Wx1old = Wx1;
	Wx2old = Wx2;
	nonlocal.NR = 2*nonlocal.NR;
	if (comparison)
	  Wx1bold = Wx1b;
	  Wx2bold = Wx2b;
	end
      end
    end
  end
  
end % while (~converged)
  
% Iteration converged successfully
nonlocal.Wx1  = Wx1;
nonlocal.Wx2  = Wx2;
if (comparison)
  nonlocal.Wx1b  = Wx1b;
  nonlocal.Wx2b  = Wx2b;
end

return % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%%      S U B F U N C T I O N S
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [l,k,delta1,delta2] = getBilinearInterpParams(evalAt, h, ...
						       ci, cj, lp, kp)
% [l,k,delta1,delta2] = getBilinearInterpParams(evalAt, h, ...
%						ci, cj, lp, kp)
%
% Compute [l,k,delta1,delta2] such that evalAt(1) = h*(l+delta1))
% and evalAt(2) = h*(k+delta2)) subject to the side conditions
%   * 0 <= delta1,delta2 <= 1
%   * ci+l > 0 and cj+k > 0
%   * l < lp and k < kp
%

% (i) compute the lower left index of evalAt relative to (cj,ci),
% i.e. (cj,ci) corresponds to (0,0), and the fractions delta1 and delta2
% are in [0,1]
l = floor(evalAt(1)/h);
delta1 = (evalAt(1) - h*l)/h;
k = floor(evalAt(2)/h);
delta2 = (evalAt(2) - h*k)/h;

% (i) ensure 0<=delta1,delta2<=1
treshfac = 100;
if (delta1 < 0)
  if (delta1 > -treshfac*eps)
    delta1 = 0;
    disp ('setting delta1 to 0');
  else
    delta1
    error('bad value of delta1')
  end
end
if (delta1 > 1)
  if (delta1-1 < treshfac*eps)
    delta1 = 1;
    disp ('setting delta1 to 1');
  else
    delta1
    error('bad value of delta1')
  end
end
if (delta2 < 0)
  if (delta2 > -treshfac*eps)
    delta2 = 0;
    disp ('setting delta2 to 0');
  else
    delta2
    error('bad value of delta2')
  end
end
if (delta2 > 1)
  if (delta2-1 < treshfac*eps)
    delta2 = 1;
    disp ('setting delta2 to 1');
  else
    delta2
    error('bad value of delta2')
  end
end

% (ii) ensure limits on (k,l)
if (ci+l<=0)
  error('l too small');
end
if (cj+k<=0)
  error('k too small');
end
if (l>lp)
  error('l too large');
end
if (l==lp)
  if (delta1 < treshfac*eps)
    l = l-1;
    delta1 = 1;
    disp('reduced l');
  else
    error('out of bounds l');
  end
end
if (k>kp)
  error('k too large');
end
if (k==kp)
  if (delta2 < treshfac*eps)
    k = k-1;
    delta2 = 1;
    disp('reduced k');
  else
    error('out of bounds k');
  end
end

%disp([lp,kp, l,k, delta1, delta2]);

% (iii) test computed values:
if (norm(evalAt(1)-h*(l+delta1))>1e-15)
  error('x error')
end
if (norm(evalAt(2)-h*(k+delta2))>1e-15)
  evalAt(2)-h*(k+delta2)
  error('y error')
end

return; % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function w = getBilinearInterpWeights(evalAt, lm, lp, km, kp, h)
% w = getBilinearInterpWeights(evalAt, lm, lp, km, kp, h)
%
% This function returns a matrix of interpolation weights for bilinear
% interpolation of data given on the uniform grid with coordinates
% x1 = h*[-lm:lp] and x2 = h*[-km:kp] at the SINGLE point evalAt
% within the rectangle defined by x1 and x2. 
%
% The data is assumed to be stored in a matrix with increasing row index
% correspoding increasing x2 coordinates and increasing column index
% corresponding to increasing x1 coordinates. The weight matrix w returned
% assumes the same ordering. The weight matrix has at most 4 non-zero
% entries but is returned as a full matrix.
%

% Setup the grid
x1 = h*[-lm:lp];
x2 = h*[-km:kp];
[XX1,XX2] = meshgrid(x1,x2);

% Initialise weight matrix to zero matrix
w = zeros(km+kp+1,lm+lp+1);

% This is a sparse zero matrix for the basis function
w0 = sparse(km+kp+1,lm+lp+1);

% loop over all basis functions
for j = 1:km+kp+1
  for i = 1:lm+lp+1
    % set up basis matrix
    Wji = w0; Wji(j,i) = 1;
    % write computed weight to weight matrix
    w(j,i) = interp2(XX1, XX2, Wji, evalAt(1), evalAt(2), '*linear', NaN);
  end
end

return
% End of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

