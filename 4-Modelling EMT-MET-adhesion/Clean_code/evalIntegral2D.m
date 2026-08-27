function [A11, A22] = evalIntegral2D(G, mask, methodId)
% [A11, A22] = evalIntegral2D(G, mask, methodId)
%
% This function evaluates the first component of the non-local term
% in 2D in the centre of all right cell faces (returned as matrix
% A11) and its second component in the centre of all upper cell
% faces  (returned as matrix A22) for a periodic (in both
% directions) grid function G of dimension mask.x2Dir.N2 times
% mask.x2Dir.N1 (alternatively, G can be a vector of length
% mask.x2Dir.N2*mask.x2Dir.N1, which is formed columenwise from
% that matrix). For the evaluation three different algorithms can
% be employed and are selected by methodId 
%   methodId == 1 ... summation based
%   methodId == 2 ... FFT-based I
%   methodId == 3 ... FFT-based II
% Usually, methodId == 3 is the most efficient.
% The non-local term is encapsulated in the data structure mask,
% which must contain two subfields x1Dir and x2Dir corresponding to
% the right and upper cell faces, respectively. Using function
% setupIntegralRule2D, this structure is set up as
%   mask.x1Dir = setupIntegralRule2D('x1',...);
%   mask.x2Dir = setupIntegralRule2D('x2',...);
% These calls also prepare data within the returned structures
% which is required for the FFT-based evaluation schemes.
%
% Alf Gerisch (alf.gerisch@mathematik.uni-halle.de)
% Version 1.0 of January 15, 2007
% Version 1.1 of Sep 808, 2008
%  * allowed for G to be a matrix and not as before just a vector.
%

if (size(G) ~= [mask.x2Dir.N2, mask.x2Dir.N1])
  % supplied G must be a vector, so try to reshape to a matrix of
  % appropriate size
  G = reshape(G, mask.x2Dir.N2, mask.x2Dir.N1);
end

switch methodId
 case 1  
  % Non-FFT Variant
  switch mask.BCs
   case 'pp'
    % extend for easy computation of integral (periodic BCs)
    km = max(mask.x1Dir.km,mask.x2Dir.km);
    kp = max(mask.x1Dir.kp,mask.x2Dir.kp);
    lm = max(mask.x1Dir.lm,mask.x2Dir.lm);
    lp = max(mask.x1Dir.lp,mask.x2Dir.lp);
    Ge = [
	G(end-km+1:end,end-lm+1:end) G(end-km+1:end,:) G(end-km+1:end,1:lp)
	G(        :   ,end-lm+1:end) G                 G(        :   ,1:lp)
	G(       1:kp ,end-lm+1:end) G(        1:kp,:) G(       1:kp ,1:lp)
	 ];
    % evaluate integral
    A11=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
    A22=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
    yindices1 = km +[-mask.x1Dir.km:mask.x1Dir.kp];
    yindices2 = km +[-mask.x2Dir.km:mask.x2Dir.kp];
    xindices1 = lm +[-mask.x1Dir.lm:mask.x1Dir.lp];
    xindices2 = lm +[-mask.x2Dir.lm:mask.x2Dir.lp];
    for i=1:mask.x2Dir.N1
      for j=1:mask.x2Dir.N2
	tmp = mask.x1Dir.Wx1.*Ge(j+yindices1,i+xindices1);
	A11(j,i)= sum(tmp(:));
	tmp = mask.x2Dir.Wx2.*Ge(j+yindices2,i+xindices2);
	A22(j,i)= sum(tmp(:));
      end
    end
   case {'pppp', 'zzzz'}
    % for A11
    km = mask.x1Dir.km;
    kp = mask.x1Dir.kp;
    lm = mask.x1Dir.lm;
    lp = mask.x1Dir.lp;
    switch mask.BCs
     case 'pppp'
      Gext = [
	  G(end-km+1:end,end-lm:end) G(end-km+1:end,:) G(end-km+1:end,1:lp)
	  G(        :   ,end-lm:end) G                 G(        :   ,1:lp)
	  G(       1:kp ,end-lm:end) G(        1:kp,:) G(       1:kp ,1:lp)
	     ];
     case 'zzzz'
      Gext = zeros(mask.x1Dir.Gext_z, mask.x1Dir.Gext_s);
      Gext((km+1):(km+mask.x1Dir.N2), (lm+2):(lm+1+mask.x1Dir.N1)) = G;
     otherwise
      error('not supported.')
    end
    %[size(Gext) mask.x1Dir.Gext_z mask.x1Dir.Gext_s]
    A11=nan(mask.x1Dir.A_z, mask.x1Dir.A_s);
    yindices = km+[-km:kp];
    xindices = lm+1+[-lm:lp]-1;
    for i=1:mask.x1Dir.A_s
      for j=1:mask.x1Dir.A_z
%	[i,j]
	tmp = mask.x1Dir.Wx1.*Gext(j+yindices,i+xindices);
	A11(j,i)= sum(tmp(:));
      end
    end
    
    % for A22
    km = mask.x2Dir.km;
    kp = mask.x2Dir.kp;
    lm = mask.x2Dir.lm;
    lp = mask.x2Dir.lp;
    switch mask.BCs
     case 'pppp'
      Gext = [
	  G(end-km:end,end-lm+1:end) G(end-km:end,:) G(end-km:end,1:lp)
	  G(      :   ,end-lm+1:end) G               G(      :   ,1:lp)
	  G(     1:kp ,end-lm+1:end) G(      1:kp,:) G(     1:kp ,1:lp)
	     ];
     case 'zzzz'
      Gext = zeros(mask.x2Dir.Gext_z, mask.x2Dir.Gext_s);
      Gext((km+2):(km+1+mask.x2Dir.N2), (lm+1):(lm+mask.x2Dir.N1)) = G;
     otherwise
      error('not supported.')
    end

    %[size(Gext) mask.x2Dir.Gext_z mask.x2Dir.Gext_s]
    A22=nan(mask.x2Dir.A_z, mask.x2Dir.A_s);
    yindices = km+1+[-km:kp]-1;
    xindices = lm+[-lm:lp];
    for i=1:mask.x2Dir.A_s
      for j=1:mask.x2Dir.A_z
	%[i,j],[size(Gext) mask.x2Dir.Gext_z mask.x2Dir.Gext_s]
	tmp = mask.x2Dir.Wx2.*Gext(j+yindices,i+xindices);
	A22(j,i)= sum(tmp(:));
      end
    end

   otherwise
    error('evalIntegralRule2D:: unsupported value of BCs.');
  end

 case 2  
  % FFT Variant I:  
  switch mask.BCs
   case 'pp'
    % compute fft of G
    gFFT = fft(G);
    % construct A11 and A22 matrices
    A11=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
    A22=zeros(mask.x2Dir.N2, mask.x2Dir.N1);
    for i=1:mask.x2Dir.N1
      A11(:,i) = sum(gFFT.*mask.x1Dir.Vx1FFT.viFFT(:,[i:-1:1 end:-1:i+1]), 2);
      A22(:,i) = sum(gFFT.*mask.x2Dir.Vx2FFT.viFFT(:,[i:-1:1 end:-1:i+1]), 2);
    end
    A11 = ifft(A11);
    A22 = ifft(A22);
   otherwise
    error('evalIntegralRule2D:: unsupported value of BCs.');
  end
  
 case 3
  switch mask.BCs
   case 'pp'
    gFFT2 = fft2(G);
    A11 = ifft2(mask.x1Dir.Vx1FFT.viFFT2 .* gFFT2);
    A22 = ifft2(mask.x2Dir.Vx2FFT.viFFT2 .* gFFT2);
   case {'pppp', 'zzzz', 'szpp'}
    % for A11
    km = mask.x1Dir.km;
    kp = mask.x1Dir.kp;
    lm = mask.x1Dir.lm;
    lp = mask.x1Dir.lp;
    switch mask.BCs
     case 'pppp'
      Gext = [
	  G(end-km+1:end,end-lm:end) G(end-km+1:end,:) G(end-km+1:end,1:lp)
	  G(        :   ,end-lm:end) G                 G(        :   ,1:lp)
	  G(       1:kp ,end-lm:end) G(        1:kp,:) G(       1:kp ,1:lp)
	     ];
     case 'zzzz'
      Gext = zeros(mask.x1Dir.Gext_z, mask.x1Dir.Gext_s);
      Gext((km+1):(km+mask.x1Dir.N2), (lm+2):(lm+1+mask.x1Dir.N1)) = G;
     case 'szpp'
      % initialise all zeros
      Gext = zeros(mask.x1Dir.Gext_z, mask.x1Dir.Gext_s);
      % set G matrix in the middle
      Gext((km+1):(km+mask.x1Dir.N2), (lm+2):(lm+1+mask.x1Dir.N1)) = G;
      % zero-flux on right boundary is ok already
      % symmetry at left boundary
      Gext((km+1):(end-kp),1:(lm+1)) = G(:,(lm+1):-1:1);
      % periodicity at top/bottom
      Gext(1:km,:) = Gext((end-kp-km+1):(end-kp),:); 
      Gext((end-kp+1):end,:) = Gext((km+1):(km+kp),:);
     otherwise
      error('not supported.')
    end
    %[size(Gext) mask.x1Dir.Gext_z mask.x1Dir.Gext_s]
    Gbig = zeros(mask.x1Dir.l,mask.x1Dir.L);
    Gbig(1:mask.x1Dir.Gext_z,1:mask.x1Dir.Gext_s) = Gext;
    A11big = ifft2(mask.x1Dir.vi_x1_fft2  .* fft2(Gbig));
    A11 = A11big(1:mask.x1Dir.A_z,1:mask.x1Dir.A_s);
    % for A22
    km = mask.x2Dir.km;
    kp = mask.x2Dir.kp;
    lm = mask.x2Dir.lm;
    lp = mask.x2Dir.lp;
    switch mask.BCs
     case 'pppp'
      Gext = [
	  G(end-km:end,end-lm+1:end) G(end-km:end,:) G(end-km:end,1:lp)
	  G(      :   ,end-lm+1:end) G               G(      :   ,1:lp)
	  G(     1:kp ,end-lm+1:end) G(      1:kp,:) G(     1:kp ,1:lp)
	     ];
     case 'zzzz'
      Gext = zeros(mask.x2Dir.Gext_z, mask.x2Dir.Gext_s);
      Gext((km+2):(km+1+mask.x2Dir.N2), (lm+1):(lm+mask.x2Dir.N1)) = G;
     case 'szpp'
      % initialise all zeros
      Gext = zeros(mask.x2Dir.Gext_z, mask.x2Dir.Gext_s);
      % set G matrix in the middle
      Gext((km+2):(km+1+mask.x2Dir.N2), (lm+1):(lm+mask.x2Dir.N1)) = G;
      % zero-flux on right boundary is ok already
      % symmetry at left boundary
      Gext((km+2):(end-kp),1:lm) = G(:,lm:-1:1);
      % periodicity at top/bottom
      Gext(1:(km+1),:) = Gext((end-kp-km):(end-kp),:); 
      Gext((end-kp+1):end,:) = Gext((km+2):(km+kp+1),:);
      otherwise
      error('not supported.')
    end
    %[size(Gext) mask.x2Dir.Gext_z mask.x2Dir.Gext_s]
    Gbig = zeros(mask.x2Dir.l,mask.x2Dir.L);
    Gbig(1:mask.x2Dir.Gext_z,1:mask.x2Dir.Gext_s) = Gext;
    A22big = ifft2(mask.x2Dir.vi_x2_fft2  .* fft2(Gbig));
    A22    = A22big(1:mask.x2Dir.A_z,1:mask.x2Dir.A_s);
   otherwise
    error('evalIntegralRule2D:: unsupported value of BCs.');
  end
 otherwise
  error('evalIntegral2D::Matrix-Vector product methodId non-existent.');
end

return; % end of function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



  

