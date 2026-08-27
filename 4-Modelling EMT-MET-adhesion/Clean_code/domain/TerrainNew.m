function [zm, xm, ym] = TerrainNew_skata(p_resolution, p_peaks)
% p_resolution:: the final grid will have 2^p_resolution cells per direction
% p_peaks:: the landscape will exhibit 2^p_peaks per direction
    
    % the random peaks of the terrain
    zm=rand(2^p_peaks,2^p_peaks);
    
    s = size(zm,1); %(leave it as is for further refinement in the future)
    % Make the initial x and y meshes.
    lnsp=linspace(0,1,s);
    [xm, ym] = meshgrid(lnsp, lnsp);

    
    % Refine for n iterations.
    for k = p_peaks+1:p_resolution
        % Distance between neighboring components.
        d = 0.5^k * s;
        % Make new meshes.
        lnsp=linspace(0,1,2^k);
        [xm_new, ym_new] = meshgrid(lnsp, lnsp);
%         xm_new
%         size(xm), size(ym), size(zm), pause
        % Interpolate and add random variation to the new mesh.
        zm = interp2(xm, ym, zm, xm_new, ym_new) + d*rand(size(xm_new));%2^k/s);
%         size(zm)
        
        % What was old becomes new.
        xm = xm_new;
        ym = ym_new;
    end
end