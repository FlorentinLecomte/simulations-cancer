% script to visualize the results 

% close all;
%scH = 10^-6.4;
F = figure(1);
clf;
if finished
    x = domP.a:numP.dx:domP.b;
    xc = (x(1:end-1)+x(2:end))/2;
    if domP.spDim == 2
        y = domP.c:numP.dy:domP.d;
        yc = (y(1:end-1)+y(2:end))/2;
    end
    steps = numP.nos+2;
end

%% 1D visualisation
if domP.spDim == 1
    for i=1:steps-1
        figure(1);
        plot(xc, S(:,:,i));
        title(['t=',num2str(t(i))]);
        legend(modP.compNames);
        pause(0.05);
        if i == 1
            pause;
        end
    end    
%% 2D visualisation
elseif domP.spDim == 2
    initial_step = 1;%steps-20;
    %for i=initial_step:steps-1
    for i = round(linspace(1, steps-1, 11))    
    clf        
    tiles = tiledlayout(3,3);
    title(tiles, ['t=', num2str(t(i))])

    %---cancer cells
    for level=[1, ceil(0.5*(1+numP.alpha_num_levels)), numP.alpha_num_levels]%[1,5,9]
        ax1 = nexttile;
        imagesc(xc, yc, ((numP.da*C(:,:,level,i).'))); 
        set(gca,'YDir','normal');
        axis square;
        colormap(ax1,(hot));
        %xlim([-4,4]); ylim([-4,4]);
        title(['CCs, a=',num2str(numP.alevels(level))]);
        colorbar; 
        % clim([1e-15,1e+0]); set(gca, 'ColorScale','log')
%         axis([0 5 -5 0]);
    end
    for level=[1, ceil(0.5*(1+numP.alpha_num_levels)), numP.alpha_num_levels]%[1,5,9]
        ax2 = nexttile;
        isolines = 20;%10.^[0:-1:-15];
        contour(xc, yc, numP.da*C(:,:,level,i).',isolines);
        %imagesc(xc, yc, C(:,:,level,i).'); 
        set(gca,'YDir','normal'); 
        axis square;
        %xlim([-4,4]); ylim([-4,4]);
        title(['CCs, a=',num2str(numP.alevels(level))]);
        colorbar; colormap(ax2,hot);
    end

        %---growth factors
        axgf = nexttile;
        % imagesc(xc, yc, S(:,:,1,i).');    % plot ths 2nd matrix component
        contour(xc, yc, S(:,:,1,i).',20);    % contour of the 2nd matrix component
        set(gca,'YDir','normal');
        axis square;
        %xlim([-4,4]); ylim([-4,4]);
        title('GFs');
        colorbar; colormap(axgf,flipud(summer))
        
        %---ECM
        axecm = nexttile;
        imagesc(xc, yc, S(:,:,2,i).');    % plot ths 3rd matrix component
        set(gca,'YDir','normal');
        axis square;
        %xlim([-4,4]); ylim([-4,4]);
        title('ECM');
        colorbar; colormap(axecm,'bone');

        axecm2 = nexttile;
        contour(xc, yc, S(:,:,2,i).',50);    % plot ths 3rd matrix component
        set(gca,'YDir','normal');
        axis square;
        %xlim([-4,4]); ylim([-4,4]);
        title('ECM');
        colorbar; colormap(axecm2,'bone') 
        
        drawnow;
        pause(0.005);
        %if i == 1
        %    pause;
        %end
        % saveas(F,['results\AFTER_REVIEW_saved_experiment_uniECM_PAPER\experiment_uniECM_',num2str(i),'.png'],'png')
        %if mod(i,50)==0
        %    saveas(F,['results\saved_figures\experiment_',num2str(i),'.png'],'png')
        %end
    end
    saveas(F,['results\saved_figures\experiment','.png'],'png')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if finished
    clearvars x xc y yc steps ssm F i;
end
