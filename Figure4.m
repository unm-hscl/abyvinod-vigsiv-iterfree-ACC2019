close all, clear,clc
load('Figure3a-0.2.mat')
        
    %% Plotting trajectories of each method: 

        plot_markersize = 10;
        plot_fontSize = 11;
        fig3 = figure(1);
        subplot(2,1,1)
        hold on
        h1 = plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        hold on
        
        h2 = plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        h3 = scatter(x0(1),x0(2),5*plot_markersize,'bo','filled');
        h4 = scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'md');
        h5 = scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'b^');
        h6 = scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'r*');
        h7 = scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    10*plot_markersize, 'ks');


       
        
        xlabel('Position ($z_1$)')
        ylabel('Position ($z_2$)')
%         title('\textbf{Trajectory}')
        [hleg, hobj, hout, mout] = legend([h1 h2 h3 h4 h5 h6 h7],{'Safe set',...
            'Target set','Initial position',...
            ['Piecewise affine' newline 'approach - QP'],...
            ['Iterative risk' newline 'allocation (IRA)'],...
            ['Piecewise affine' newline 'approach - MIQP'],...
            ['Particle control' newline sprintf('(PC), %i particles',...
            N)]},...
            'Location','EastOutside','FontSize',11);
                set(gca,'FontSize',11)
                hobj(10).Children.MarkerSize = 10;
                hobj(11).Children.MarkerSize = 10;
                hobj(12).Children.MarkerSize = 10;
                hobj(13).Children.MarkerSize = 10;
                hobj(14).Children.MarkerSize = 10;
        
        box on;
        axis tight
        axis([-1 0.2 -1 0.05])
        set(gca,'FontSize',plot_fontSize)
        set(fig3,'Units','centimeters');
        set(fig3,'Position',[0 0 15 10]);
        
        handaxes3 = axes('position', [0.175 0.75 0.125 0.125]);
        hold on
        grid on
        plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        scatter(x0(1),x0(2),50,'bo','filled');
        scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    5, 'bx');
        scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    5, 'ks');
        scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    5, 'md');
        scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    5, 'r*');
        
  
%%

load('Figure3b-0.6.mat')


 %% Plotting trajectories of each method: 

        plot_markersize = 10;
        plot_fontSize = 11;
        subplot(2,1,2)
        hold on
        plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        scatter(x0(1),x0(2),5*plot_markersize,'bo','filled');
        scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'md');
        h5 = scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'bx');
        h6 = scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    10*plot_markersize, 'r*');
        h7 = scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    10*plot_markersize, 'ks');


        
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        
        xlabel('Position ($z_1$)')
        ylabel('Position ($z_2$)')
%         title('\textbf{Trajectory}')
        set(gca,'FontSize',plot_fontSize)
%         legend([h1 h2 h3 h4 h5 h6 h7],{'Safe Set',...
%             'Target Set','Initial Position',...
%             'Piecewise affine approach - QP',...
%             'Ono2008 IRA Method ',...
%             'Piecewise affine approach - MIQP',...
%             sprintf('Blackmore11 PC Method, %i Particles',...
%             N)},...
%             'Location','SouthOutside','FontSize',plot_fontSize);
        
        box on;
        axis tight
        axis([-1.15 0.1 -1.15 0.05])
        set(gca,'XTick',sort([-1.15,get(gca,'XTick')]));
        set(gca,'YTick',sort([-1.15,get(gca,'YTick')]));
        set(gca,'FontSize',plot_fontSize)
        set(fig3,'Units','centimeters');
        set(fig3,'Position',[0 0 15 10]);
        
        
        handaxes3 = axes('position', [0.175 0.295 0.125 0.125]);
        hold on
        grid on
        plot(safe_set.slice([3,4], slice_at_vx_vy), 'color', 'y','alpha',0.5);
        plot(target_set.slice([3,4], slice_at_vx_vy), 'color', 'g','alpha',0.5);
        scatter(x0(1),x0(2),50,'bo','filled');
        scatter(ono_opt_mean_X(1:4:end),...
                    ono_opt_mean_X(2:4:end),...
                    5, 'bx');
        scatter(blackmore_opt_mean_X(1:4:end,1),...
                    blackmore_opt_mean_X(2:4:end,1),...
                    5, 'ks');
        scatter(onopwl_opt_mean_X(1:4:end),...
                    onopwl_opt_mean_X(2:4:end),...
                    5, 'md');
        scatter(pwa_opt_mean_X(1:4:end),...
                    pwa_opt_mean_X(2:4:end),...
                    5, 'r*');
        fig3 = tightfig(fig3);
        
        hgexport(fig3,'Figure4',hgexport('factorystyle'),'Format', 'png')
        hgexport(fig3,'Figure4',hgexport('factorystyle'),'Format', 'eps')
        saveas(gcf,'Figures/Figure4.fig','fig');
    
        
    