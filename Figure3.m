close all, clear,clc
load('Figure1a-0.2-Run3.mat')
        
    %% Plotting trajectories of each method: 

        plot_markersize = 10;
        plot_fontSize = 11;
        fig1 = figure(1);
        subplot(2,1,1)
        hold on
        polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
            [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
        P = Polyhedron('V',polyvertex);
        h1 = plot(P,'alpha',0.1);
        h11 = scatter(1,mean_x(1),plot_markersize*10,'filled');
        h2 = plot(2:(T+1),xtarget(1:2:end),'go','MarkerSize',...
            plot_markersize,'LineWidth',2);
        h3 = plot(2:(T+1),onopwl_opt_mean_X(1:2:end),'md',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h4 = plot(2:(T+1),ono_opt_mean_X(1:2:end),'b^',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h5 = plot(2:(T+1),pwa_opt_mean_X(1:2:end),'r*',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h6 = plot(2:(T+1),blackmore_opt_mean_X(1:2:end,1),...
            'ks','MarkerSize',plot_markersize);
        xlabel('Time Step')
        ylabel('Position (x)')
        
        set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
        set(groot, 'defaultLegendInterpreter','latex');
        set(groot, 'defaulttextInterpreter','latex');
        

%         title('\textbf{Trajectory}')
        [hleg, hobj, hout, mout] = legend([h1 h11 h2 h3 h4 h5 h6],{'Target tube',...
            'Initial position','Target trajectory',...
            ['Piecewise affine' newline 'approach - QP'],...
            ['Iterative risk' newline 'allocation (IRA)'],...
            ['Piecewise affine' newline 'approach - MIQP'],...
            ['Particle control' newline sprintf('(PC), %i particles',N)]},...
            'Location','EastOutside','FontSize',11);
                hobj(9).Children.MarkerSize = 10;

        box on;
        set(gca,'FontSize',plot_fontSize)
%         set(fig1,'Units','centimeters');
  
%%

load('Figure1b-0.4-Run3.mat')


        plot_markersize = 10;
        plot_fontSize = 11;
        subplot(2,1,2)
        hold on
        polyvertex =[1:T+1,1:T+1;[-gbig(1),-gbig(1:2:end)'],...
            [gbig(1),gbig(1:2:end)']]'; % Note: This needs MPT to run!!
           Polyhedron('V',polyvertex);
              plot(P,'alpha',0.1);
            scatter(1,mean_x(1),plot_markersize*10,'filled');
            plot(2:(T+1),xtarget(1:2:end),'go','MarkerSize',...
            plot_markersize,'LineWidth',2);
        h5 = plot(2:(T+1),pwa_opt_mean_X(1:2:end),'r*',...
            'LineWidth',1,'MarkerSize',plot_markersize);
        h6 = plot(2:(T+1),blackmore_opt_mean_X(1:2:end,1),...
            'ks','MarkerSize',plot_markersize);
        
%         set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%         set(groot, 'defaultLegendInterpreter','latex');
%         set(groot, 'defaulttextInterpreter','latex');
        
        xlabel('Time Step')
        ylabel('Position (x)')
        box on
        set(gca,'FontSize',plot_fontSize)
        
        set(fig1,'Units','centimeters');
        set(fig1,'Position',[0 0 15 10]);
        fig1 = tightfig(fig1);
        pause
        tightfig
        print('Figure3.pdf','-dpdf')
        hgexport(fig1,'Figure3',hgexport('factorystyle'),'Format', 'png')
        hgexport(fig1,'Figure3',hgexport('factorystyle'),'Format', 'eps')
        saveas(gcf,'Figures/Figure1b.fig','fig');
        
    