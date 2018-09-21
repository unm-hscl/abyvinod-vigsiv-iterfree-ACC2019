    clear
    close all
    clc
    load('../Figure1a-0.2-Run3.mat','T','gbig','polyvertex','P','mean_x','xtarget','ono_opt_mean_X','onopwl_opt_mean_X','pwa_opt_mean_X','blackmore_opt_mean_X','N');
    plot_markersize = 10;
    plot_fontSize = 10;
    fig1 = figure(1);
    clf
    subplot(2,2,1);
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

    xlabel('Time step')
    ylabel('Position (x)')
%         title('\textbf{Trajectory}')
    leg=legend([h1 h11 h2 h3 h4 h5 h6],{'Target Tube',...
        'Initial state','Target Trajectory',...
        ['Piecewise affine' newline 'approach - QP'],...
        ['Iterative risk' newline 'allocation (IRA)'],...
        ['Piecewise affine' newline 'approach - MIQP'],...
        ['Particle control' newline sprintf('(PC) %i Particles',N)]},'FontSize',plot_fontSize);
    box on;
    set(gca,'FontSize',plot_fontSize);
    set(leg,'Position', [0.48 0.0 0.25 0.9]);
    load('../Figure1b-0.6-Run3.mat','T','gbig','polyvertex','P','mean_x','xtarget','ono_opt_mean_X','onopwl_opt_mean_X','pwa_opt_mean_X','blackmore_opt_mean_X');
    subplot(2,2,3);
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
    xlabel('Time step')
    ylabel('Position (x)')
    box on;
    set(gca,'FontSize',plot_fontSize);