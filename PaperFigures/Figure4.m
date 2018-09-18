    clear
    close all
    clc
    
    plot_markersize = 9;
    plot_fontSize = 9;
    fig2 = figure(2);
    clf
    hold on
    h1 = plot(T_array,PWAWD,'md',...
        'LineWidth',1,'MarkerSize',plot_markersize);
    h2 = plot(T_array,OnoTTS,'bx','MarkerSize',...
        plot_markersize,'LineWidth',2);
    h3 = plot(T_array,PWAWOD,'r*',...
        'LineWidth',1,'MarkerSize',plot_markersize);
    h4 = plot(T_array,BlackmoreTTS,'ks',...
        'LineWidth',1,'MarkerSize',plot_markersize);

    xlabel('\textbf{Time Horizon}')
    ylabel('\textbf{Time to Solve (seconds)}')
    legend([h1 h2 h3 h4],{'Piecewise affine approach - QP',...
        'Ono2008 IRA Method',...
        'Piecewise linear approach - MIQP',...
        sprintf('Blackmore11 PC Method, %i Particles',N)},...
        'Location','SouthOutside','FontSize',plot_fontSize);
    box on;
    set(gca,'FontSize',plot_fontSize)
    set(gca, 'YScale', 'log')