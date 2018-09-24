close all, clear, clc

xi = -1:0.25:1;
yi = -xi.^2;
dyi = -2*xi;

sdpvar x f
f = interp1(xi,yi,x,'graph',dyi);
plot(f >= -1,[x;f]);hold on
f = interp1(xi,yi,x,'graph');
plot(f >= -1,[x;f],'yellow');

xifine = -1:0.001:1;
yifine = -xifine.^2;
l = plot(xifine,yifine,'g','LineWidth',2)
plot(xi,yi,'k*','MarkerSize',10)
legend('$l^+_f(x)$','$l^+_f(x)$','f(x)','Location','EastOutside')

 plot_markersize = 10;
 plot_fontSize = 10;
 xlabel('x')
        ylabel('y')
        box on
        set(gca,'FontSize',plot_fontSize)
        
        set(gca,'Units','centimeters');
        set(gca,'Position',[0 0 4.4 2.2]);
        tightfig