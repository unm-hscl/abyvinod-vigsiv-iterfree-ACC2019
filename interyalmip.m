close all, clear, clc

g_matlabfun = @(x) -x.^2;
xi = -1:0.25:1;
xii = -1:0.001:1;
yi = -xi.^2;
dyi = -2*xi;
y = g_matlabfun(xii);

sdpvar x f
f = interp1(xi,yi,x,'graph',dyi);
plot(f >= -1,[x;f],'yellow');hold on
actual_poly = Polyhedron('V',[xii', y';-1,-1]);
plot(actual_poly,'color','m'); 
f = interp1(xi,yi,x,'graph');
plot(f >= -1,[x;f],'cyan');

xifine = -1:0.001:1;
yifine = -xifine.^2;
% l = plot(xifine,yifine,'g','LineWidth',2)

plot(xi,yi,'k*','MarkerSize',10)
legend('$l^+_f(x)$','f(x)','$l^-_f(x)$','Location','EastOutside')

 plot_markersize = 10;
 plot_fontSize = 10;
 axis([0 0.4 -0.15 0])
 xlabel('x')
        ylabel('y')
        box on
        set(gca,'FontSize',plot_fontSize)
        
        set(gca,'Units','centimeters');
        set(gca,'Position',[0 0 4.4 2.2]);
        tightfig