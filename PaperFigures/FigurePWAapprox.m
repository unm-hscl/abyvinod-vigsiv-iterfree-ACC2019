clear
clc

g_matlabfun = @(x) -x.^2;
lb = 0;
ub = 1;
desired_accuracy = 0.1;
hessian_monotone = 'mono-inc';

[PWA_overapprox_m,...
 PWA_overapprox_c,...
 PWA_underapprox_m,...
 PWA_underapprox_c,...
 knots_underapprox] = getPWAOverAndUnderApprox(lb,...
    ub,...
    desired_accuracy,...
    g_matlabfun,...
    hessian_monotone);

x = lb:1e-3:ub;
y = g_matlabfun(x);
over_y_all = PWA_overapprox_m'*x + PWA_overapprox_c';
over_y = min(over_y_all);
under_y_all = PWA_underapprox_m'*x + PWA_underapprox_c';
under_y = min(under_y_all);

%% Plot all of them
figure(1);
clf
plot(x,y,'Linewidth',2);
hold on
plot(x,over_y,'--','Linewidth',2);
plot(x,under_y,':','Linewidth',2);
scatter(knots_underapprox,g_matlabfun(knots_underapprox),30,'r','filled');
mid_points = (knots_underapprox(2:end)+knots_underapprox(1:end-1))/2;
scatter(mid_points,g_matlabfun(mid_points),30,'c','filled');
leg=legend('$f(x)$','$\ell_f^+(x)$','$\ell_f^-(x)$');
set(leg,'interpreter','latex');
set(gca,'FontSize',20)
axis equal
axis off

%% Second plot
fig = figure(2);
clf
hold on
over_approx = Polyhedron('H',[-PWA_overapprox_m               , 0,-1;
                              ones(1,length(PWA_overapprox_m)),-1, 0;
                              PWA_overapprox_c                  1  0]');
plot(over_approx,'color','y');                 
actual_poly = Polyhedron('V',[x', y';
                              0,-1]);
plot(actual_poly,'color','b');                 
under_approx = Polyhedron('V',[knots_underapprox             , 0;
                               g_matlabfun(knots_underapprox),-1]');
plot(under_approx,'color','r');
leg=legend('$\ell_f^+(x)\geq s$','$f(x)\geq s$','$\ell_f^-(x)\geq s$','Location','EastOutside');
set(leg,'interpreter','latex');
set(gca,'FontSize',20)
set(gca,'XTick',[0:0.2:1]);
xlabel('x','interpreter','latex');
ylabel('s','interpreter','latex');
axis([0 1 -1 0.25])
box on
grid on
% axis off

plot_markersize = 10;
plot_fontSize = 10;

box on
set(gca,'FontSize',plot_fontSize)

set(gca,'Units','centimeters');
set(gca,'Position',[0 0 4.4 2.2]);
tightfig

hgexport(fig,'fit',hgexport('factorystyle'),'Format', 'png')
hgexport(fig,'fit',hgexport('factorystyle'),'Format', 'eps')
saveas(fig,'Figures/fit.fig','fig');