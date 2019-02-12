%% abyvnod-vigsiv-iterfree-ACC2019: Figure 1
    % This plot just shows an example of -x^2 with our algorithm
    % REQUIRED DEPENDENCIES: - MATLAB symbolic toolbox
    %                        - MATLAB global optimization toolbox


clc, clear, close all


%% norminv(1-delta)
maxlierror=1e-2;
g = @(z) -z.^2;
fun_monotone = 'mono-dec';
lower_bound = maxlierror;
upper_bound = 4; 
no_of_test_points = 100;
function_handle = @(z) g(z);

x_iter=linspace(lower_bound,upper_bound,no_of_test_points);
y_iter_true = function_handle(x_iter);     % Only th

[PWA_overapprox_m, PWA_overapprox_c,...
    PWA_underapprox_m, PWA_underapprox_c] = getPWAOverAndUnderApprox(...
    lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);

% Added a negative sign since function_handle is the negation of g(z)
% Also, because of the image across x-axis due to negation, the definitions
% of approximations get flipped
if length(PWA_overapprox_m) > 1
    y_iter_underapprox  = -max(-PWA_overapprox_m' *x_iter - PWA_overapprox_c');
    y_iter_overapprox = -max(-PWA_underapprox_m'*x_iter - PWA_underapprox_c');
else
    y_iter_underapprox  = PWA_overapprox_m' *x_iter + PWA_overapprox_c';
    y_iter_overapprox = PWA_underapprox_m'*x_iter + PWA_underapprox_c';
end    
err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%% Plotting

fig1 = figure(1);
hold on
polyvertex1 =[x_iter,x_iter(1);y_iter_overapprox,y_iter_overapprox(end)]'; % Note: This needs MPT to run!!
P1 = Polyhedron('V',polyvertex1);
polyvertex2 =[x_iter,x_iter(1);y_iter_true,y_iter_true(1)]'; % Note: This needs MPT to run!!
P2 = Polyhedron('V',polyvertex2);
polyvertex3 =[x_iter,x_iter(1);y_iter_underapprox,y_iter_underapprox(end)]'; % Note: This needs MPT to run!!
P3 = Polyhedron('V',polyvertex3);
h1 = plot(P3,'color','y');
hold on
h2 = plot(P2,'color','b');
h3 = plot(P1,'color','r');
axis([0.009 0.275 -0.11 0.025])
% yticks([-1 -0.5 0])
set(fig1,'Units','centimeters');
set(fig1,'Position',[0 0 12 4]);
grid on
box on
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaulttextInterpreter','latex');
legend({'$l_f^+\geq s$','$f(x)\geq s$','$l_f^-\geq s$'},'Location','eastoutside','FontSize',12)
xlabel('x')
ylabel('s')
set(gca,'FontSize',12)
fig1 = tightfig(fig1);



