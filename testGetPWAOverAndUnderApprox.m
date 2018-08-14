clc;clear;close all;

underlineDelta = 1e-5;
Delta = 0.5;
K=5;
DeltaMax = 0.7;

no_of_test_points = 1000;
legend_str = {'Bounds on approx. error','$\ell_f^+(x) - f(x)$','$f(x) - \ell_f^-(x)$'};
title_str2 = '';%'PWA approximation quality of';
errorboundMarkerSize = 10;
errorMarkerSize = 5;
approxMarkerSize = 5;

%% norminv(1-delta)
maxlierror=1e-3;
g = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
fun_monotone = 'mono-inc';
lower_bound = underlineDelta;
upper_bound = Delta; 
function_handle = @(z) -g(z);

x_iter=linspace(lower_bound,upper_bound,no_of_test_points);
y_iter_true = -function_handle(x_iter);     % Only th

tic
[PWA_overapprox_m, PWA_overapprox_c, PWA_underapprox_m, PWA_underapprox_c] = getPWAOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc
% Added a negative sign since function_handle is the negation of g(z)
% Also, because of the image across x-axis due to negation, the definitions
% of approximations get flipped
if length(PWA_overapprox_m) > 1
    y_iter_underapprox  = max(-PWA_overapprox_m' *x_iter - PWA_overapprox_c');
    y_iter_overapprox = max(-PWA_underapprox_m'*x_iter - PWA_underapprox_c');
else
    y_iter_underapprox  = -(PWA_overapprox_m' *x_iter + PWA_overapprox_c');
    y_iter_overapprox = -(PWA_underapprox_m'*x_iter + PWA_underapprox_c');
end    
err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(1); 
clf
subplot(3,1,1);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-','MarkerSize',approxMarkerSize);
plot(x_iter,y_iter_underapprox,'gd-','MarkerSize',approxMarkerSize);
% leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
% set(leg,'interpreter','latex');
xlim([x_iter(1),x_iter(end)]);  
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=\Phi^{-1}(1-x)$','interpreter','latex');
figure(2)
clf
subplot(3,1,1);
h1=plot(x_iter,maxlierror*ones(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
hold on;
h2=plot(x_iter,zeros(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
h3=plot(x_iter,err_overapprox,'ro-','MarkerSize',errorMarkerSize);
h4=plot(x_iter,err_underapprox,'gd-','MarkerSize',errorMarkerSize);
xlim([x_iter(1),x_iter(end)]);  
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
% leg=legend([h1 h3 h4],legend_str);
% set(leg,'interpreter','latex');
% title(strcat(title_str2,' of $f(x)=-\Phi^{-1}(1-x)$'),'interpreter','latex');
title(sprintf('%s $f(x)=\\Phi^{-1}(1-x)$ $x\\in[\\underline{\\delta},\\Delta]$, $\\underline{\\delta}=%1.1e$, $\\Delta=%1.2f$, $\\eta=%1.2f$  with %d inequalities',title_str2, underlineDelta, Delta, maxlierror, length(PWA_overapprox_m)),'interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(PWA_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(PWA_underapprox_m));

%% log(Phi(x))
maxlierror=1e-2;
logphi = @(z) log(normcdf(z));
fun_monotone = 'mono-inc';
lower_bound = -K;
upper_bound = norminv(exp(-maxlierror)); 
function_handle = logphi;

x_iter=linspace(lower_bound,upper_bound,no_of_test_points);
y_iter_true = function_handle(x_iter);

tic
[PWA_overapprox_m, PWA_overapprox_c, PWA_underapprox_m, PWA_underapprox_c] = getPWAOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc

if length(PWA_overapprox_m) > 1
    y_iter_overapprox  = min(PWA_overapprox_m' *x_iter + PWA_overapprox_c');
    y_iter_underapprox = min(PWA_underapprox_m'*x_iter + PWA_underapprox_c');
else
    y_iter_overapprox  = PWA_overapprox_m' *x_iter + PWA_overapprox_c';
    y_iter_underapprox = PWA_underapprox_m'*x_iter + PWA_underapprox_c';
end    

err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(1)
subplot(3,1,2);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-','MarkerSize',approxMarkerSize);
plot(x_iter,y_iter_underapprox,'gd-','MarkerSize',approxMarkerSize);
leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
set(leg,'interpreter','latex','Location','EastOutside');
xlim([x_iter(1),x_iter(end)]);  
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=\log\Phi(x)$','interpreter','latex');
figure(2)
subplot(3,1,2);
h1=plot(x_iter,maxlierror*ones(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
hold on;
h2=plot(x_iter,zeros(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
h3=plot(x_iter,err_overapprox,'ro-','MarkerSize',errorMarkerSize);
h4=plot(x_iter,err_underapprox,'gd-','MarkerSize',errorMarkerSize);
xlim([x_iter(1),x_iter(end)]);  
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
leg=legend([h1 h3 h4],legend_str);
set(leg,'interpreter','latex','Location','EastOutside');
title(sprintf('%s $f(x)=\\log\\Phi(x)$ $x\\in[-K,\\Phi^{-1}(e^{-\\eta})]$, $K=%1.2f$, $\\eta=%1.2f$ with %d inequalities',title_str2, K, maxlierror, length(PWA_overapprox_m)),'interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(PWA_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(PWA_underapprox_m));

%% log(1-x)
maxlierror=1e-2;
logOneMinusZ = @(z) log(1-z);
fun_monotone = 'mono-dec';
lower_bound = log(1-DeltaMax);
upper_bound = -maxlierror; 
function_handle = logOneMinusZ;

x_iter=linspace(lower_bound,upper_bound,no_of_test_points);
y_iter_true = function_handle(x_iter);

tic
[PWA_overapprox_m, PWA_overapprox_c, PWA_underapprox_m, PWA_underapprox_c] = getPWAOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc
if length(PWA_overapprox_m) > 1
    y_iter_overapprox  = min(PWA_overapprox_m' *x_iter + PWA_overapprox_c');
    y_iter_underapprox = min(PWA_underapprox_m'*x_iter + PWA_underapprox_c');
else
    y_iter_overapprox  = PWA_overapprox_m' *x_iter + PWA_overapprox_c';
    y_iter_underapprox = PWA_underapprox_m'*x_iter + PWA_underapprox_c';
end    
err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(1); 
subplot(3,1,3);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-','MarkerSize',approxMarkerSize);
plot(x_iter,y_iter_underapprox,'gd-','MarkerSize',approxMarkerSize);
% leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
% set(leg,'interpreter','latex');
xlim([x_iter(1),x_iter(end)]);  
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=\log(1-x)$','interpreter','latex');
figure(2)
subplot(3,1,3);
h1=plot(x_iter,maxlierror*ones(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
hold on;
h2=plot(x_iter,zeros(length(x_iter),1),'b--','MarkerSize',errorboundMarkerSize);
h3=plot(x_iter,err_overapprox,'ro-','MarkerSize',errorMarkerSize);
h4=plot(x_iter,err_underapprox,'gd-','MarkerSize',errorMarkerSize);
xlim([x_iter(1),x_iter(end)]);  
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
% leg=legend([h1 h3 h4],legend_str);
% set(leg,'interpreter','latex');
title(sprintf('%s $f(x)=\\log(1-x)$ $x\\in[\\log(1-\\Delta),0]$, $\\Delta=%1.2f$, $\\eta=%1.2f$  with %d inequalities',title_str2,DeltaMax, maxlierror, length(PWA_overapprox_m)),'interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(PWA_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(PWA_underapprox_m));

%% Save figures
figure(2);
saveas(gcf,'Figures/ErrorPlots.png','png');
saveas(gcf,'Figures/ErrorPlots.fig','fig');
