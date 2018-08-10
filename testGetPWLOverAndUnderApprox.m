clc, clear, close
clear;close all;

Delta = 0.2;

%% norminv(1-delta)
maxlierror=1e-2;
K=6;
g = @(z) sqrt(2)* erfinv(2*(1 - z) -1 );
fun_monotone = 'mono-inc';
lower_bound = 1e-5;
upper_bound = Delta; 
function_handle = @(z) -g(z);

x_iter=lower_bound:1e-7:upper_bound;
y_iter_true = -function_handle(x_iter);     % Only th

tic
[cdf_overapprox_m, cdf_overapprox_c, cdf_underapprox_m, cdf_underapprox_c] = getPWLOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc
% Added a negative sign since function_handle is the negation of g(z)
y_iter_underapprox  = -min(cdf_overapprox_m' *x_iter + cdf_overapprox_c');
y_iter_overapprox = -min(cdf_underapprox_m'*x_iter + cdf_underapprox_c');

err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(1); 
clf
subplot(2,1,1);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-'); 
plot(x_iter,y_iter_underapprox,'gd-'); 
leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
set(leg,'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=-\Phi^{-1}(1-x)$','interpreter','latex');
subplot(2,1,2);
plot(x_iter,maxlierror*ones(length(x_iter),1));
hold on;
plot(x_iter,err_overapprox,'ro-');
plot(x_iter,err_underapprox,'gd-');
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
leg=legend('Max allowable error','$\ell_f^+(x) - f(x)$','$f(x) - \ell_f^-(x)$');
set(leg,'interpreter','latex');
title('Approximation quality','interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(cdf_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(cdf_underapprox_m));

%% log(Phi(x))
maxlierror=1e-2;
K=6;
logphi = @(z) log(normcdf(z));
fun_monotone = 'mono-inc';
lower_bound = -K;
upper_bound = K; 
function_handle = logphi;

x_iter=lower_bound:1e-4:upper_bound;
y_iter_true = function_handle(x_iter);

tic
[cdf_overapprox_m, cdf_overapprox_c, cdf_underapprox_m, cdf_underapprox_c] = getPWLOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc

y_iter_overapprox  = min(cdf_overapprox_m' *x_iter + cdf_overapprox_c');
y_iter_underapprox = min(cdf_underapprox_m'*x_iter + cdf_underapprox_c');

err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(2); 
clf
subplot(2,1,1);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-'); 
plot(x_iter,y_iter_underapprox,'gd-'); 
leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
set(leg,'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=\log\Phi(x)$','interpreter','latex');
subplot(2,1,2);
plot(x_iter,maxlierror*ones(length(x_iter),1));
hold on;
plot(x_iter,err_overapprox,'ro-');
plot(x_iter,err_underapprox,'gd-');
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
leg=legend('Max allowable error','$\ell_f^+(x) - f(x)$','$f(x) - \ell_f^-(x)$');
set(leg,'interpreter','latex');
title('Approximation quality','interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(cdf_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(cdf_underapprox_m));

%% log(1-x)
maxlierror=1e-3;
logOneMinusZ = @(z) log(1-z);
fun_monotone = 'mono-dec';
lower_bound = log(1-Delta);
upper_bound = log(normcdf(K)); 
function_handle = logOneMinusZ;

x_iter=lower_bound:1e-6:upper_bound;
y_iter_true = function_handle(x_iter);

tic
[cdf_overapprox_m, cdf_overapprox_c, cdf_underapprox_m, cdf_underapprox_c] = getPWLOverAndUnderApprox(lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);
toc

y_iter_overapprox  = min(cdf_overapprox_m' *x_iter + cdf_overapprox_c');
y_iter_underapprox = min(cdf_underapprox_m'*x_iter + cdf_underapprox_c');

err_overapprox=y_iter_overapprox-y_iter_true;
err_underapprox=y_iter_true-y_iter_underapprox;

%%
figure(3); 
clf
subplot(2,1,1);
plot(x_iter,y_iter_true,'bs-'); 
hold on; 
plot(x_iter,y_iter_overapprox,'ro-'); 
plot(x_iter,y_iter_underapprox,'gd-'); 
leg=legend('True','$\ell_f^+(x)$','$\ell_f^-(x)$');
set(leg,'interpreter','latex');
xlabel('$x$','interpreter','latex');
ylabel('$f(x)$','interpreter','latex');
title('Approximation of $f(x)=\log(1-x)$','interpreter','latex');
subplot(2,1,2);
plot(x_iter,maxlierror*ones(length(x_iter),1));
hold on;
plot(x_iter,err_overapprox,'ro-');
plot(x_iter,err_underapprox,'gd-');
ylim([-0.5*maxlierror,1.5*maxlierror]);  
xlabel('$x$','interpreter','latex');
ylabel('Error');
leg=legend('Max allowable error','$\ell_f^+(x) - f(x)$','$f(x) - \ell_f^-(x)$');
set(leg,'interpreter','latex');
title('Approximation quality','interpreter','latex');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_overapprox), max(err_overapprox), length(cdf_overapprox_m));
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'], min(err_underapprox), max(err_underapprox), length(cdf_underapprox_m));
