clc, clear, close all


%% norminv(1-delta)
maxlierror=1e-5;
g = @(z) -z.^2;
fun_monotone = 'mono-dec';
lower_bound = 1E-3;
upper_bound = 2; 
no_of_test_points = 1000;
function_handle = @(z) -g(z);

x_iter=linspace(lower_bound,upper_bound,no_of_test_points);
y_iter_true = -function_handle(x_iter);     % Only th

[PWA_overapprox_m, PWA_overapprox_c,...
    PWA_underapprox_m, PWA_underapprox_c] = getPWAOverAndUnderApprox(...
    lower_bound,upper_bound,maxlierror,function_handle,fun_monotone);

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