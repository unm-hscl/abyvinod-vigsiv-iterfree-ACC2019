tic
clc, clear, close
clear;close all;
Delta = 0.2;
errorub = 1E-3;
errorlb = 9E-4;
h = 0.025;
[onopwl_invcdf_approx_m, onopwl_invcdf_approx_c,max_error_onopwl,...
        onopwl_lb_deltai] = RolleLerp(Delta,h,errorlb,errorub);
x=1.5E-6:1e-6:0.2;
y=max(onopwl_invcdf_approx_m'*x+onopwl_invcdf_approx_c');
y_true = norminv(1-x);
err=y-y_true;
figure(); plot(x,y,'ro-'); hold on; plot(x,y_true,'b*'); legend('PWL','True')
xlabel('x'),ylabel('norminv(1-x)');title('Approximation');
figure(); plot(x,err,'ro-');hold on;plot(x,1e-3*ones(length(x),1));
xlim([-0.2,0.5]);ylim([-1e-5,1.5e-2]);  xlabel('x'); ylabel('PWL-true curve');
title('Approximation quality');
if abs(max(err)) > abs(max_error_onopwl)
disp('Predicted error is smaller than the original error');
end
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'],...
min(err), max(err), length(onopwl_invcdf_approx_m));

toc