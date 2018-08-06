tic
clc, clear, close
clear;close all;
Delta = 0.2;
maxlierror = 1E-3;
lbdelta = 1E-3;
[onopwl_invcdf_approx_m, onopwl_invcdf_approx_c,...
        onopwl_lb_deltai] = RolleLerpClosedForm(Delta,lbdelta,maxlierror);
x=onopwl_lb_deltai:1e-6:0.2;
y=max(onopwl_invcdf_approx_m'*x+onopwl_invcdf_approx_c');
y_true = norminv(1-x);
err=y-y_true;
%%
figure(); plot(x,y,'ro-'); hold on; plot(x,y_true,'b*'); legend('PWL','True')
xlabel('x'),ylabel('norminv(1-x)');title('Approximation');
figure(); plot(x,err,'ro-');hold on;plot(x,maxlierror*ones(length(x),1));
xlim([-1E-5,Delta]);ylim([-1e-5,1.5*maxlierror]);  xlabel('x'); ylabel('PWL-true curve');
title('Approximation quality');
fprintf(['Min error: %1.4e | Max error: %1.4e | No. of ineq: %d\n'],...
min(err), max(err), length(onopwl_invcdf_approx_m));

toc