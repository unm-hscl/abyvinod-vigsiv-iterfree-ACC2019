function [cdf_approx_m, cdf_approx_c,lb_x] = RolleLerpClosedForm(Delta,lbdelta,maxlierror)

    myeps = 1e-10;
% Compute a piecewise-linear overapproximation of norminv(1-x) utilizing
% Rolle's theorem. 
    xe(1) = Delta;
    hval(1) = 0; % To be overwritten
    i = 1;
    syms x
    g = sqrt(2)* erfinv(2*(1 - x) -1 );
    G=diff(g,2);
    syms h
    syms y
    Hk=subs(G,y-h)*h^2/8 - maxlierror;
while xe(i)> lbdelta
    Hyk=subs(Hk,y,xe(i));
    fun = matlabFunction(Hyk);
    hval_intv = [myeps min(0.1,xe(i)-myeps)]; % initial interval
    options = optimset('Display','off'); % show iterations
    [hval(i)] = fzero(fun,hval_intv,options);
    xe(i+1) = xe(i) - hval(i);
    i = i+1;
end
    
    xe = fliplr(xe);
    for indx_x = 1:length(xe)
        if indx_x == length(xe)
            y_2 = norminv(1-xe(end));
            y_1 = norminv(1-xe(end-1));
            x_2 = xe(end);
            x_1 = xe(end-1);
            % y=mx+c where m and c are computed from secant end points
            cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
            cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
        else
            y_2 = norminv(1-xe(indx_x+1));
            y_1 = norminv(1-xe(indx_x));
            x_2 = xe(indx_x+1);
            x_1 = xe(indx_x);
            % y=mx+c where m and c are computed from secant end points
            cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
            cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
        end
    end 
    lb_x = min(xe);


end
