function [cdf_overapprox_m, cdf_overapprox_c, cdf_underapprox_m, cdf_underapprox_c] = getPWLOverAndUnderApprox(lbdelta,ubdelta,maxlierror)

    myeps = 1e-10;
    %% Initialization
    knots_underapprox(1) = lbdelta;
    knots_overapprox(1) = lbdelta; % To be overwritten
    hval(1) = 0;                   % To be overwritten
    j = 1;
    search_interval = [myeps 0.1];
    fzero_options = optimset('Display','off'); % show iterations
                
    %% Symbolic function definitions
    syms x;
    syms h;
    % TODO: Get g as an argument
    % Function definition for Phi^{-1}(1-x)
    g = -sqrt(2)* erfinv(2*(1 - x) -1 );
    % First derivative
    g1diff=diff(g,1);
    % Second derivative: Function in (x,h) that needs to be solved for h given x,
    % d^2/dh^2 g(x+h) == 8 * eta^-;
    g2diff=diff(g,2);
    g_maxerror=subs(g2diff,x+h)*h^2 + 8*maxlierror;
    
    % Iterate till we reach the end
    while knots_underapprox(j)< ubdelta
        % Set up the function to be solved by fzero
        g_maxerror_at_j = subs(g_maxerror,x,knots_underapprox(j));
        g_maxerror_at_j_mf = matlabFunction(g_maxerror_at_j);
        % Set up the search interval to not exceed ubdelta
        search_interval(2) = min(search_interval(2), ubdelta - knots_underapprox(j));
        % Solve for h
%         fprintf('%1.4e\n',search_interval);
%         fprintf('%1.4e\n',g_maxerror_at_j_mf(search_interval));
        [hval(j)] = fzero(g_maxerror_at_j_mf,search_interval,fzero_options);
        % Ensure that the knot chosen does not exceed
        knots_underapprox(j+1) = min(knots_underapprox(j) + hval(j),ubdelta);
        
        %% Construction of the underapproximation of a concave function --- Lagrange interpolation
        y_2 = norminv(1-knots_underapprox(j+1));
        y_1 = norminv(1-knots_underapprox(j));
        x_2 = knots_underapprox(j+1);
        x_1 = knots_underapprox(j);
        % Lagrange linear interpolation
        cdf_underapprox_m(j) = (y_2 - y_1)/(x_2 - x_1);
        cdf_underapprox_c(j) = y_1 - cdf_underapprox_m(j) * x_1;
        
%         %% Construction of the overapproximation of a concave function --- first-order Taylor series
%         % Set up the gradient function
%         g1diff_at_j = subs(g1diff - cdf_underapprox_m(j),x,knots_underapprox(j));
%         g1diff_at_j_mf = matlabFunction(g1diff_at_j);
%         % Solve for Taylor series operating point
%         knots_overapprox(j) = fzero(g1diff_at_j_mf,[x_1 x_2],fzero_options);        
%         cdf_overapprox_m(j) = double(subs(g1diff,knots_overapprox(j)));
%         cdf_overapprox_c(j) = double(subs(g,knots_overapprox(j))) - double(subs(g1diff,knots_overapprox(j))) * knots_overapprox(j);
        
        %% Increment j
        j = j+1;       
    end
end
