function [PWA_overapprox_m,...
          PWA_overapprox_c,...
          PWA_underapprox_m,...
          PWA_underapprox_c] = getPWAOverAndUnderApprox(lb,ub,eta,g_matlabfun,funtype)
    warning('getPWLOverAndUnderApprox does not check for monotonicity and concavity claims.');

    %% Initialization
    knots_underapprox(1) = lb;
    j = 1;
    search_interval_min = 0;
    search_interval_max = ub - lb;
    fzero_options = optimset('Display','off'); % show iterations
                
    %% Symbolic function definitions
    % TODO: Get g as an argument
    % Function definition for Phi^{-1}(1-x)
    % First derivative
    syms x;
    syms h;
    syms hmax;
    g = g_matlabfun(x);
    g1diff=diff(g,1);
    % Second derivative: Function in (x,h) that needs to be solved for h given x,
    % d^2/dh^2 g(x+h) == 2 * eta^-;
    g2diff=diff(g,2);
    % Due to concavity of g, the error becomes more negative for larger h
    
    % Iterate till we reach the end
    while knots_underapprox(j)< ub

        if strcmpi(funtype,'mono-inc')
            hval = sqrt(-8*eta/double(subs(g2diff,x,knots_underapprox(j))));
            knots_underapprox(j+1) = knots_underapprox(j) + hval;
        elseif strcmpi(funtype,'mono-dec')
            g_maxerror=subs(g2diff,x+h)*h^2 + 8*eta;
            g_maxerror_at_j = subs(g_maxerror,x,knots_underapprox(j));
            g_maxerror_at_j_mf = @(z) double(subs(g_maxerror_at_j,z));
            % Set up the search interval to not exceed ubdelta
            search_interval = [search_interval_min  min(search_interval_max,ub-knots_underapprox(j))];
            % Solve for h
            try
                hval = fzero(g_maxerror_at_j_mf,search_interval,fzero_options);
                knots_underapprox(j+1) = min(knots_underapprox(j) + hval,ub);
            catch
                % If errored, then the search interval doesn't see a sign flip
                % => we have reached the end
%                 hval = ub - knots_underapprox(j);
                knots_underapprox(j+1) = ub;
            end
%             fprintf('%2d. Searched in [%1.4e,%1.4e] to get %1.4e\n',j, knots_underapprox(j)+search_interval,knots_underapprox(j+1));
        else
            disp('Untested');
%             %% Use these steps in a general case --- Untested
%             g_maxerror=subs(g2diff,x+hmax)*h^2 + 8*eta;
%             g_maxerror_at_j = subs(g_maxerror,x,knots_underapprox(j));
%             g_maxerror_at_j_mf = matlabFunction(g_maxerror_at_j);
%             % Set up the search interval to not exceed ubdelta
%             search_interval = [search_interval_min  search_interval_max];
%             % Solve for h
%             try
%                 max_soln = fmincon(@(v) g_maxerror_at_j_mf(v(1),v(2)),[1,1],[],[],[],[],search_interval_min*ones(2,1),search_interval_max*ones(2,1));
%                 hval(j) = max_soln(2);
%                 knots_underapprox(j+1) = min(knots_underapprox(j) + hval(j),ubdelta);
%             catch
%                 % If errored, then the search interval doesn't see a sign flip
%                 % => we have reached the end
%                 hval(j) = ubdelta - knots_underapprox(j);
%                 knots_underapprox(j+1) = ubdelta;
%             end
%             fprintf('%2d. Searched in [%1.4e,%1.4e] to get %1.4e\n',j, knots_underapprox(j)+search_interval,knots_underapprox(j+1));
        end
        
        %% Construction of the underapproximation of a concave function --- Lagrange interpolation
        y_2 = double(subs(g,knots_underapprox(j+1)));
        y_1 = double(subs(g,knots_underapprox(j)));
        x_2 = knots_underapprox(j+1);
        x_1 = knots_underapprox(j);
        % Lagrange linear interpolation
        PWA_underapprox_m(j) = (y_2 - y_1)/(x_2 - x_1);
        PWA_underapprox_c(j) = y_1 - PWA_underapprox_m(j) * x_1;
        
        %% Construction of the overapproximation of a concave function --- first-order Taylor series
        % Set up the gradient function
        g1diff_at_j = g1diff - PWA_underapprox_m(j);
        % NOTE: Using matlabFunction fails for large -K in \log\Phi (see
        % the end of the code)
        g1diff_at_j_mf = @(z) double(subs(g1diff_at_j,z));
        % Set up the search interval as [x(j), x(j+1)]
        search_interval = [x_1  x_2];
        % Search for hvalgrad such that f'(x(j)+hvalgrad) = c_{j} ---
        % existence guaranteed by mean value theorem
        [x_grad_match] = fzero(g1diff_at_j_mf,search_interval,fzero_options);
        %cdf_underapprox_m(j) is the same as double(subs(g1diff,x_grad_match));
        PWA_overapprox_m(j) = PWA_underapprox_m(j);
        PWA_overapprox_c(j) = double(subs(g,x_grad_match)) - PWA_overapprox_m(j) * x_grad_match;        
        
        %% Increment j
        j = j+1;       
    end
end
% syms z;
% logphi1diff = diff(log(normcdf(z)),1); 
% logphi1diff_mf = matlabFunction(logphi1diff); 
% figure();
% subplot(2,1,1);
% plot(-20:1e-2:10,logphi1diff_mf(-20:1e-2:10));
% subplot(2,1,2);
% plot(-20:1e-2:10,double(subs(logphi1diff,-20:1e-2:10)))
        