function [cdf_approx_m, cdf_approx_c,errorub,lb_x] = RolleLerp()

% Compute a piecewise-linear overapproximation of norminv(1-x) utilizing
% Rolle's theorem. 

% Define the second derivative of error function as a standin for norminv.
% initialize x0 and x1 along with iteration, i, and the precision to move
% along: 
    x(1) = 0.2;
    h = 0.025;
    i = 1;
    errorub = 1E-2;
    errorlb = errorub/10;
    oLfl = 0;
% while the resulting values do not produce obsure values, continue the PW
% approximation: 
    while erfinv(x(i)) ~= Inf || isnan(erfinv(x(i)))~= 1
% Reactivate the error checker to go into the the internal while loop:
        active = 1;
        x(i+1) = x(i)-h;

% If with the spacing we end up with x(i+1) being negative, then switch it
% to zero: 
        if abs(x(i)-x(i+1)) <= 0 
            break;
            x = x(1:end-1);
        end
% While the current error produced in the interval is not greater than
% error specified, continue!
        while active == 1 
            xt = linspace(x(i),x(i+1),200); 
            for j = 1:length(xt)
                secderiv(j) = 1/8*(x(i+1)-x(i))^2*sqrt(2)*...
                    (2*pi*exp(2*erfinv(2*(1-xt(j)) - 1)^2)*...
                    erfinv(2*(1-xt(j)) - 1));
            end
            [max_sec_err ind2] = max(secderiv);
            if  errorlb <= max_sec_err && max_sec_err <= errorub
                active = 0;
                max_sec_err_ar(i) = max_sec_err;

            elseif sum(secderiv) ~= 0
                active = 1;
                x(i+1) = xt(ind2-1);
            else
                break;
            end
            
        end
        
        

    
%      fprintf('With point %i  x(i+1) = %1.4f\n',i,x(i+1))
    if length(x) > 1
        if isnan(x) == 1
            break;
        elseif abs(x(i)-x(i+1)) == 0
            x = x(1:end-1);
            break;
        end
    end
    i = i+1;

    end
    
    
    
    
    h = 0.0002;
    errorub = 1E-2;
    errorlb = errorub/10;
    oLfl = 0;
% while the resulting values do not produce obsure values, continue the PW
% approximation: 
    while erfinv(x(i)) ~= Inf || isnan(erfinv(x(i)))~= 1
% Reactivate the error checker to go into the the internal while loop:
        active = 1;
        x(i+1) = x(i)-h;

% If with the spacing we end up with x(i+1) being negative, then switch it
% to zero: 
        if abs(x(i)-x(i+1)) <= 0 
            break;
            x = x(1:end-1);
        end
% While the current error produced in the interval is not greater than
% error specified, continue!
        while active == 1 
            xt = linspace(x(i),x(i+1),100); 
            for j = 1:length(xt)
                secderiv(j) = 1/8*(x(i+1)-x(i))^2*sqrt(2)*...
                    (2*pi*exp(2*erfinv(2*(1-xt(j)) - 1)^2)*...
                    erfinv(2*(1-xt(j)) - 1));
            end
            [max_sec_err ind2] = max(secderiv);
            prev_max_sec = max_sec_err;
            if   max_sec_err <= errorub &&...
                    max_sec_err == prev_max_sec
                active = 0;
                max_sec_err_ar(i) = max_sec_err;

            elseif sum(secderiv) ~= 0
                active = 1;
                x(i+1) = xt(ind2-1);
            else
                break;
            end
            
        end
        
        

    
%      fprintf('With point %i  x(i+1) = %1.6f\n',i+1,x(i+1))
    if length(x) > 1
        if isnan(x) == 1
            break;
        elseif abs(x(i)-x(i+1)) == 0
            x = x(1:end-1);
            break;
        end
    end
    i = i+1;

    end
    
    x = fliplr(x);
    for indx_x = 1:length(x)
        if indx_x == length(x)
            y_2 = norminv(1-x(end));
            y_1 = norminv(1-x(end-1));
            x_2 = x(end);
            x_1 = x(end-1);
            % y=mx+c where m and c are computed from secant end points
            cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
            cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
        else
            y_2 = norminv(1-x(indx_x+1));
            y_1 = norminv(1-x(indx_x));
            x_2 = x(indx_x+1);
            x_1 = x(indx_x);
            % y=mx+c where m and c are computed from secant end points
            cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
            cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
        end
    end 
    lb_x = min(x);

end
