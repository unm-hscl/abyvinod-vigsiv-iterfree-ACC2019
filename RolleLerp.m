function [cdf_approx_m, cdf_approx_c,errorub,lb_x] = RolleLerp(Delta,h,errorlb,errorub)

% Compute a piecewise-linear overapproximation of norminv(1-x) utilizing
% Rolle's theorem. 

% Define the second derivative of error function as a standin for norminv.
% initialize x0 and x1 along with iteration, i, and the precision to move
% along: 
    x(1) = Delta ;
    i = 1;
    oLfl = 0;
    errorlb_reset = errorlb;
% while the resulting values do not produce obsure values, continue the PW
% approximation: 
    while erfinv(x(i)) ~= Inf || isnan(erfinv(x(i)))~= 1
% Reactivate the error checker to go into the the internal while loop:
        errorlb = errorlb_reset;
        active = 1;
        if length(x) > 1
            h = abs(x(i) - x(i-1));
        end
        
        x(i+1) = x(i)-h;

% If with the spacing we end up with x(i+1) being negative, then switch it
% to zero: 
        if abs(x(i)-x(i+1)) <= 0 
            break;
            x = x(1:end-1);
        end
% While the current error produced in the interval is not greater than
% error specified, continue!
        xt = linspace(x(i),x(i+1),100); 
        secderiv = zeros(1,length(xt));
        while active == 1 

            for j = 1:length(xt)
%                 secderiv(j) = 1/8*(x(i+1)-x(i))^2*sqrt(2)*...
%                     (2*pi*exp(2*erfinv(2*(1-xt(j)) - 1)^2)*...
%                     erfinv(2*(1-xt(j)) - 1));
                secderiv(j) = 1/8*(x(i+1)-x(i))^2*...
                    -2*2^(1/2)*pi*exp(2*erfinv(2*xt(j) - 1).^2)*...
                    erfinv(2*xt(j) - 1);
            end
            [max_sec_err ind1] = max(secderiv);
            [min_sec_err ind2] = min(secderiv);
            prev_max_sec = max_sec_err;
            prev_min_sec = min_sec_err;
            if  abs(min_sec_err - errorub) <= 2E-6
                active = 0;
                ind3 = find(secderiv==errorlb);
                max_sec_err_ar(i) = max_sec_err;
                x(i+1) = xt(ind3); 

            elseif sum(secderiv) ~= 0  &&...
                    min_sec_err >= errorlb && min_sec_err < errorub
                active = 1;
                errorlb = secderiv(ind2+1);
                xt = xt(ind2+1:end);
                secderiv = secderiv(ind2+1:end);
                secderiv = secderiv(1:ind1-1);
            end
            
        end
        
        

    
     fprintf('With point %i  x(i+1) = %1.4f\n',i,x(i+1))
    if length(x) > 1
        if isnan(x) == 1
            break;
        elseif abs(x(i)-x(i+1)) <= 0 || x(i+1) <= 1E-6 
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
