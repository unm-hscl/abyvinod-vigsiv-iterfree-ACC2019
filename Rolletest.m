clc, clear, close

% Define the second derivative of error function as a standin for norminv.
% initialize x0 and x1 along with iteration, i, and the precision to move
% along: 
    x(1) = 0.5;
    h = 0.1;
    i = 1;
    precision = 1E-4; 
    error = 1E-4;
% while the resulting values do not produce obsure values, continue the PW
% approximation: 
    while erfinv(x(i)) ~= Inf || isnan(erfinv(x(i)))~= 1
% Reactivate the error checker to go into the the internal while loop:
        active = 1;
% Let initial x1 be shifted by some spacing h: 
        x(i+1) = x(i)-h;
% If with the spacing we end up with x(i+1) being negative, then switch it
% to zero: 
        if x(i+1) < 0 
            x(i+1) = 0;
        end
% While the current error produced in the interval is not greater than
% error specified, continue!
        while active == 1
            xt = x(i):-precision:x(i+1); 
            y = zeros(1,length(xt));
            for j = 1:length(xt)
                y(j) = norminv(1-x(i))+(norminv(1-x(i+1))-...
                    norminv(1-x(i+1)))/(x(i+1)-x(i))*(xt(j)-norminv(1-x(i)));
                secderiv(j) = sqrt(2)*...
                    (2*pi*exp(2*erfinv(2*xt(j) - 1)^2)*...
                    erfinv(2*(1-xt(j)) - 1));
            end
            errorx =  abs(norminv(1-xt)-y);
            [max_error ind1] = max(errorx); 
            [max_sec ind2] = max(secderiv);
            max_err = 1/8*(x(i+1)-x(i))^2*max_sec;
            if max_error >= max_err && max_err >= error
                x(i+1) = xt(ind1-1);
                active = 1;
            else
                active = 0;
            end
            
        end

    
    fprintf('With point %i  x(i+1) = %1.4f\n',i,x(i+1))
    if isnan(x) == 1
        break;
    elseif x(i) == 0
        break;
    end
    i = i+1;

    end
    
    for indx_x = 1:length(x)
        y_2 = norminv(1-x(indx_x+1));
        y_1 = norminv(1-x(indx_x));
        x_2 = x(indx_x + 1);
        x_1 = x(indx_x);
        % y=mx+c where m and c are computed from secant end points
        cdf_approx_m(indx_x) = (y_2 - y_1)/(x_2 - x_1);
        cdf_approx_c(indx_x) = y_1 - cdf_approx_m(indx_x) * x_1;
    end 