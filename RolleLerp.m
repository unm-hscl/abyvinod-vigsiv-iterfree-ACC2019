function [x] = RolleLerp()
% Compute a piecewise-linear overapproximation of norminv(1-x) utilizing
% Rolle's theorem. 

% Define the second derivative of error function as a standin for norminv.
% initialize x0 and x1 along with iteration, i, and the precision to move
% along: 
    x(1) = 0.5;
    x(2) = 0.4999;
    i = 1;
    precision = 1E-6; 
    error = 1E-4;
% while the resulting values do not produce obsure values, continue the PW
% approximation: 
    while erfinv(x(i)) ~= Inf || isnan(erfinv(x(i)))~= 1
% while the current error produced in the interval is not greater than
% 10E-4 continue!
        while error >= max_sec && max_sec >= error
            x(i+1) = xt(ind);
            xt = x(i):-precision:x(i+1); 
            for j = 1:length(xt)
                y(j) = norminv(1-x(i))+(norminv(1-x(i+1))-...
                    norminv(1-x(i+1)))/(x(i+1)-x(i))*(xt(j)-norminv(1-x(i)));
                secderiv(j) = sqrt(2)*...
                    (2*pi*exp(2*erfinv(2*(1-j) - 1)^2)*...
                    erfinv(2*(1-j) - 1));
            end
            error =  abs(norminv(1-xt)-y); 
            [max_sec ind] = max(secderiv);
            max_err = 1/8*(x(i+1)-x(i))^2*max_sec;
            
        end

    i = i+1;
    if isnan(x) == 1
        break;
    end
    end

end
