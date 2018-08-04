function [cdf_approx_m, cdf_approx_c,varargout] = RolleLerp()
% Compute a piecewise-linear overapproximation of norminv(1-x) utilizing
% Rolle's theorem. 

% Define the second derivative of error function as a standin for norminv.
% This is done by evaluating the symbolic expression as follows: 

syms x; d = diff(erfinv(2*x-1),x,1)






end
