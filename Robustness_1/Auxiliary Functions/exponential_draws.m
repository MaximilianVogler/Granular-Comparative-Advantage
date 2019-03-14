function [R] = exponential_draws(sigmaT,rudraws,rvdraws)
% This function generates draws from a two-sided exponential distribution.

Finv = -sigmaT*log(1-rudraws);

if Finv <0
    error('The inverse of F needs to be non-negative.')
end

R = 1/(sqrt(2))*sign(rvdraws-0.5).*Finv;

end