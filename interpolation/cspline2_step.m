% calculate step of strictly causal cubic spline in Meinsma et al. formulation [1]
% input:
%   u --- next sample;
%   u1 --- current sample;
%   fs --- sampling rate [Hz];
%   mem --- memory bus (previous pre-filter output).
% output:
%   coefs --- polynomial coefficients (from the highest degree to the lowest);
%   mem --- updated memory bus.
% references:
% [1] G. Meinsma and L. Mirkin, "L2 Sampled signal reconstruction with 
%     causality constraints - Part I: Setup and solutions," in IEEE 
%     Transactions on Signal Processing, vol. 60, no. 5, pp. 2260-2272, 
%     May 2012, doi: 10.1109/TSP.2012.2185228.
function [coefs, mem] = cspline2_step(u, u1, fs, mem)
    %% parameters
    alpha = sqrt(3) - 2;
    w11_1 = mem(1);
    w12_1 = mem(2);

    %% prefilter
    % w1 (future by 1 sample)
    w11 = alpha*w11_1 + (4-sqrt(3))*u + alpha*(2+sqrt(3))*u1;
    w12 = alpha*w12_1 + (3-sqrt(3))*u + alpha*(3+sqrt(3))*u1;
    
    % w2 (current)
    w21_1 = 6*fs^3*(-4+3*sqrt(3))*u + ...
            alpha*fs^3*(6*sqrt(3))*w11 + ...
            alpha*fs^3*(3-3*sqrt(3))*w12 + ...
            fs^3*(6-6*sqrt(3))*w11_1 + ...
            fs^3*(3*sqrt(3)-3)*w12_1;
    w22_1 = 0;

    %% output
    % current (curve between u1 and u)
    coefs = [-w21_1/6,... 
            (w21_1 + w22_1)/(2*fs),...
            (w21_1 + w22_1)/(2*sqrt(3)*fs^2) + w12_1*fs,...
             w11_1 - w12_1];

    %% write to memory
    mem = [w11; w12];
end