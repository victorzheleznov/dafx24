% evaluate polynomial using Horner's rule
% input:
%   x --- input value or vector;
%   coefs --- polynomial coefficients (from the highest degree to the lowest).
% output:
%   y --- output.
function y = horner(x, coefs)
    y = coefs(1).*x;
    for i = 2:length(coefs)-1
        y = (y + coefs(i)).*x;
    end
    y = y + coefs(end);
end