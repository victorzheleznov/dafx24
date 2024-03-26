% initialise strictly causal cubic spline in Meinsma et al. formulation
% input:
%   u1 --- initial value.
% output:
%   mem --- initial coefficients for pre-filter.
function mem = cspline2_init(u1)
    mem = [(4-sqrt(3))*u1; (3-sqrt(3))*u1];
end