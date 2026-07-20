% Split a data series into its stable (solid) and unstable (dashed) parts for
% plotting, by masking the other part out with NaN. The stable mask is extended
% one point each side so the segment spanning a stability change (a Hopf or fold)
% is still drawn solid, rather than dropped from both parts and left as a gap.
% Shared by the equilibrium, limit-cycle and Hopf-curve plotters.
function [solid, dashed] = split_stable(y, stable)
    stx = stable | [stable(2:end), false] | [false, stable(1:end-1)];
    solid = y; solid(~stx) = NaN;
    dashed = y; dashed(stable) = NaN;
end
