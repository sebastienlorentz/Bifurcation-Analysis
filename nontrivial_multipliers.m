% Floquet multipliers with the trivial one removed. A periodic orbit always carries
% a multiplier at +1 (the flow direction); the rest govern stability, and a real one
% crossing +1 is a fold of cycles. The trivial one is identified as the multiplier
% nearest 1, which is also the accuracy monitor: its drift away from 1 bounds where
% the variational integration is still trustworthy.
function nt = nontrivial_multipliers(mu)
    [~, it] = min(abs(mu - 1));
    nt = mu([1:it-1, it+1:end]);
end
