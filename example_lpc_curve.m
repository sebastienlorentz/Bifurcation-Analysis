% Two-parameter fold-of-cycles (LPC) curve of the SIS model in (R0, gamma). The
% LPC is the locus of folds of limit cycles: periodic orbits whose 
% Floquet multiplier is +1. It is born at the generalized Hopf point
% where the Hopf curve switches super/subcritical, and bounds the region of
% bistability between a stable cycle and a stable equilibrium. Here we draw the
% fold/Hopf/transcritical curves first (with the Bautin point), then seed an LPC
% from a fold of cycles on a one-parameter limit-cycle branch and overlay the
% continued LPC curve in green.

lambda = 1.5;
mu = 0.1;

syms S I R0 gamma

u = [S; I];
p = [R0; gamma];
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0); % f at the disease-free equilibrium (I = 0), = 1
beta = R0*(mu + gamma + f0); % normalized transmission (= raw beta * N0)

F = [mu - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

pBox = [0.5 2; 0 6];   % [R0 range; gamma range]
uStart = [1; 0];       % disease-free equilibrium (normalized), valid on every slice
stateBox = [0 3; 0 3]; % physical region (fractions), stops nonphysical runaways
invIdx = 2;            % I is the invariant ("infected") coordinate

% --- codim-1 curves + codim-2 points -----------------------------------------
bif = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);
plot_bifurcation_curves(bif, {'R_0', '\gamma'}, pBox);

% --- fold-of-cycles (LPC) curve ----------------------------------------------
% the fold of cycles is born at the GH point, but is fragile to seed right
% there: the cycle amplitude vanishes at the GH, so the fold sits at tiny
% amplitude the shooting corrector collapses back onto the equilibrium. Instead
% seed from robust slices well inside the subcritical region (below the GH
% gamma), where the fold sits at a healthy amplitude, and let the continuation
% carry the curve back up to the GH point. Falls back to a spread of slices
% if no generalized Hopf point was detected.
if ~isempty(bif.gh_points)
    gGH = bif.gh_points(2, 1);
    gammaSlices = gGH - [0.6 0.4 0.9 0.3 1.2 0.2]; % subcritical side, robust first
else
    gammaSlices = linspace(pBox(2, 1) + 0.5, pBox(2, 2) - 0.5, 11);
end
gammaSlices = gammaSlices(gammaSlices > pBox(2, 1) & gammaSlices < pBox(2, 2));

lcsys = lpc_collocation_system(F, u, p);
[x0, T0, p0] = find_lpc_seed(F, u, p, gammaSlices, uStart, pBox(1, :), stateBox);
% The curve is born at the GH point (top) and runs down the subcritical side.
% Orthogonal collocation follows the whole fold-of-cycles curve. Unlike single
% shooting it stays robust as the cycle stretches toward a homoclinic (T -> inf),
% so the arc runs well past where shooting stalled, at a fraction of the cost.
lpc = trace_lpc_collocation(lcsys, x0, T0, p0, pBox, 1e-2, 1000);

plot(lpc.p(1, :), lpc.p(2, :), 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', 'LPC (fold of cycles)');
legend('Location', 'best');
