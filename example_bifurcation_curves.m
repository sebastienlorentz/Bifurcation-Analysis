% Two-parameter bifurcation diagram of the SIS model in (R0, gamma): the fold, Hopf and
% transcritical curves of the equilibria, with the Bogdanov-Takens and generalized Hopf
% points where they meet. Seeds come from one-parameter sweeps at fixed gamma slices.

lambda = 1.5;
mu = 0.1;

syms S I R0 gamma

u = [S; I];
p = [R0; gamma];
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0);     % f at the disease-free equilibrium (I = 0), = 1
beta = R0*(mu + gamma + f0); % normalized transmission (= raw beta * N0)

F = [mu - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

% the box is tall enough in gamma to capture the crest of the Hopf curve, where it meets
% the fold curve at a Bogdanov-Takens point, and wide enough in R0 to show the
% transcritical line at R0 = 1
pBox = [0.5 2; 0 6];   % [R0 range; gamma range]
uStart = [1; 0];       % disease-free equilibrium (normalized), valid on every slice
stateBox = [0 3; 0 3]; % physical region (fractions), stops nonphysical runaways
invIdx = 2;            % I is the invariant ("infected") coordinate (the I = 0 subspace)

out = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);
plot_bifurcation_curves(out, {'R_0', '\gamma'}, pBox);
