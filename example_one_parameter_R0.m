% One-parameter bifurcation analysis of the SIS model in R0: the equilibrium
% branches (with stability, folds, Hopfs and branch points), then the limit
% cycles born at each Hopf, continued by single shooting and overlaid on the
% same diagram. 

lambda = 1.5;
mu = 0.1;
gamma = 3;

syms S I R0

u = [S; I];
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0);     % f at the disease-free equilibrium (I = 0), = 1

% normalized transmission via the basic reproduction number (= raw beta * N0)
beta = R0*(mu + gamma + f0);

F = [mu - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

% build the handles (Jacobian + second derivatives) from the symbols
sys = continuation_system(F, u, R0);

% trace every branch, switching at branch points automatically. The state box
% stops a branch once it leaves the physical region (S, I >= 0); without it the
% nonphysical S -> inf, I -> -inf tail would run all the way to max_steps.
R0Range = [0, 4];
stateBox = [0 3; 0 1]; % [S; I] min/max (normalized: disease-free s = 1)
out = trace_branches(sys, [1; 0], R0Range(end), R0Range, 1e-4, 1e6, true, stateBox);

plot_branches(out, 'R_0', R0Range, [1 2], {'S', 'I'});

% every Hopf spawns a branch of periodic orbits; continue them all by single
% shooting (seeded from the Hopfs in out) and overlay their S/I max/min
% envelopes on the diagram above, plus a companion period figure. Each branch
% stops itself at the Hopf and at the homoclinic approach.
lcs = trace_limit_cycles(sys, out);
plot_limit_cycles(lcs, 'R_0', [1 2], {'S', 'I'}); % same idx as plot_branches
