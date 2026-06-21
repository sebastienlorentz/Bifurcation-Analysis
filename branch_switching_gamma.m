lambda = 1.5;
mu = 0.1;
beta = 0.17;

syms S I gamma

u = [S; I];
f = 1/(1+I^2);
F = [lambda - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

% build the handles (Jacobian + second derivatives) from the symbols
sys = continuation_system(F, u, gamma);

% trace every branch, switching at branch points automatically. The state box
% stops a branch once it leaves the physical region (S, I >= 0); without it the
% nonphysical I < 0 tail would run on needlessly.
% [15; 0] is the disease-free state (S = lambda/mu = 15, I = 0), valid for any gamma.
gammaRange = [0, 3];
stateBox = [0 15; 0 15]; % [S; I] min/max
out = trace_branches(sys, [15; 0], gammaRange(end), gammaRange, 1e-4, 1e6, true, stateBox);

plot_branches(out, '\gamma', gammaRange);
