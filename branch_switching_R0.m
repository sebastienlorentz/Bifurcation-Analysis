lambda = 1.5;
mu = 0.1;
gamma = 1;

syms S I R0

u = [S; I];
f = 1/(1+I^2);
f0 = subs(f, I, 0); % f at the disease-free equilibrium (I = 0)

% transmission rate expressed via the basic reproduction number:
% R0 = beta*lambda / (mu*(mu+gamma+f0))
beta = R0*mu*(mu + gamma + f0)/lambda;

F = [lambda - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

% build the handles (Jacobian + second derivatives) from the symbols
sys = continuation_system(F, u, R0);

% trace every branch, switching at branch points automatically. The state box
% stops a branch once it leaves the physical region (S, I >= 0); without it the
% nonphysical S -> inf, I -> -inf tail would run all the way to max_steps.
R0Range = [0, 4];
stateBox = [0 15; 0 15]; % [S; I] min/max
out = trace_branches(sys, [15; 0], R0Range(end), R0Range, 1e-4, 1e6, true, stateBox);

plot_branches(out, 'R_0', R0Range);
