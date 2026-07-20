% Phase portrait of the SIS model at fixed parameters: a speed heatmap with a
% unit direction field, and forward/backward trajectories, both marking the
% classified equilibria. phase_portrait is dimension-general; for a system with
% 3+ states pass idx to pick the 2 to plot (e.g. [2 3]) and it projects onto them.

lambda = 1.5;
mu = 0.1;
gamma = 3;
R0 = 1.29;

syms S I
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0); % treatment term at the disease-free equilibrium (I = 0), = 1
beta = R0*(mu + gamma + f0); % normalized transmission (= raw beta * N0)
F = [mu - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

u = [S; I];
box = [0 1.5; 0 1.5]; % [S; I] min/max (fractions)

% draw both figures; idx defaults to [1 2] and state_labels to the symbol names
phase_portrait(F, u, box, [1 2], {'S', 'I'});
