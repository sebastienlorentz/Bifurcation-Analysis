% Time series of the 3-D SIRS model from a given initial condition: integrate the field
% and plot the states in idx against time. Reads the approach to an attractor, and the
% period directly off the trace when that attractor is a cycle.

lambda = 1.5;
mu = 0.1;
gamma = 3;
R0 = 1.29;
delta = 0.4; % waning immunity: R -> S

syms S I R

u = [S; I; R];
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0);     % f at the disease-free equilibrium (I = 0), = 1
beta = R0*(mu + gamma + f0); % normalized transmission (= raw beta * N0)

F = [mu - beta*S*I - mu*S + delta*R;
     beta*S*I - (mu + gamma + f)*I;
     (gamma + f)*I - (mu + delta)*R];

rhs = matlabFunction(F, 'Vars', {u}); % single n-vector field handle @(u)

u0 = [10 5 0]/15; % [S I R] fractions of N0, one row per initial condition
tspan = [0 500];

% idx picks which states to draw (here S, I); omit for all states
plot_time_series(rhs, u0, tspan, u, [1 2]);
