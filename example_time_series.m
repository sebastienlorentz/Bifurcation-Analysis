% Time series of the 3-D SIRS model from a given initial condition: integrate the
% field and plot the states in idx against time. The period (and its divergence
% near a global bifurcation) is visible directly. A stepping stone toward the
% global-bifurcation goal.

lambda = 1.5;
mu = 0.1;
gamma = 3;
R0 = 1.29;
delta = 0.4; % waning immunity: R -> S

syms S I R

u = [S; I; R];
f = 1/(1+I^2); % saturating recovery (limited treatment capacity)
f0 = subs(f, I, 0); % f at the disease-free equilibrium (I = 0)

% transmission rate expressed via the basic reproduction number
beta = R0*mu*(mu + gamma + f0)/lambda;

F = [lambda - beta*S*I - mu*S + delta*R;
     beta*S*I - (mu + gamma + f)*I;
     (gamma + f)*I - (mu + delta)*R];

rhs = matlabFunction(F, 'Vars', {u}); % single n-vector field handle @(u)

u0 = [10 5 0]; % [S I R], one row per initial condition
tspan = [0 500];

% idx picks which states to draw (here S, I); omit for all states
plot_time_series(rhs, u0, tspan, u, [1 2]);
