lambda = 1.5;
mu = 0.1;

syms S I R0 gamma

u = [S; I];
p = [R0; gamma];
f = 1/(1+I^2);
f0 = subs(f, I, 0); % f at the disease-free equilibrium (I = 0)
beta = R0*mu*(mu + gamma + f0)/lambda;

F = [lambda - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

% box is tall enough in gamma to capture the Hopf 'wave' (its crest is near
% gamma = 4, where the Hopf curve meets the fold curve at a Bogdanov-Takens
% point) and wide enough in R0 to show the transcritical line at R0 = 1
pBox = [0.5 2; 0 6]; % [R0 range; gamma range]
uStart = [15; 0]; % disease-free equilibrium, valid on every slice
stateBox = [0 50; 0 50]; % physical region, stops nonphysical runaways
invIdx = 2; % I is the invariant ("infected") coordinate (I = 0 subspace)

out = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);

plot_bifurcation_curves(out, {'R_0', '\gamma'}, pBox);
