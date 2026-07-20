% One-parameter bifurcation analysis of a 3-D SIRS model in R0, to exercise the
% dimension-general Hopf detection (det(2J (.) I) = 0). It extends the 2-D SIS
% driver with a recovered compartment R: susceptibles are infected, recover into
% R at the saturating rate gamma + f(I), and lose immunity back to S at rate
% delta. The extra compartment lifts the flow out of the plane, so the 2-D
% trace(J) = 0 Hopf test no longer applies.

lambda = 1.5;
mu = 0.1;
gamma = 1.5;
delta = 0.1; % waning immunity: R -> S

% these values sit just past a Generalized Hopf point: the branch has
% two supercritical Hopfs (R0 ~ 1.68 and 2.00) bounding a stable limit-cycle
% bubble. Stronger waning (e.g. delta = 0.4) stabilises the focus and leaves no
% Hopf at all, so nothing bifurcates.

syms S I R R0

u = [S; I; R];
f = 1/(1+I^2); % saturating recovery (limited treatment capacity)
f0 = subs(f, I, 0); % f at the disease-free equilibrium (I = 0)

% transmission rate expressed via the basic reproduction number
beta = R0*mu*(mu + gamma + f0)/lambda;

F = [lambda - beta*S*I - mu*S + delta*R;
     beta*S*I - (mu + gamma + f)*I;
     (gamma + f)*I - (mu + delta)*R];

% build the handles (Jacobian + second/third derivatives) from the symbols
sys = continuation_system(F, u, R0);

% trace every branch, switching at branch points automatically. The state box
% stops a branch once it leaves the physical region (S, I, R >= 0).
R0Range = [0, 4];
stateBox = [0 20; 0 20; 0 20]; % [S; I; R] min/max
out = trace_branches(sys, [15; 0; 0], R0Range(end), R0Range, 1e-4, 1e6, true, stateBox);

% plot_branches draws the state components listed in the index vector; here all
% three (S, I, R) vs R0. Stability, folds, Hopfs (styled by criticality) and
% branch points are all marked. Pass e.g. [1 3] to draw only S and R.
plot_branches(out, 'R_0', R0Range, [1 2 3], {'S', 'I', 'R'});

% every Hopf spawns a branch of periodic orbits; continue them all by single
% shooting and overlay their S/I/R max/min envelopes (same idx as plot_branches)
% on the diagram above, plus a companion period figure. The two Hopfs bound the
% same cycle bubble, so 150 continuation steps per side already cover it fully;
% the default (1e4) would grind on long past the bubble. Args after out are
% forwarded to trace_limit_cycle as (dlambda, step_size, max_steps).
lcs = trace_limit_cycles(sys, out, [], [], 150);
plot_limit_cycles(lcs, 'R_0', [1 2 3], {'S', 'I', 'R'});

% phase portrait at a fixed R0, projected onto 2 of the 3 states. The off-plane
% state is held at the in-window equilibrium when the field is drawn; trajectories
% integrate in full 3-D. Here the (S, I) plane; pass e.g. [2 3] for (I, R).
R0fixed = 1.8; % mid-bubble, where the stable cycle is largest
Fnum = subs(F, R0, R0fixed);
phase_portrait(Fnum, u, stateBox, [1 2], {'S', 'I', 'R'});
phase_portrait(Fnum, u, stateBox, [2 3], {'S', 'I', 'R'});
