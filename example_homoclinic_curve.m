% Two-parameter homoclinic curve of the SIS model in (R0, gamma). A homoclinic orbit
% leaves the saddle along its unstable manifold and returns along the stable one; the
% locus of such connections is a codim-1 curve, the global bifurcation that terminates
% the family of limit cycles. Here we draw the fold/Hopf/transcritical curves and the
% fold-of-cycles (LPC) curve first, then seed the homoclinic BVP from a moderate-period
% orbit on the LPC arc and overlay the continued homoclinic curve in red.

lambda = 1.5;
mu = 0.1;

syms S I R0 gamma

u = [S; I];
p = [R0; gamma];
N0 = lambda/mu;         % disease-free population; states are fractions of it
f = 1/(1 + N0^2*I^2);   % treatment saturation carries the population scale
f0 = subs(f, I, 0);
beta = R0*(mu + gamma + f0);

F = [mu - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

pBox = [0.5 2; 0 6];
uStart = [1; 0];
stateBox = [0 3; 0 3];
invIdx = 2;

% --- codim-1 curves + codim-2 points -----------------------------------------
bif = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);
plot_bifurcation_curves(bif, {'R_0', '\gamma'}, pBox);

% --- fold-of-cycles (LPC) curve ----------------------------------------------
gGH = bif.gh_points(2, 1);
gammaSlices = gGH - [0.6 0.4 0.9 0.3 1.2 0.2];
gammaSlices = gammaSlices(gammaSlices > pBox(2, 1) & gammaSlices < pBox(2, 2));

lcsys = lpc_collocation_system(F, u, p);
[x0, T0, p0] = find_lpc_seed(F, u, p, gammaSlices, uStart, pBox(1, :), stateBox);
lpc = trace_lpc_collocation(lcsys, x0, T0, p0, pBox, 1e-2, 1000);
plot(lpc.p(1, :), lpc.p(2, :), 'Color', [0 0.6 0], 'LineWidth', 2, 'DisplayName', 'LPC (fold of cycles)');

% --- homoclinic curve ---------------------------------------------------------
% The LPC arc runs asymptotically alongside the homoclinic curve, so any far-tail
% orbit is a valid seed. Take a moderate period rather than the largest: at T ~ 200+
% almost all the arclength is the crawl past the saddle and a uniform mesh cannot
% resolve the fast excursion, so the landing Newton fails.
[~, k] = min(abs(lpc.T - 60));
hsys = homoclinic_collocation_system(F, u, p);
Xseed = reshape(lpc.W(:, k), hsys.n, hsys.Np);
% enough steps for the curve to reach the BT point and stop there on its own
hom = trace_homoclinic_curve(hsys, Xseed, lpc.T(k), lpc.p(:, k), pBox, 1e-2, 800);

% purple: red is the Hopf curve, green the LPC, so the homoclinic needs its own hue
plot(hom.p(1, :), hom.p(2, :), 'Color', [0.49 0.18 0.56], 'LineWidth', 2, 'DisplayName', 'homoclinic');
legend('Location', 'best');

% the saddle quantity lu + ls sets which side the cycles live on: negative means the
% bifurcating cycle is stable, positive unstable.
fprintf('homoclinic curve: %d points, R0 %.4f..%.4f, gamma %.4f..%.4f\n', ...
        numel(hom.T), min(hom.p(1, :)), max(hom.p(1, :)), min(hom.p(2, :)), max(hom.p(2, :)));
fprintf('truncation interval T: %.2f .. %.2f\n', min(hom.T), max(hom.T));
fprintf('saddle quantity: %+.5f .. %+.5f\n', min(hom.lu + hom.ls), max(hom.lu + hom.ls));
fprintf('arms ended: fwd = %s, bwd = %s\n', hom.fwd_reason, hom.bwd_reason);
