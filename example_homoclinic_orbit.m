% A single homoclinic orbit of the SIS model in the (S, I) phase plane: the loop, the
% saddle q it leaves and returns to, the ball of radius eps the boundary conditions pin
% the two ends to, and the equilibrium it encircles. The seeding chain is the same as
% example_homoclinic_curve (bifurcation curves, LPC arc, moderate-period orbit); only a
% few continuation steps are taken here, since one orbit is wanted rather than a curve.
%
% Writes docs/figures/homoclinic_orbit.pdf.

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

pBox = [0.5 2; 0 6];   % [R0 range; gamma range]
uStart = [1; 0];       % disease-free equilibrium (normalized), valid on every slice
stateBox = [0 3; 0 3]; % physical region (fractions), stops nonphysical runaways
invIdx = 2;            % I is the invariant ("infected") coordinate

% --- seed the homoclinic off the fold-of-cycles arc ---------------------------
bif = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);

gGH = bif.gh_points(2, 1);
gammaSlices = gGH - [0.6 0.4 0.9 0.3 1.2 0.2];
gammaSlices = gammaSlices(gammaSlices > pBox(2, 1) & gammaSlices < pBox(2, 2));

lcsys = lpc_collocation_system(F, u, p);
[x0, T0, p0] = find_lpc_seed(F, u, p, gammaSlices, uStart, pBox(1, :), stateBox);
lpc = trace_lpc_collocation(lcsys, x0, T0, p0, pBox, 1e-2, 1000);

% moderate period, not the largest available: at T ~ 200+ almost all the arclength is
% the crawl past the saddle and a uniform mesh cannot resolve the fast excursion
[~, k] = min(abs(lpc.T - 60));
hsys = homoclinic_collocation_system(F, u, p);
Xseed = reshape(lpc.W(:, k), hsys.n, hsys.Np);
hom = trace_homoclinic_curve(hsys, Xseed, lpc.T(k), lpc.p(:, k), pBox, 1e-2, 20);

% --- pick the clearest loop of the few traced ---------------------------------
n = hsys.n; Np = hsys.Np;
amps = arrayfun(@(j) max(vecnorm(reshape(hom.X(:, j), n, Np) - hom.q(:, j))), ...
                1:size(hom.X, 2));
[amp, j] = max(amps);

X   = reshape(hom.X(:, j), n, Np);
q   = hom.q(:, j);
pj  = hom.p(:, j);
e   = hom.eps(j);
circ = q + e*[cos(linspace(0, 2*pi, 200)); sin(linspace(0, 2*pi, 200))];

% the equilibrium the loop encircles, found from the loop's centroid
Fn = matlabFunction(F, 'Vars', {u, p});
focus = fsolve(@(z) Fn(z(:), pj), mean(X, 2), optimset('Display', 'off'));
have_focus = all(focus > min(X, [], 2)) && all(focus < max(X, [], 2)) && norm(focus - q) > 1e-3;

fprintf('R0 = %.4f, gamma = %.4f, T = %.1f\n', pj(1), pj(2), hom.T(j));
fprintf('saddle q = (%.4f, %.4f), lu = %+.4f, ls = %+.4f\n', q(1), q(2), hom.lu(j), hom.ls(j));
fprintf('eps = %.3e, loop amplitude = %.3e, ratio = %.3f\n', e, amp, e/amp);

% --- draw ---------------------------------------------------------------------
purple = [0.49 0.18 0.56];   % the homoclinic colour used throughout
grey   = [0.2 0.2 0.2];

% Pin the canvas: a bare figure inherits the last window's size, and since
% exportgraphics crops tight, the PDF's physical size then varies run to run. LaTeX
% scales it to a fixed width, so that variation shows up as the fonts changing size.
figure('Units', 'inches', 'Position', [1 1 10 5]);
hold on; grid on;

plot(X(1, :), X(2, :), 'Color', purple, 'LineWidth', 1.8, 'DisplayName', 'Homoclinic orbit');

% an arrowhead a quarter of the way round, to give the loop a direction
kq = round(0.25*Np);
d = X(:, kq+1) - X(:, kq);
d = 0.02*d/norm(d);          % fixed short arrow, not scaled by the step size
quiver(X(1, kq), X(2, kq), d(1), d(2), 0, 'Color', purple, 'LineWidth', 1.5, ...
       'MaxHeadSize', 3, 'HandleVisibility', 'off');

plot(circ(1, :), circ(2, :), ':', 'Color', grey, 'LineWidth', 1, ...
     'DisplayName', sprintf('Truncation ball, \\epsilon = %.2e', e));

% the untruncated orbit takes infinite time over these last two stretches, so the BVP
% cannot represent them: sketch them in to show the loop does close on q
plot([q(1) X(1, 1)], [q(2) X(2, 1)], '--', 'Color', purple, 'LineWidth', 1, ...
     'DisplayName', 'Asymptotic tails (t \rightarrow \pm\infty)');
plot([X(1, end) q(1)], [X(2, end) q(2)], '--', 'Color', purple, 'LineWidth', 1, ...
     'HandleVisibility', 'off');

plot(q(1), q(2), 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 8, ...
     'LineStyle', 'none', 'DisplayName', 'Saddle');
if have_focus
    plot(focus(1), focus(2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6, ...
         'LineStyle', 'none', 'DisplayName', 'Endemic equilibrium');
end

xlabel('S'); ylabel('I');
title(sprintf('Homoclinic orbit at R_0 = %.3f, \\gamma = %.3f', pj(1), pj(2)), 'FontSize', 14);
legend('Location', 'northeast', 'FontSize', 12);
axis equal; axis tight; box on; % equal aspect so the truncation ball reads as a circle
xlim(xlim + 0.04*diff(xlim)*[-1 1]); ylim(ylim + 0.15*diff(ylim)*[-1 1]);

export_figure('homoclinic_orbit');
