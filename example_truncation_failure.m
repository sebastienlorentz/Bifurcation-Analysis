% What goes wrong when the homoclinic truncation distance eps is held constant instead
% of pinned to the loop amplitude. The homoclinic curve ends at a Bogdanov-Takens point,
% where the loop shrinks to nothing; a constant eps then swallows the whole orbit, and
% the continuation sails past BT still converging to tolerance. This runs one seed twice,
% once with the constant eps and once with the shipped system, and draws both curves plus
% the eps/amplitude ratio that diagnoses the difference.
%
% The broken run is not a copy of the old code: the fix was a single residual row, so
% that one row is overwritten on the current system (see broken_pin below) and cannot
% drift away from what ships.
%
% Writes three PDFs to docs/figures/, assembled into one figure in the report.

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

% --- the seed, and the BT point both runs are measured against -----------------
% bt comes off the Hopf curve where omega passes through zero, so it is independent of
% either homoclinic run and can be used to judge them.
bif = trace_bifurcation_curves(F, u, p, pBox, uStart, stateBox, 7, invIdx);
bt = bif.bt_points(:, 1);

gGH = bif.gh_points(2, 1);
gammaSlices = gGH - [0.6 0.4 0.9 0.3 1.2 0.2];
gammaSlices = gammaSlices(gammaSlices > pBox(2, 1) & gammaSlices < pBox(2, 2));

lcsys = lpc_collocation_system(F, u, p);
[x0, T0, p0] = find_lpc_seed(F, u, p, gammaSlices, uStart, pBox(1, :), stateBox);
lpc = trace_lpc_collocation(lcsys, x0, T0, p0, pBox, 1e-2, 1000);

[~, k] = min(abs(lpc.T - 60));  % moderate period seed, as in example_homoclinic_curve
hsys = homoclinic_collocation_system(F, u, p);
ix = hsys.idx; n = hsys.n; Np = hsys.Np;
Xseed = reshape(lpc.W(:, k), n, Np);
Tseed = lpc.T(k);
pseed = lpc.p(:, k);

% --- after the fix: the shipped system ----------------------------------------
fixed = trace_homoclinic_curve(hsys, Xseed, Tseed, pseed, pBox, 1e-2, 800);

% --- before the fix: eps frozen at its landed value ---------------------------
% Landing was never the bug, so land with the correct system, then freeze eps at the
% value it landed on and continue with no bracket on T, so the blow-up shows in full.
[Z0, delta, kfar] = hsys.seed_from_cycle(Xseed, Tseed, pseed, []);
prob = hsys.make_problem(reshape(Z0(ix.X), n, Np), delta, kfar);
Z0 = land(prob, Z0, ix);
assert(norm(prob.residual(Z0)) < 1e-8, 'landed point does not satisfy the homoclinic BVP');

pts = follow_curve(broken_pin(prob, ix, Z0(ix.e)), Z0, pBox, 1e-2, 800);
broke.p = pts(ix.p, :);
broke.e = pts(ix.e, :);
broke.X = pts(ix.X, :);
broke.q = pts(ix.q, :);

% --- the diagnostic: eps against the loop amplitude it is supposed to truncate --
amp = @(X, q) arrayfun(@(j) max(vecnorm(reshape(X(:, j), n, Np) - q(:, j))), 1:size(X, 2));
r_bad = broke.e ./ amp(broke.X, broke.q);
r_fix = fixed.eps ./ amp(fixed.X, fixed.q);

fprintf('BT point (independent):     R0 = %.4f, gamma = %.4f\n', bt(1), bt(2));
fprintf('closest approach, before:   %.4f  (%d points)\n', min(vecnorm(broke.p - bt)), size(broke.p, 2));
fprintf('closest approach, after:    %.4f  (%d points)\n', min(vecnorm(fixed.p - bt)), size(fixed.p, 2));
fprintf('eps/amplitude, before: %.3f .. %.3f\n', min(r_bad), max(r_bad));
fprintf('eps/amplitude, after:  %.3f .. %.3f\n', min(r_fix), max(r_fix));

% --- the two curve panels ------------------------------------------------------
% Colour carries the before/after distinction, not line style: dashed means
% unstable/subcritical everywhere else in this project.
grey   = [0.4 0.4 0.4];
purple = [0.49 0.18 0.56];
win_x = bt(1) + [-0.12 0.12];
win_y = bt(2) + [-0.6  0.6];

% Titles carry no numbers; the closest approaches are printed above and quoted in the
% surrounding prose, and repeating them crowds an already small panel.
curve_panel(bif, broke.p, grey, win_x, win_y, 'Before the fix: \epsilon constant');
export_figure('failure_homoclinic_curve_before');

curve_panel(bif, fixed.p, purple, win_x, win_y, 'After the fix: \epsilon = \delta \times amplitude');
export_figure('failure_homoclinic_curve_after');

% --- the ratio panel -----------------------------------------------------------
% Truncation is only honest while the ball is a small dimple on a large loop, so this
% quantity has an absolute ceiling at 1: the ball is the whole orbit and there is
% nothing left outside it to have truncated.
%
% x is arclength along each curve, normalised per curve. Point index is arbitrary and
% hides where BT sits, gamma is not monotone and folds back, and raw shared arclength
% leaves the after-fix line ending at 45% of the width, which reads as a broken plot.
% The consequence is that reading across at a fixed x compares different points, so the
% BT markers anchor the two lines. Both runs share one axes: the content here is the
% contrast between them, and splitting the pair left the after panel empty.
s_bad = arclen(broke.p); s_bad = s_bad/s_bad(end);
s_fix = arclen(fixed.p); s_fix = s_fix/s_fix(end);
[~, kb] = min(vecnorm(broke.p - bt));
[~, kf] = min(vecnorm(fixed.p - bt));

% same 7in canvas width as the curve panels: LaTeX scales by width, so equal canvas
% width plus equal font size is what makes the three panels match on the page
figure('Units', 'inches', 'Position', [1 1 7 4.5]);
hold on; grid on;
plot(s_bad, r_bad, 'Color', grey, 'LineWidth', 1.5, 'DisplayName', 'Before the fix');
plot(s_fix, r_fix, 'Color', purple, 'LineWidth', 1.5, 'DisplayName', 'After the fix');
yline(1, 'k:', 'LineWidth', 1.2, 'DisplayName', 'ball covers the whole orbit');
plot(s_bad(kb), r_bad(kb), 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 7, ...
     'LineStyle', 'none', 'DisplayName', 'closest approach to BT');
plot(s_fix(kf), r_fix(kf), 'ks', 'MarkerFaceColor', 'm', 'MarkerSize', 7, ...
     'LineStyle', 'none', 'HandleVisibility', 'off');
xlim([0 1]); ylim([0 1.08]);
xlabel('Fraction of the curve traced (by arclength)');
ylabel('\epsilon / loop amplitude');
% two lines at 20pt: on one line the title is wider than the 7in canvas and would
% squeeze the axes to make room
title({'Fraction of the orbit swallowed', 'by the truncation ball'}, 'FontSize', 20);
legend('Location', 'northwest', 'FontSize', 16);
export_figure('failure_homoclinic_ratio');

% Cumulative arclength along a traced curve in the parameter plane.
function s = arclen(P)
    s = [0, cumsum(vecnorm(diff(P, 1, 2)))];
end

% One homoclinic curve on its fold/Hopf/BT context, zoomed and centred on BT.
% plot_bifurcation_curves draws exactly this context but opens its own figure over the
% whole box, so the context is redrawn here in the same colours.
function curve_panel(bif, P, col, win_x, win_y, ttl)
    % pinned canvas, as in example_homoclinic_orbit: exportgraphics crops tight, so an
    % inherited window size would drift and LaTeX's width scaling would show that up as
    % the fonts changing size
    figure('Units', 'inches', 'Position', [1 1 7 5.5]);
    hold on; grid on;
    for i = 1:numel(bif.fold_curves)
        c = bif.fold_curves{i};
        a = legend_arg('Fold', i == 1);
        plot(c.p(1, :), c.p(2, :), 'k', 'LineWidth', 1.5, a{:});
    end
    for i = 1:numel(bif.hopf_curves)
        c = bif.hopf_curves{i};
        [Xs, Xu] = split_stable(c.p(1, :), c.l1 < 0);
        [Ys, Yu] = split_stable(c.p(2, :), c.l1 < 0);
        a = legend_arg('Hopf (supercritical)', i == 1);
        plot(Xs, Ys, 'r', 'LineWidth', 1.5, a{:});
        a = legend_arg('Hopf (subcritical)', i == 1);
        plot(Xu, Yu, 'r--', 'LineWidth', 1.5, a{:});
    end
    plot(bif.bt_points(1, :), bif.bt_points(2, :), 'ks', 'MarkerFaceColor', 'm', ...
         'MarkerSize', 7, 'LineStyle', 'none', 'DisplayName', 'Bogdanov-Takens');
    plot(P(1, :), P(2, :), 'Color', col, 'LineWidth', 1.5, 'DisplayName', 'Homoclinic');
    xlim(win_x); ylim(win_y);
    xlabel('R_0'); ylabel('\gamma');
    title(ttl, 'FontSize', 20);
    legend('Location', 'southeast', 'FontSize', 16);
end

% Legend arguments for a set of curves sharing one entry: name the first, hide the rest.
function a = legend_arg(name, first)
    if first
        a = {'DisplayName', name};
    else
        a = {'HandleVisibility', 'off'};
    end
end

% Replace the truncation-scale row, the last residual row of the homoclinic system,
% with the pre-fix e - e0. Everything else is untouched.
function bad = broken_pin(prob, ix, e0)
    bad.residual = @(Z) set_last(prob.residual(Z), Z(ix.e) - e0);
    bad.jacobian = @(Z) set_last_row(prob.jacobian(Z), ix);
end

function R = set_last(R, v)
    R(end) = v;
end

function J = set_last_row(J, ix)
    J(end, :) = 0;
    J(end, ix.e) = 1;
end

% Damped Newton onto the homoclinic BVP with p1 pinned, keeping T away from zero. Same
% landing trace_homoclinic_curve does internally; it is private to that file, and the
% assert above catches any drift between the two.
function Z = land(prob, Z0, ix)
    Z = Z0;
    T0 = Z0(ix.T);
    pin = zeros(1, ix.N); pin(ix.p(1)) = 1;
    for it = 1:60
        R = [prob.residual(Z); Z(ix.p(1)) - Z0(ix.p(1))];
        if norm(R) < 1e-9; return; end
        dZ = [prob.jacobian(Z); pin]\R;
        s = 1;
        while Z(ix.T) - s*dZ(ix.T) < 0.25*T0 && s > 1e-3; s = s/2; end
        Z = Z - s*dZ;
    end
    error('example_truncation_failure:land', 'could not land (|R| = %.2e)', norm(R));
end
