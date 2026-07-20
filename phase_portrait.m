% Phase portrait of an n-D system at fixed parameters, projected onto the two
% states in idx. Draws two figures: a speed heatmap with a unit direction field,
% and a set of forward/backward trajectories with a coarse direction field; both
% mark the classified equilibria. F_sym is the symbolic vector field (n-vector),
% u_sym its states, box an n-by-2 of per-state [min max]. idx (optional) names the
% 2 states to plot (default [1 2]); state_labels (optional) names all n states.
% For n > 2 the off-plane states are held at the in-window equilibrium nearest the
% centre of the plotted box when the 2-D field is drawn, so the slice passes
% through a fixed point; trajectories still integrate in full n-D and are projected.
function phase_portrait(F_sym, u_sym, box, idx, state_labels)
    u_sym = u_sym(:);
    n = numel(u_sym);
    if nargin < 4 || isempty(idx)
        idx = [1 2];
    end
    assert(numel(idx) == 2, 'phase_portrait draws a 2-D projection: idx must name exactly 2 states.');
    if nargin < 5 || isempty(state_labels)
        state_labels = arrayfun(@(k) char(u_sym(k)), 1:n, 'UniformOutput', false);
    end

    xrange = box(idx(1), :);
    yrange = box(idx(2), :);
    labels = state_labels(idx);

    % equilibria (full n-D) plus the reference point whose off-plane coordinates
    % fix the slice the 2-D field is drawn on
    [eqpts, cls] = equilibria(F_sym, u_sym, box);
    eq0 = slice_reference(eqpts, box, idx);

    % component handles of all n states, so the plotted plane can be swept with
    % the off-plane coordinates held at eq0
    fi = matlabFunction(F_sym(idx(1)), 'Vars', num2cell(u_sym).');
    fj = matlabFunction(F_sym(idx(2)), 'Vars', num2cell(u_sym).');
    rhs = matlabFunction(F_sym, 'Vars', {u_sym});

    % integration horizon for the trajectories: scale to the slowest linear decay
    % at the equilibria, so trajectories actually reach the attractors instead of
    % stalling on a slow manifold (e.g. a stiff system with fast infection and
    % slow demography settles orders of magnitude later than a fixed span would show)
    [Tf, Tb] = trajectory_time(F_sym, u_sym, eqpts);

    field_portrait(fi, fj, eq0, idx, eqpts, cls, xrange, yrange, labels);
    trajectory_portrait(rhs, fi, fj, eq0, idx, eqpts, cls, box, xrange, yrange, labels, Tf, Tb);
end

% Forward/backward integration horizons from the equilibrium eigenvalues.
% Forward runs ~8 of the slowest timescale so trajectories reach the attractors.
% Backward runs only ~8 of the fastest timescale: in a stiff system the fast
% mode blows up in reverse time, so a long backward span just produces runaway
% lines; a short one still sketches the local outflow near repellers/saddles.
function [Tf, Tb] = trajectory_time(F_sym, u_sym, eqpts)
    if isempty(eqpts)
        Tf = 500; Tb = 500;
        return
    end
    Jfun = matlabFunction(jacobian(F_sym, u_sym), 'Vars', {u_sym});
    slow = inf; fast = 0;
    for r = 1:size(eqpts, 1)
        re = abs(real(eig(Jfun(eqpts(r, :).'))));
        re = re(re > 1e-9); % ignore zero modes (centres, conserved directions)
        if ~isempty(re)
            slow = min(slow, min(re));
            fast = max(fast, max(re));
        end
    end
    if ~isfinite(slow)
        Tf = 500;
    else
        Tf = min(8/slow, 1e6); % ~8 slow times, capped so a near-zero rate can't run away
    end
    if fast == 0
        Tb = Tf;
    else
        Tb = min(2/fast, Tf); % ~2 fast times: enough to sketch local outflow,
                              % short enough not to paint the fast foliation
    end
end

function field_portrait(fi, fj, eq0, idx, eqpts, cls, xrange, yrange, labels)
    [X, Y, U, V, speed] = plane_field(fi, fj, eq0, idx, xrange, yrange, 150);

    figure;
    hold on;

    % speed heatmap: log scale so the slow regions near equilibria read as
    % distinct dark patches instead of being crushed by the fast corners
    logspeed = log10(speed + eps);
    pcolor(X, Y, logspeed);
    colormap(turbo)
    shading interp;
    cb = colorbar;
    cb.Label.String = 'log_{10} speed';
    % clamp to a few decades below the fastest speed, so a single near-zero
    % grid point (an equilibrium landing on a node) can't blow out the scale
    hi = max(logspeed(:));
    clim([hi - 4, hi]);

    % unit-length arrows (direction only); the heatmap carries the magnitude
    arrow_skip = 5; % draw every arrow_skip-th arrow so they don't swamp the map
    s = 1:arrow_skip:size(X, 1);
    quiver(X(s, s), Y(s, s), U(s, s), V(s, s), 0.2, 'w');

    draw_equilibria(eqpts, cls, idx, 4);
    finalize_axes(labels, xrange, yrange, 'Phase portrait');
end

function trajectory_portrait(rhs, fi, fj, eq0, idx, eqpts, cls, box, xrange, yrange, labels, Tf, Tb)
    odefun = @(t, y) rhs(y);

    figure;
    hold on;

    % seed a grid of initial conditions on the plotted plane (off-plane states
    % held at eq0) and integrate both ways in time; trajectories flow into
    % attractors (forward, long horizon) and away from repellers (backward, short)
    nseed = 9;
    [Sx, Sy] = meshgrid(linspace(xrange(1), xrange(2), nseed), linspace(yrange(1), yrange(2), nseed));
    seeds = repmat(eq0(:).', numel(Sx), 1);
    seeds(:, idx(1)) = Sx(:);
    seeds(:, idx(2)) = Sy(:);

    opts = odeset('RelTol', 1e-6, ...
                  'Events', @(t, y) leave_box(t, y, box));
    % ode15s (stiff) because the forward horizon spans many fast timescales: a
    % non-stiff solver would take a fast-timescale step across the whole slow drift
    for k = 1:size(seeds, 1)
        [~, Yf] = ode15s(odefun, [0 Tf], seeds(k, :)', opts);
        [~, Yb] = ode15s(odefun, [0 -Tb], seeds(k, :)', opts);
        plot_trajectory(Yf, idx, box);
        plot_trajectory(Yb, idx, box);
    end

    % direction field: normalized arrows on a coarse grid show flow direction
    [Ax, Ay, U, V] = plane_field(fi, fj, eq0, idx, xrange, yrange, 10);
    quiver(Ax, Ay, U, V, 0.3, 'Color', [0.6 0.6 0.6], 'MaxHeadSize', 0.2);

    draw_equilibria(eqpts, cls, idx, 5);
    finalize_axes(labels, xrange, yrange, 'Phase portrait (trajectories)');
end

% Plot one trajectory, dropping any samples that ran outside the (padded) box.
% A stiff solver failing near a saddle can emit a huge final point; drawing it
% would streak a spurious line across the plot.
function plot_trajectory(Y, idx, box)
    pad = 0.1*(box(:, 2) - box(:, 1));
    lo = (box(:, 1) - pad).';
    hi = (box(:, 2) + pad).';
    inside = all(Y >= lo & Y <= hi, 2);
    Y(~inside, :) = NaN; % NaN breaks the line instead of streaking to the runaway point
    plot(Y(:, idx(1)), Y(:, idx(2)), 'Color', [0.3 0.4 0.8]);
end

function [X, Y, U, V, speed] = plane_field(fi, fj, eq0, idx, xrange, yrange, n)
    % grid of unit-length flow vectors (NaN at equilibria) plus raw speed, over
    % the plotted plane with the off-plane states held at eq0
    [X, Y] = meshgrid(linspace(xrange(1), xrange(2), n), ...
                      linspace(yrange(1), yrange(2), n));
    args = field_args(eq0, idx, X, Y);
    U = fi(args{:}).*ones(size(X));
    V = fj(args{:}).*ones(size(X));
    speed = hypot(U, V);
    U = U ./ speed;
    V = V ./ speed;
end

% Argument list for the n-state component handles: the two plotted states vary
% over the grid (X, Y), every off-plane state is held at its eq0 value.
function args = field_args(eq0, idx, X, Y)
    args = num2cell(eq0(:).');
    args{idx(1)} = X;
    args{idx(2)} = Y;
end

function [eqpts, cls] = equilibria(F_sym, u_sym, box)
    % equilibria: solve F = 0, keep real in-box points, classify by Jacobian
    n = numel(u_sym);
    Jfun = matlabFunction(jacobian(F_sym, u_sym), 'Vars', {u_sym});
    sol = solve(F_sym == 0, u_sym, 'Real', true);
    cols = arrayfun(@(k) double(sol.(char(u_sym(k)))), 1:n, 'UniformOutput', false);
    P = [cols{:}]; % one candidate equilibrium per row

    in = true(size(P, 1), 1);
    for k = 1:n
        in = in & imag(P(:, k)) == 0 & real(P(:, k)) >= box(k, 1) & real(P(:, k)) <= box(k, 2);
    end
    eqpts = real(P(in, :));

    cls = zeros(size(eqpts, 1), 1); % 1 stable, 2 unstable, 3 saddle
    for r = 1:size(eqpts, 1)
        re = real(eig(Jfun(eqpts(r, :).')));
        if all(re < 0)
            cls(r) = 1;
        elseif all(re > 0)
            cls(r) = 2;
        else
            cls(r) = 3;
        end
    end
end

% Reference point for the field slice: the in-window equilibrium nearest the
% centre of the plotted box, or the box centre if no equilibrium lies in it.
function eq0 = slice_reference(eqpts, box, idx)
    c = mean(box, 2);
    if isempty(eqpts)
        eq0 = c;
        return
    end
    d = sum((eqpts(:, idx) - c(idx).').^2, 2);
    [~, r] = min(d);
    eq0 = eqpts(r, :).';
end

function draw_equilibria(eqpts, cls, idx, markersize)
    % project the classified equilibria onto the plotted plane
    % filled circle = stable, open circle = unstable, cross = saddle
    specs = {1, 'o', 'k', 'stable'; 2, 'o', 'none', 'unstable'; 3, 'x', 'none', 'saddle'};
    hLeg = gobjects(0); labels = {};
    for c = 1:size(specs, 1)
        m = cls == specs{c, 1};
        if any(m)
            hLeg(end+1) = plot(eqpts(m, idx(1)), eqpts(m, idx(2)), specs{c, 2}, ...
                'MarkerFaceColor', specs{c, 3}, 'MarkerEdgeColor', 'k', ...
                'Color', 'k', 'MarkerSize', markersize, 'LineWidth', 1);
            labels{end+1} = specs{c, 4};
        end
    end
    if ~isempty(hLeg)
        legend(hLeg, labels, 'Location', 'best');
    end
end

function finalize_axes(labels, xrange, yrange, titleStr)
    xlabel(labels{1});
    ylabel(labels{2});
    title(titleStr);
    xlim(xrange);
    ylim(yrange);
    daspect([1 1 1]); % equal scale on both axes
end

function [val, term, dir] = leave_box(~, y, box)
% stop integrating once the trajectory leaves the (slightly padded) box
    pad = 0.05*(box(:, 2) - box(:, 1));
    inside = all(y >= box(:, 1) - pad) && all(y <= box(:, 2) + pad);
    val = double(inside) - 0.5; % crosses zero when the trajectory exits
    term = 1;
    dir = -1;
end
