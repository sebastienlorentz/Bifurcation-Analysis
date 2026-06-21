% Find and trace every fold and Hopf curve of a two-parameter system in a box.
% Seeds come from 1-parameter sweeps at several fixed p2 slices; each fold/Hopf
% detected seeds a 2-parameter curve, which is then continued and de-duplicated.
function out = trace_bifurcation_curves(F_sym, u_sym, p_sym, p_box, u_start, state_box, n_slices, inv_idx)
    u_sym = u_sym(:);
    p_sym = p_sym(:);
    n = numel(u_sym);
    if nargin < 6 || isempty(state_box)
        state_box = [-inf(n, 1), inf(n, 1)];
    end
    if nargin < 7 || isempty(n_slices)
        n_slices = 15;
    end
    if nargin < 8
        inv_idx = [];
    end

    seed_step = 1e-3; % step for the 1-parameter seeding sweeps
    seed_max = 1e5;
    curve_step = 1e-3; % step for tracing each 2-parameter curve
    curve_max = 1e5;

    % build the two-parameter fold and Hopf problems once
    fsys = fold_continuation_system(F_sym, u_sym, p_sym);
    hsys = hopf_continuation_system(F_sym, u_sym, p_sym);

    % seed by sweeping 1-parameter continuations at several p2 slices. The sweep
    % runs over a p1 range extended below the box so that a branch point sitting
    % on the box edge (e.g. a transcritical at R0=1) is still crossed and the
    % branches beyond it get traced; seeds are then kept only inside the box.
    seed_lo = p_box(1, 1) - diff(p_box(1, :));
    seed_hi = p_box(1, 2);
    fold_seeds = zeros(n + 2, 0);
    hopf_seeds = zeros(n + 2, 0);
    bp_seeds = zeros(n + 2, 0);
    for pv = linspace(p_box(2, 1), p_box(2, 2), n_slices)
        sys = continuation_system(subs(F_sym, p_sym(2), pv), u_sym, p_sym(1));
        out = trace_branches(sys, u_start, seed_hi, [seed_lo seed_hi], seed_step, seed_max, true, state_box);
        fold_seeds = [fold_seeds, seeds_from(out.folds, n, pv, p_box(1, :))]; 
        hopf_seeds = [hopf_seeds, seeds_from(out.hopfs, n, pv, p_box(1, :))]; 
        bp_seeds = [bp_seeds, seeds_from({out.branch_points}, n, pv, p_box(1, :))]; 
    end

    % continue each seed into a curve, skipping seeds already on a traced curve
    out.fold_curves = continue_seeds(fold_seeds, n, p_box, @(u0, p0) trace_fold_curve(fsys, u0, p0, p_box, curve_step, curve_max));
    out.hopf_curves = continue_seeds(hopf_seeds, n, p_box, @(u0, p0) trace_hopf_curve(hsys, u0, p0, p_box, curve_step, curve_max));

    % transcritical curves need the invariant-subspace structure, so they are
    % only traced when the invariant ("infected") state index is supplied
    out.transcritical_curves = {};
    if ~isempty(inv_idx)
        tsys = transcritical_continuation_system(F_sym, u_sym, p_sym, inv_idx);
        out.transcritical_curves = continue_seeds(bp_seeds, n, p_box, @(u0, p0) trace_transcritical_curve(tsys, u0, p0, p_box, curve_step, curve_max));
    end

    % first Lyapunov coefficient along each Hopf curve: its sign is the Hopf
    % criticality (super/sub), and a sign change is a generalized Hopf (Bautin)
    lsys = lyapunov_system(F_sym, u_sym, p_sym);
    for i = 1:numel(out.hopf_curves)
        c = out.hopf_curves{i};
        l1 = zeros(1, size(c.p, 2));
        for k = 1:size(c.p, 2)
            l1(k) = lyapunov_coefficient(lsys, c.u(:, k), c.p(:, k));
        end
        out.hopf_curves{i}.l1 = l1;
    end
    out.gh_points = find_gh_points(out.hopf_curves);

    % Bogdanov-Takens points: codim-2 points where a Hopf curve meets a fold
    % curve, found as the zeros of the Hopf frequency along the Hopf curves
    out.bt_points = find_bt_points(out.hopf_curves);
end

% Generalized Hopf points = where the first Lyapunov coefficient changes
% sign along a Hopf curve (the Hopf switches super/subcritical). Detection uses
% the finite quantity l1*|omega|
function gh = find_gh_points(hopf_curves)
    gh = zeros(2, 0);
    for i = 1:numel(hopf_curves)
        test = hopf_curves{i}.l1.*abs(hopf_curves{i}.omega);
        P = hopf_curves{i}.p;
        for k = find(test(1:end-1).*test(2:end) < 0)
            s = test(k) / (test(k) - test(k+1));
            gh(:, end+1) = (1 - s)*P(:, k) + s*P(:, k+1); 
        end
    end
end

% Bogdanov-Takens points = where the Hopf frequency omega passes through zero
% Located by interpolating the parameter point at the sign change of omega.
function bt = find_bt_points(hopf_curves)
    bt = zeros(2, 0);
    for i = 1:numel(hopf_curves)
        w = hopf_curves{i}.omega;
        P = hopf_curves{i}.p;
        for k = find(w(1:end-1).*w(2:end) < 0)
            s = w(k) / (w(k) - w(k+1));
            bt(:, end+1) = (1 - s)*P(:, k) + s*P(:, k+1); 
        end
    end
end

% Collect [state; p1; p2] seed columns from a cell of detected bifurcation points,
% keeping only those whose p1 lies inside the box.
function seeds = seeds_from(cells, n, pv, p1_range)
    seeds = zeros(n + 2, 0);
    for i = 1:numel(cells)
        B = cells{i};
        for j = 1:size(B, 2)
            p1 = B(end, j);
            if p1 >= p1_range(1) && p1 <= p1_range(2)
                seeds(:, end+1) = [B(1:n, j); p1; pv]; 
            end
        end
    end
end

% Continue every seed, skipping any that already lies on a traced curve.
function curves = continue_seeds(seeds, n, p_box, tracer)
    curves = {};
    for s = seeds
        u0 = s(1:n);
        p0 = s(n+1:n+2);
        if on_existing_curve(p0, curves, p_box)
            continue
        end
        curves{end+1} = tracer(u0, p0); 
    end
end

% True if parameter point p0 already lies on one of the traced curves (distance
% measured relative to the box size so both parameters count equally).
function tf = on_existing_curve(p0, curves, p_box)
    tf = false;
    scale = [diff(p_box(1, :)); diff(p_box(2, :))];
    for i = 1:numel(curves)
        if min(vecnorm((curves{i}.p - p0) ./ scale)) < 0.02
            tf = true;
            return
        end
    end
end
