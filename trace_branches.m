% Trace every equilibrium branch of F(u,lambda)=0, switching at branch points.
function out = trace_branches(sys, x0, lambda0, lambda_range, step_size, max_steps, decreasing, state_range)
    % state_range (optional) is an n-by-2 [min max] box; a branch stops on leaving
    % it. Default unbounded. Cuts off nonphysical branches that would hit max_steps.
    if nargin < 8 || isempty(state_range)
        state_range = [-inf(numel(x0), 1), inf(numel(x0), 1)];
    end

    dup_tol = 1e-4; % points closer than this count as the same one
    dir_tol = 1e-3; % cosine > 1-dir_tol means same direction
    win_arc = 10*step_size; % arclength a sign change must persist over to count

    % the continuation problem in the augmented vector z = [state; lambda]:
    % residual is F and its Jacobian is [Fu, Flambda]. A branch stops once
    % lambda leaves its range or the state leaves the allowed box.
    prob.residual = @(z) sys.F(z(1:end-1), z(end));
    prob.jacobian = @(z) [sys.Fu(z(1:end-1), z(end)), sys.Flambda(z(1:end-1), z(end))];
    % box_tol lets a branch sit exactly on a boundary (e.g. the disease-free
    % I = 0 branch) without a tiny fsolve residual pushing it out of the box
    box_tol = 1e-6;
    stop = @(z) z(end) < lambda_range(1) || z(end) > lambda_range(2) || any(z(1:end-1) < state_range(:, 1) - box_tol) || any(z(1:end-1) > state_range(:, 2) + box_tol);

    % land on the curve before we start
    x0 = fsolve(@(u) sys.F(u, lambda0), x0(:), optimset('Display', 'off'));
    t0 = initial_tangent(x0, lambda0, sys.Fu, sys.Flambda, decreasing);

    out.branches = {};
    out.folds = {};
    out.hopfs = {};
    out.stable = {};
    all_bps = []; % branch points we've already found
    done_seeds = []; % [point; tangent] seeds we've already launched

    y0 = [x0; lambda0];
    queue = {{y0, t0}, {y0, -t0}}; % go both ways along the arc

    while ~isempty(queue)
        seed = queue{1};
        queue(1) = [];
        y0 = seed{1};
        tau = seed{2};

        % don't re-follow a seed we've already done. A seed is a repeat if it
        % sits at the same point AND points the same way as one we've launched.
        m = numel(y0);
        if ~isempty(done_seeds)
            seed_pts = done_seeds(1:m, :);
            seed_dirs = done_seeds(m+1:end, :);
            same_point = vecnorm(seed_pts-y0) < dup_tol;
            same_dir = tau'*seed_dirs > 1-dir_tol;
            if any(same_point & same_dir)
                continue
            end
        end
        done_seeds(:, end+1) = [y0; tau];

        [pts, tans, arcl] = pseudo_arclength(prob, y0, tau, step_size, max_steps, stop);
        dfu = zeros(1, size(pts, 2)); % det(Fu) along the branch, for detection
        for k = 1:size(pts, 2)
            dfu(k) = det(sys.Fu(pts(1:end-1, k), pts(end, k)));
        end

        [bp_pts, bp_idx] = find_branch_switches(sys, pts, tans, dfu, arcl, win_arc);

        if ~isempty(bp_idx)
            [k1, w] = min(bp_idx);
            bp = bp_pts(:, w);
            tau_arr = tans(:, k1); % how we came into the bp

            % only keep the arc up to that first branch point
            pts = pts(:, 1:k1);
            tans = tans(:, 1:k1);
            arcl = arcl(1:k1);

            % handle each branch point once. Shoot off down the other curves
            % through it (both signs), but skip the way we came in.
            if isempty(all_bps) || all(vecnorm(all_bps - bp) >= dup_tol)
                all_bps(:, end+1) = bp;
                T = branch_point_tangents(sys, bp(1:end-1), bp(end));
                for c = 1:size(T, 2)
                    for s = [1, -1]
                        d = s*T(:, c);
                        if d'*tau_arr < -(1 - dir_tol) % straight back the way we came
                            continue
                        end
                        queue{end+1} = {bp, d};
                    end
                end
            end
        end

        [folds, hopfs] = detect_bifurcations(pts, tans, sys.Fu, arcl, win_arc);

        % stable = all Jacobian eigenvalues have negative real part; the tolerance
        % stops points on a bifurcation flickering stable/unstable
        n = size(pts, 1) - 1;
        stable = false(1, size(pts, 2));
        for k = 1:size(pts, 2)
            ev = eig(sys.Fu(pts(1:n, k), pts(end, k)));
            stable(k) = all(real(ev) < 1e-9);
        end

        out.branches{end+1} = pts;
        out.folds{end+1} = folds;
        out.hopfs{end+1} = hopfs;
        out.stable{end+1} = stable;
    end
    out.branch_points = all_bps;

    % first Lyapunov coefficient at each detected Hopf (sign = criticality;
    % sys carries Fu/Fuu/Fuuu so it doubles as a lyapunov_system)
    out.hopf_l1 = cell(size(out.hopfs));
    for i = 1:numel(out.hopfs)
        H = out.hopfs{i};
        l1 = zeros(1, size(H, 2));
        for k = 1:size(H, 2)
            l1(k) = lyapunov_coefficient(sys, H(1:end-1, k), H(end, k));
        end
        out.hopf_l1{i} = l1;
    end
end

% First tangent to the curve at the starting point.
function t = initial_tangent(x0, lambda0, Fx, Flambda, decreasing)
    A = Fx(x0, lambda0);
    b = -Flambda(x0, lambda0);
    xdot = A \ b;
    t = [xdot; 1];
    t = t / norm(t);
    if decreasing
        t = -t;
    end
end

% Locate folds and Hopf points along an already-traced branch. Works in any
% dimension: the Hopf test function is det(2J (.) I), the determinant of the
% bialternate product, whose eigenvalues are the pairwise sums lambda_i+lambda_j.
% It vanishes when a pair reaches the imaginary axis (Hopf) or a real pair sums
% to zero (neutral saddle); the two are told apart from the eigenvalues at the
% crossing. For n = 2 the bialternate product is the 1-by-1 matrix [trace(J)],
% so this reduces to the old trace(J) = 0 test.
function [folds, hopfs] = detect_bifurcations(points, tangents, Fu, svec, win_arc)
    N = size(points, 2);
    n = size(points, 1) - 1;
    t_lam = tangents(end, :); % lambda-component of the tangent

    psi = zeros(1, N); % det(2J (.) I), the general Hopf test function
    for k = 1:N
        psi(k) = det(bialt(Fu(points(1:n, k), points(n+1, k))));
    end

    % fold = lambda-component of the tangent flips sign
    fold_br = find(t_lam(1:N-1).*t_lam(2:N) < 0);
    fold_br = filter_sustained(t_lam, fold_br, svec, win_arc);
    folds = interp_bif(points, t_lam, fold_br);

    % Hopf candidates: the test function flips sign
    hopf_br = find(psi(1:N-1).*psi(2:N) < 0);
    hopf_br = filter_sustained(psi, hopf_br, svec, win_arc);

    % keep crossings where the pair reaching the axis is complex (a Hopf); drop
    % real pairs summing to zero (neutral saddles)
    keep = false(1, numel(hopf_br));
    for i = 1:numel(hopf_br)
        k = hopf_br(i);
        keep(i) = hopf_pair_is_complex(Fu(points(1:n, k), points(n+1, k)));
    end
    hopf_br = hopf_br(keep);
    hopfs = interp_bif(points, psi, hopf_br);
end

% True if the eigenvalue pair nearest to summing to zero is complex, i.e. the
% sign change in det(2J (.) I) is a Hopf rather than a neutral saddle.
function tf = hopf_pair_is_complex(J)
    ev = eig(J);
    n = numel(ev);
    S = ev + ev.'; % S(i,j) = lambda_i + lambda_j
    S(1:n+1:end) = inf; % ignore the diagonal (i = j)
    [~, w] = min(abs(S(:)));
    [i, ~] = ind2sub([n, n], w);
    tf = abs(imag(ev(i))) > 1e-6;
end

% Bialternate product 2A (.) I: an m-by-m matrix, m = n(n-1)/2, indexed by pairs
% (p,q) with p > q. Its eigenvalues are the sums lambda_i + lambda_j over i > j,
% so it is singular exactly at a Hopf or a neutral saddle (Kuznetsov, Elements
% of Applied Bifurcation Theory).
function B = bialt(A)
    n = size(A, 1);
    idx = nchoosek(1:n, 2); % rows [q p] with q < p
    m = size(idx, 1);
    B = zeros(m, m);
    for a = 1:m
        q = idx(a, 1);
        p = idx(a, 2);
        for b = 1:m
            s = idx(b, 1);
            r = idx(b, 2);
            if r == q
                v = -A(p, s);
            elseif r ~= p && s == q
                v = A(p, r);
            elseif r == p && s == q
                v = A(p, p) + A(q, q);
            elseif r == p && s ~= q
                v = A(q, s);
            elseif s == p
                v = -A(q, r);
            else
                v = 0;
            end
            B(a, b) = v;
        end
    end
end

% Drop sign-change brackets that don't persist, i.e. numerical noise.
function kept = filter_sustained(test, brackets, svec, win_arc)
% Keep only brackets whose sign persists for at least win_arc of arclength on both
% sides, rejecting single-point noise flips. svec is the cumulative arclength.
    kept = [];
    for i = 1:numel(brackets)
        k = brackets(i);
        lo = find(svec <= svec(k) - win_arc, 1, 'last');
        hi = find(svec >= svec(k+1) + win_arc, 1, 'first');
        if isempty(lo) || isempty(hi)
            continue
        end
        if sign(test(lo)) == sign(test(k)) && sign(test(hi)) == sign(test(k + 1))
            kept(end+1) = k;
        end
    end
end

% Linearly interpolate a bifurcation point within each sign-change bracket.
function bp = interp_bif(points, tau, brackets)
    bp = zeros(size(points, 1), numel(brackets));
    for i = 1:numel(brackets)
        k = brackets(i);
        s = tau(k) / (tau(k)-tau(k+1));
        bp(:, i) = (1 - s)*points(:, k) + s*points(:, k+1);
    end
end

% Find branch points (det(Fu)=0 without a fold) along a branch.
function [bp_pts, bp_idx] = find_branch_switches(sys, points, tangents, detFu, svec, win_arc)
    n = size(points, 1)-1;
    N = numel(detFu);
    t_lam = tangents(end, :);

    % branch point: det(Fu) flips sign but t_lambda doesn't. At a fold t_lambda
    % flips too, so this condition leaves folds out.
    det_cross = find(detFu(1:N-1).*detFu(2:N) < 0);
    no_fold = t_lam(det_cross).*t_lam(det_cross+1) > 0;
    brackets = det_cross(no_fold);

    brackets = filter_sustained(detFu, brackets, svec, win_arc);

    nb = numel(brackets);
    bp_pts = zeros(n+1, nb);
    for i = 1:nb
        k = brackets(i);
        [x_star, lam_star] = localise_bp(points(:, k), points(:, k+1), detFu(k), detFu(k+1), sys.F, sys.Fu, 1e-10, 30);
        bp_pts(:, i) = [x_star; lam_star];
    end
    bp_idx = brackets;
end

% Pin down a branch point between two stored points.
function [x_star, lam_star] = localise_bp(p1, p2, d1, d2, F, Fx, tol, max_iter)
% Find det(Fx) = 0 between p1 and p2 by regula falsi, Newton-correcting back onto
% the solution curve before testing det at each guess.
    n = length(p1)-1;
    s_lo = 0;
    s_hi = 1;
    dlo = d1;
    dhi = d2;

    x = p1(1:n);
    lam = p1(end);
    for iter = 1:max_iter
        % linearly interpolate s to where det should hit zero
        s = s_lo + dlo*(s_hi-s_lo) / (dlo-dhi);
        y = (1-s)*p1 + s*p2;
        x = y(1:n);
        lam = y(end);

        % Newton step back onto F = 0 (lambda held fixed)
        for k = 1:20
            r = F(x, lam);
            if norm(r) < tol
                break
            end
            x = x - Fx(x, lam) \ r;
        end

        d = det(Fx(x, lam));
        if abs(d) < tol
            break
        end

        % keep the half of the bracket that still straddles the sign change
        if sign(d) == sign(dlo)
            s_lo = s;
            dlo = d;
        else
            s_hi = s;
            dhi = d;
        end
    end

    x_star = x;
    lam_star = lam;
end

% Tangents of the two curves crossing a simple branch point.
function T = branch_point_tangents(sys, x_star, lam_star)
% From the branching equation a*alpha^2 + 2b*alpha + c = 0:
%   transcritical -> two finite roots
%   pitchfork -> a = 0, one smooth root plus the null vector phi (vertical in lambda)
    n = length(x_star);

    Fu = sys.Fu(x_star, lam_star);
    Fl = sys.Flambda(x_star, lam_star);

    [U, ~, V] = svd(Fu);
    phi = V(:, end); % right null vector
    psi = U(:, end); % left null vector

    % particular solution v from the bordered system
    B = [Fu, phi; phi', 0];
    sol = B \ [-Fl; 0];
    v = sol(1:n);

    Huu = sys.Fuu(x_star, lam_star);
    Hul = sys.Fulambda(x_star, lam_star);
    Hll = sys.Flambdalambda(x_star, lam_star);

    a = psi'*bilin(Huu, phi, phi);
    b = psi'*(Hul*phi + bilin(Huu, phi, v));
    c = psi'*(Hll + Hul*v + bilin(Huu, v, v));

    if abs(a) <= 1e-6*max([abs(b), abs(c), 1])
        % pitchfork: smooth branch, plus the bifurcating one along phi
        alpha = -c / (2*b);
        tau1 = [alpha*phi + v; 1];
        tau2 = [phi; 0];
    else
        disc = max(b^2 - a*c, 0);
        alpha1 = (-b + sqrt(disc)) / a;
        alpha2 = (-b - sqrt(disc)) / a;
        tau1 = [alpha1*phi + v; 1];
        tau2 = [alpha2*phi + v; 1];
    end

    T = [tau1 / norm(tau1), tau2 / norm(tau2)];
end
