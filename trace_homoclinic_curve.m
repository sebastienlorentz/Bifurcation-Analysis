% Trace a homoclinic curve in the (p1,p2) plane by orthogonal collocation. hsys comes
% from homoclinic_collocation_system. The seed is a large-period cycle profile (X, T, p)
% taken off the LPC arc: it is cut at the saddle passage, landed onto the homoclinic
% BVP with p1 pinned, then followed both ways with pseudo_arclength. Returns the curve
% in (p1,p2), the truncation interval and saddle along it, and the per-arm reasons.
% Seed with a moderate period (T ~ 60 for the SIS model); the largest-T orbit on the
% LPC arc is under-resolved on a uniform mesh and will not land.
function out = trace_homoclinic_curve(hsys, X, T, p, p_box, step_size, max_steps, delta)
    if nargin < 6 || isempty(step_size); step_size = 1e-2; end
    if nargin < 7 || isempty(max_steps); max_steps = 400; end
    if nargin < 8; delta = []; end
    ix = hsys.idx; n = hsys.n; Np = hsys.Np;

    [Z0, delta, kfar] = hsys.seed_from_cycle(X, T, p(:), delta);
    prob = hsys.make_problem(reshape(Z0(ix.X), n, Np), delta, kfar);

    % land on the curve with p1 held: a homoclinic is codim-1, so exactly one
    % parameter stays free and the system is square.
    Z0 = land(prob, Z0, ix, T);

    % Bracket T from both sides. Below: the truncated BVP admits a degenerate T -> 0
    % solution, the orbit shrinking onto the saddle. Above: the curve ends at a
    % Bogdanov-Takens point, where the saddle's eigenvalues vanish and the transit
    % time diverges, so the parameters converge on BT while T marches off to infinity.
    % Stopping on lu -> 0 does not work: it is approached asymptotically and never
    % reaches any floor, leaving the continuation to crawl at fixed parameters.
    T_lo = 0.25*Z0(ix.T);
    T_hi = 10*Z0(ix.T);
    extra_stop = @(Z) Z(ix.T) < T_lo || Z(ix.T) > T_hi;

    [pts, fwd_reason, bwd_reason, fwd_n, bwd_n] = follow_curve(prob, Z0, p_box, step_size, max_steps, extra_stop);

    out.p = pts(ix.p, :);            % the homoclinic curve in (p1, p2)
    out.T = pts(ix.T, :);            % truncation interval along the curve
    out.q = pts(ix.q, :);            % the saddle along the curve
    out.lu = pts(ix.lu, :);          % the unstable eigenvalue, -> 0 at a BT point
    out.X = pts(ix.X, :);            % full orbit profile along the curve
    out.eps = pts(ix.e, :);          % truncation distance, shrinking with the loop
    out.delta = delta;

    % leading stable eigenvalue, for the saddle quantity lu + ls whose sign decides
    % whether the bifurcating cycle is stable. Not part of the system (the stable
    % subspace is only ever used implicitly), so it is recovered here.
    out.ls = zeros(1, size(pts, 2));
    for k = 1:size(pts, 2)
        ev = eig(hsys.Fu(pts(ix.q, k), pts(ix.p, k)));
        out.ls(k) = max(real(ev(real(ev) < 0)));
    end
    out.fwd_reason = fwd_reason; out.fwd_steps = fwd_n;
    out.bwd_reason = bwd_reason; out.bwd_steps = bwd_n;
end

% Damped Newton onto the homoclinic BVP with p1 pinned, keeping T away from zero.
function Z = land(prob, Z0, ix, T0)
    Z = Z0;
    pin = zeros(1, ix.N); pin(ix.p(1)) = 1;
    for it = 1:60
        R = [prob.residual(Z); Z(ix.p(1)) - Z0(ix.p(1))];
        if norm(R) < 1e-9; return; end
        dZ = [prob.jacobian(Z); pin]\R;
        s = 1;
        while Z(ix.T) - s*dZ(ix.T) < 0.25*T0 && s > 1e-3
            s = s/2;
        end
        Z = Z - s*dZ;
    end
    error('trace_homoclinic_curve:land', 'could not land on the homoclinic (|R| = %.2e)', norm(R));
end
