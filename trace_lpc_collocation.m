% Trace a fold-of-cycles (LPC) curve in the (p1,p2) plane by orthogonal
% collocation, the robust, fast replacement for trace_lpc_curve (single shooting).
% lsys comes from lpc_collocation_system. The seed (x0, T0, p0) is a fold-of-cycles
% orbit from find_lpc_seed (section point, period, parameter pair). Lands on the LPC
% curve, then follows it both ways with pseudo_arclength on the collocation extended
% system. Returns the curve in (p1,p2), the period along it, and the per-arm reasons.
function out = trace_lpc_collocation(lsys, x0, T0, p0, p_box, step_size, max_steps)
    if nargin < 6 || isempty(step_size); step_size = 1e-2; end
    if nargin < 7 || isempty(max_steps); max_steps = 200; end
    p0 = p0(:);
    n = lsys.n; Np = lsys.Np;

    % reference profile (fixes the integral phase condition) from the seed orbit
    Xref = lsys.seed_profile(x0(:), T0, p0);
    prob = lsys.make_problem(Xref);
    nW = prob.nW;

    % land the cycle exactly at p0
    W = [Xref(:); T0];
    for it = 1:20
        Rc = prob.orbit.residual([W; p0]);
        if norm(Rc) < 1e-11; break; end
        Jo = prob.orbit.jacobian([W; p0]);
        W = W - Jo(:, 1:nW)\Rc;
    end

    % null vector of the cycle Jacobian: the fold direction, seed for V
    Jo = prob.orbit.jacobian([W; p0]);
    [~, ~, Vv] = svd(Jo(:, 1:nW));
    V = Vv(:, end);
    Z0 = [W; V; p0];

    % land on the LPC curve holding p2 fixed (analytic Jacobian, bordered row)
    prow = [zeros(1, 2*nW), 0, 1];
    land = @(Z) deal([prob.residual(Z); Z(end) - p0(2)], [prob.jacobian(Z); prow]);
    opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true, ...
                        'FunctionTolerance', 1e-11, 'StepTolerance', 1e-12, 'MaxIterations', 100);
    Z0 = fsolve(land, Z0, opts);

    % terminate cleanly where the cycle collapses (a Bautin point: the fold of
    % cycles is born there at zero amplitude, so the extended system degenerates).
    % Amplitude is the peak-to-peak of the invariant state (index 2) over the profile.
    nX = n*Np;
    amp0 = ptp(Z0, nX, n);
    amp_floor = 0.02*amp0;
    extra_stop = @(Z) ptp(Z, nX, n) < amp_floor;

    % follow the curve both ways (params are the last two entries of Z)
    [pts, fwd_reason, bwd_reason, fwd_n, bwd_n] = follow_curve(prob, Z0, p_box, step_size, max_steps, extra_stop);

    out.p = pts(2*nW+1:2*nW+2, :);   % the LPC curve in (p1, p2)
    out.T = pts(nW, :);              % period along the curve
    out.x = pts(1:n, :);            % first node of the cycle along the curve
    out.fwd_reason = fwd_reason; out.fwd_steps = fwd_n;
    out.bwd_reason = bwd_reason; out.bwd_steps = bwd_n;
end

% Peak-to-peak amplitude of the invariant state (index 2) over the cycle profile.
function a = ptp(Z, nX, n)
    I = Z(2:n:nX);
    a = max(I) - min(I);
end

