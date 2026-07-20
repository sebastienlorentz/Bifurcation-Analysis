% Trace a fold-of-cycles (LPC) curve through the (p1,p2) parameter plane from a
% seed periodic orbit. Two-parameter analogue of trace_fold_curve for cycles: it
% lands on the LPC curve and follows it with pseudo_arclength on the
% lpc_continuation_system. The seed is an orbit at a fold of cycles (a nontrivial
% Floquet multiplier = +1), e.g. the turning point of a one-parameter limit-cycle
% branch from trace_limit_cycle. x0/T0 are that orbit's section point and period,
% p0 the parameter pair there.
function out = trace_lpc_curve(lcsys, x0, T0, p0, p_box, step_size, max_steps)
    n = lcsys.n;
    x0 = x0(:);
    p0 = p0(:);
    if nargin < 6 || isempty(step_size); step_size = 1e-3; end
    if nargin < 7 || isempty(max_steps); max_steps = 1e4;  end

    % Poincare section fixed at the seed cycle's I-level, held along the curve
    prob = lcsys.make_problem(x0(2));

    % land on the LPC curve: hold p2 fixed and solve the LPC equations for
    % (x0, T, p1) so periodicity, section and det(M)-1 all vanish. Pass the
    % analytic Jacobian (bordered with the fixed-p2 row) and tight tolerances so
    % the seed lands well below the corrector's tolerance, a numerical Jacobian
    % over the ODE-integrated residual stalls near 1e-4 and stops the arc dead.
    land = @(z) deal([prob.residual(z); z(end) - p0(2)], ...
                     [prob.jacobian(z); zeros(1, n+2), 1]);
    opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true, ...
                        'FunctionTolerance', 1e-12, 'StepTolerance', 1e-14, 'MaxIterations', 200);
    z0 = fsolve(land, [x0; T0; p0], opts);

    % tangent, box stop and both arcs are the shared two-parameter curve tail
    [pts, fwd_reason, bwd_reason, fwd_n, bwd_n] = follow_curve(prob, z0, p_box, step_size, max_steps);

    out.x = pts(1:n, :);       % section point of the cycle along the curve
    out.T = pts(n+1, :);       % period along the curve
    out.p = pts(n+2:n+3, :);   % the LPC curve in the (p1, p2) plane
    % why each arm ended, and its length, so a caller can check the curve hit
    % its true terminus ('nonfinite'/'stop') rather than running out of steps
    out.fwd_reason = fwd_reason;   out.fwd_steps = fwd_n;
    out.bwd_reason = bwd_reason;   out.bwd_steps = bwd_n;
end
