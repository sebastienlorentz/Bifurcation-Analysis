% Trace a transcritical (branch-point) curve through the (p1,p2) parameter plane,
% starting from a detected branch point. Two-parameter analogue of trace_branches
function out = trace_transcritical_curve(tsys, u0, p0, p_box, step_size, max_steps)
    n = tsys.n;
    u0 = u0(:);
    p0 = p0(:);
    prob = tsys.make_problem();

    % land exactly on the curve: hold p2 fixed and solve for (u, p1)
    opts = optimset('Display', 'off');
    land = @(z) [prob.residual(z); z(end) - p0(2)];
    z0 = fsolve(land, [u0; p0], opts);

    % initial tangent = null direction of the augmented Jacobian
    t0 = null(prob.jacobian(z0));
    t0 = t0(:, 1) / norm(t0(:, 1));

    % stop once either parameter leaves its box
    stop = @(z) z(end-1) < p_box(1, 1) || z(end-1) > p_box(1, 2) || z(end) < p_box(2, 1) || z(end) > p_box(2, 2);

    % follow the curve both ways from the start point
    fwd = pseudo_arclength(prob, z0, t0, step_size, max_steps, stop);
    bwd = pseudo_arclength(prob, z0, -t0, step_size, max_steps, stop);
    pts = [fliplr(bwd(:, 2:end)), fwd]; % drop the shared start point in bwd

    out.u = pts(1:n, :);
    out.p = pts(n+1:n+2, :);
end
