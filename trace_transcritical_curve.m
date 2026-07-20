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

    % tangent, box stop and both arcs are the shared two-parameter curve tail
    pts = follow_curve(prob, z0, p_box, step_size, max_steps);

    out.u = pts(1:n, :);
    out.p = pts(n+1:n+2, :);
end
