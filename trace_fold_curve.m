% Trace a fold curve through the (p1,p2) parameter plane, starting
% from a known fold. This is the two-parameter analogue of trace_branches and
% reuses the same generic pseudo_arclength engine.
function out = trace_fold_curve(fsys, u0, p0, p_box, step_size, max_steps)
    n = fsys.n;
    u0 = u0(:);
    p0 = p0(:);

    % null vector of Fu at the start: used as the v-guess and as the fixed
    % normalization vector c in c'v = 1
    [~, ~, V] = svd(fsys.Fu(u0, p0));
    v0 = V(:, end);
    prob = fsys.make_problem(v0);

    % land exactly on the fold curve: hold p2 fixed and solve the fold equations
    % for (u, v, p1) so the residual vanishes
    land = @(z) [prob.residual(z); z(end) - p0(2)];
    z0 = fsolve(land, [u0; v0; p0], optimset('Display', 'off'));

    % tangent, box stop and both arcs are the shared two-parameter curve tail
    pts = follow_curve(prob, z0, p_box, step_size, max_steps);

    out.u = pts(1:n, :);
    out.p = pts(2*n+1:2*n+2, :);
end
