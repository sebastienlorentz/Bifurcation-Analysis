% Trace a Hopf curve through the (p1,p2) parameter plane, starting from a known
% Hopf point. Two-parameter analogue of trace_branches;
function out = trace_hopf_curve(hsys, u0, p0, p_box, step_size, max_steps)
    n = hsys.n;
    u0 = u0(:);
    p0 = p0(:);

    % the eigenvalue pair nearest the imaginary axis (with positive frequency)
    % is the Hopf pair; take its frequency and eigenvector as the start guess
    [V, D] = eig(hsys.Fu(u0, p0));
    ev = diag(D);
    cand = find(imag(ev) > 1e-6);
    [~, k] = min(abs(real(ev(cand))));
    k = cand(k);
    om0 = imag(ev(k));
    v = V(:, k);

    % scale and rotate the eigenvector so c'qR = 1 and c'qI = 0, with c = real(v)
    vR = real(v);
    vI = imag(v);
    c = vR;
    ab = [c'*vR, -c'*vI; c'*vI, c'*vR] \ [1; 0];
    q = (ab(1) + 1i*ab(2))*v;
    qR0 = real(q);
    qI0 = imag(q);

    prob = hsys.make_problem(c);

    % land exactly on the Hopf curve: hold p2 fixed and solve for the rest
    land = @(z) [prob.residual(z); z(end) - p0(2)];
    z0 = fsolve(land, [u0; qR0; qI0; om0; p0], optimset('Display', 'off'));

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
    out.omega = pts(3*n+1, :);
    out.p = pts(3*n+2:3*n+3, :);
end
