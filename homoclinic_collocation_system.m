% Homoclinic orbit continuation by orthogonal collocation, for a saddle with a 1-D
% unstable manifold in any dimension. The orbit is truncated onto tau in [0,1] with T
% the truncation interval; it leaves the saddle q along the unstable eigenvector and
% returns inside the stable subspace:
%       x'(tau) - T*F(x, p) = 0                    collocation      (n*Ntst*Ncol)
%       F(q, p) = 0                                the saddle       (n)
%       Fu(q,p)  vu - lu vu = 0,  vu'vu - 1 = 0    right eigenpair  (n+1)
%       Fu(q,p)' wu - lw wu = 0,  wu'wu - 1 = 0    left  eigenpair  (n+1)
%       x(0) - q - e vu = 0                        leaves along W^u (n)
%       wu'(x(1) - q) = 0                          lands in  W^s    (1)
%       |x(1) - q| - e = 0                         at distance e    (1)
%       e - delta*|x(kfar) - q| = 0                truncation scale (1)
% The unknown is Z = [X(:); T; q; vu; lu; wu; lw; e; p], one longer than the residual,
% leaving the single free direction pseudo_arclength needs. The count works out to
% N - m = 1 in every dimension.
%
% Two choices carry the design. First, the eigenvectors are unknowns constrained by
% their own defining equations, so nothing ever differentiates an eig call and the
% Jacobian stays analytic. Second, the return condition uses the LEFT unstable
% eigenvector rather than a basis of the stable subspace: by biorthogonality
% wu'(x(1) - q) = 0 says exactly that x(1) - q lies in the stable subspace, as one
% scalar equation, whatever that subspace's dimension. That keeps the system
% dimension-general and, because the stable eigenvectors are never formed, costs
% nothing when the stable eigenvalues are complex (a saddle-focus, where a basis of
% real stable eigenvectors does not exist and an invariant-subspace/Riccati
% formulation would otherwise be needed).
%
% The truncation distance e is an unknown held at a fixed fraction delta of the loop
% amplitude, measured at the node kfar farthest from the saddle at landing. A fixed e
% would be wrong near a Bogdanov-Takens point, where the loop shrinks to nothing and
% a constant truncation ball ends up swallowing the whole orbit.
% ponytail: 1-D unstable manifold only. A higher-dimensional W^u needs an orthonormal
% basis in place of the single vu, and there the Riccati approach does come back.
function hsys = homoclinic_collocation_system(F_sym, u_sym, p_sym, Ntst, Ncol)
    if nargin < 4 || isempty(Ntst); Ntst = 20; end
    if nargin < 5 || isempty(Ncol); Ncol = 4;  end
    d = bifurcation_derivatives(F_sym, u_sym, p_sym);
    fld = struct('n', d.n, 'nparam', 2, 'F', d.F, 'Fu', d.Fu, 'Fp', d.Fp);
    csys = orbit_collocation_system(fld, Ntst, Ncol);
    n = d.n;

    hsys.n = n;
    hsys.Fu = d.Fu;                              % for the saddle quantity, post hoc
    hsys.csys = csys;
    hsys.Np = csys.Np;
    hsys.nX = n*csys.Np;
    hsys.nColl = csys.Ntst*csys.Ncol*n;          % interior collocation rows
    hsys.idx = layout(n, hsys.nX);
    hsys.make_problem = @(Xref, delta, kfar) make_problem(Xref, delta, kfar, d, csys, hsys.idx, hsys.nColl);
    hsys.seed_from_cycle = @(X, T, p, delta) seed_from_cycle(X, T, p, delta, d, csys.Np, hsys.idx);
end

% Where each block of unknowns sits in Z. p stays last: follow_curve reads the two
% continuation parameters off the end of the vector.
function ix = layout(n, nX)
    ix.X  = 1:nX;
    ix.T  = nX + 1;
    ix.q  = nX + 1 + (1:n);
    ix.vu = nX + 1 + n + (1:n);
    ix.lu = nX + 2 + 2*n;
    ix.wu = nX + 2 + 2*n + (1:n);                % left eigenvector of the same lu
    ix.lw = nX + 3 + 3*n;
    ix.e  = nX + 4 + 3*n;
    ix.p  = nX + 4 + 3*n + (1:2);
    ix.N  = nX + 3*n + 6;                        % total unknowns
end

function prob = make_problem(Xref, delta, kfar, d, csys, ix, nColl)
    % the underlying cycle problem supplies the interior collocation rows; its
    % periodicity and phase rows are sliced off and replaced by the saddle blocks
    orbit = csys.make_problem(Xref);
    prob.residual = @(Z) res(Z, orbit, d, csys.Np, ix, nColl, delta, kfar);
    prob.jacobian = @(Z) jac(Z, orbit, d, csys.Np, ix, nColl, delta, kfar);
    prob.idx = ix;
end

function R = res(Z, orbit, d, Np, ix, nColl, delta, kfar)
    n = numel(ix.q);
    X = reshape(Z(ix.X), n, Np);
    q = Z(ix.q); vu = Z(ix.vu); lu = Z(ix.lu); wu = Z(ix.wu); lw = Z(ix.lw);
    e = Z(ix.e); p = Z(ix.p);

    Rc = orbit.residual([Z(ix.X); Z(ix.T); p]);
    Jq = d.Fu(q, p);
    a = X(:, Np) - q;                            % arrival offset from the saddle
    R = [ Rc(1:nColl);                           % collocation interior
          d.F(q, p);                             % the saddle
          Jq*vu  - lu*vu;   vu.'*vu - 1;         % right eigenpair
          Jq.'*wu - lw*wu;  wu.'*wu - 1;         % left eigenpair
          X(:, 1) - q - e*vu;                    % departure along W^u
          wu.'*a;                                % arrival lies in W^s
          norm(a) - e;                           % at distance e
          e - delta*norm(X(:, kfar) - q) ];      % truncation tracks the loop size
end

function J = jac(Z, orbit, d, Np, ix, nColl, delta, kfar)
    n = numel(ix.q);
    q = Z(ix.q); vu = Z(ix.vu); lu = Z(ix.lu); wu = Z(ix.wu); lw = Z(ix.lw);
    e = Z(ix.e); p = Z(ix.p);
    nX = numel(ix.X);
    J = zeros(ix.N - 1, ix.N);
    I = eye(n);

    % --- interior collocation rows, lifted straight off the cycle Jacobian -------
    Jc = orbit.jacobian([Z(ix.X); Z(ix.T); p]);
    rr = 1:nColl;
    J(rr, ix.X) = Jc(rr, 1:nX);
    J(rr, ix.T) = Jc(rr, nX+1);
    J(rr, ix.p) = Jc(rr, nX+2:nX+3);

    Jq = d.Fu(q, p);
    Fuu = d.Fuu(q, p);
    Fup = cat(3, d.Fup1(q, p), d.Fup2(q, p));

    % --- the saddle: F(q,p) = 0 -------------------------------------------------
    r = nColl;
    J(r+(1:n), ix.q) = Jq;
    J(r+(1:n), ix.p) = d.Fp(q, p);
    r = r + n;

    % --- eigenpairs: (Fu - l I) v = 0 and (Fu' - l I) w = 0, both normalised -----
    [r, J] = eig_block(J, r, ix, ix.vu, ix.lu, vu, lu, Jq,     Fuu, Fup, n, false);
    [r, J] = eig_block(J, r, ix, ix.wu, ix.lw, wu, lw, Jq.',   Fuu, Fup, n, true);

    % --- boundary conditions ----------------------------------------------------
    X = reshape(Z(ix.X), n, Np);
    J(r+(1:n), ix.X(1:n))            =  I;       % x(0) = q + e*vu
    J(r+(1:n), ix.q)                 = -I;
    J(r+(1:n), ix.vu)                = -e*I;
    J(r+(1:n), ix.e)                 = -vu;
    r = r + n;

    a = X(:, Np) - q;
    iXend = ix.X(nX-n+1:nX);
    J(r+1, iXend) =  wu.';                       % wu'(x(1) - q) = 0
    J(r+1, ix.q)  = -wu.';
    J(r+1, ix.wu) =  a.';
    r = r + 1;

    ah = (a/norm(a)).';
    J(r+1, iXend) =  ah;                         % |x(1) - q| = e
    J(r+1, ix.q)  = -ah;
    J(r+1, ix.e)  = -1;
    r = r + 1;

    % --- truncation scale: e = delta*|x(kfar) - q| ------------------------------
    w = X(:, kfar) - q;
    wh = (w/norm(w)).';
    J(r+1, ix.X((kfar-1)*n + (1:n))) = -delta*wh;
    J(r+1, ix.q)                     =  delta*wh;
    J(r+1, ix.e)                     =  1;
end

% Rows for one eigenpair: (M - l I) v = 0 and v'v - 1 = 0, with M either Fu or its
% transpose. The q and p derivatives are where the field tensors earn their keep;
% transposing M means transposing each derivative slice too.
function [r, J] = eig_block(J, r, ix, iv, il, v, l, M, Fuu, Fup, n, transposed)
    rr = r + (1:n);
    for b = 1:n
        B = Fuu(:, :, b);
        if transposed; B = B.'; end
        J(rr, ix.q(b)) = B*v;                    % d(M v)/dq_b
    end
    J(rr, iv) = M - l*eye(n);
    J(rr, il) = -v;
    for b = 1:2
        B = Fup(:, :, b);
        if transposed; B = B.'; end
        J(rr, ix.p(b)) = B*v;
    end
    J(r+n+1, iv) = 2*v.';                        % normalisation
    r = r + n + 1;
end

% Convert a large-period cycle profile into a homoclinic seed: locate the saddle at
% the orbit's slowest point, split the loop there, and pin the ends along the
% eigenvectors. Returns Z0 plus the truncation fraction delta and the far node kfar
% the amplitude is measured at.
function [Z0, delta, kfar] = seed_from_cycle(X, T, p, delta, d, Np, ix)
    n = numel(ix.q);
    speed = vecnorm(cell2mat(arrayfun(@(k) d.F(X(:, k), p), 1:Np, 'UniformOutput', false)));
    [~, kmin] = min(speed);                      % slowest node = the saddle passage
    q = fsolve(@(x) d.F(x, p), X(:, kmin), optimoptions('fsolve', 'Display', 'off', ...
               'FunctionTolerance', 1e-14, 'StepTolerance', 1e-14));

    % the unstable eigenpair, right and left. Only the unstable direction is formed:
    % the stable subspace is handled implicitly through wu, so a complex stable pair
    % (saddle-focus) needs no special treatment.
    Jq = d.Fu(q, p);
    [V, D] = eig(Jq);
    ev = diag(D);
    [~, iu] = max(real(ev));
    lu = real(ev(iu));
    vu = real(V(:, iu)); vu = vu/norm(vu);
    [W, DW] = eig(Jq.');
    [~, iw] = min(abs(diag(DW) - ev(iu)));       % left vector for the same eigenvalue
    lw = real(DW(iw, iw));
    wu = real(W(:, iw)); wu = wu/norm(wu);

    % cut the closed loop at the saddle passage: node kmin becomes the last node, so
    % the arc departs just after closest approach and returns to it. X(:,1) and
    % X(:,Np) coincide on a cycle, so rotate the Np-1 distinct nodes.
    rot = mod((kmin:kmin+Np-2) - 1, Np-1) + 1;
    Xh = X(:, rot);
    Xh(:, Np) = Xh(:, 1);

    % orient vu along the departing arc, then pin both ends a distance e out
    if (Xh(:, 2) - q).'*vu < 0; vu = -vu; end
    if isempty(delta); delta = 0.05; end
    [amp, kfar] = max(vecnorm(Xh - q));          % far side of the loop
    e = delta*amp;
    Xh(:, 1) = q + e*vu;

    % the arrival end is placed by projecting the incoming point onto the stable
    % subspace (i.e. removing its wu component) and rescaling to length e. This needs
    % no stable eigenvector, so it works whatever the stable spectrum looks like.
    a = Xh(:, Np-1) - q;
    a = a - (wu.'*a)*wu;
    Xh(:, Np) = q + e*a/norm(a);

    Z0 = zeros(ix.N, 1);
    Z0(ix.X)  = Xh(:);
    Z0(ix.T)  = T;
    Z0(ix.q)  = q;
    Z0(ix.vu) = vu;  Z0(ix.lu) = lu;
    Z0(ix.wu) = wu;  Z0(ix.lw) = lw;
    Z0(ix.e)  = e;
    Z0(ix.p)  = p(:);
end
