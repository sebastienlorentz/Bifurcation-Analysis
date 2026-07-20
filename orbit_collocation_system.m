% Build a periodic-orbit boundary-value problem by orthogonal collocation, the
% robust alternative to single shooting for continuing limit cycles (and, on top,
% folds of cycles). Time is rescaled to tau in [0,1] so the period T is an unknown:
%       x'(tau) = T * F(x, p),   x(0) = x(1),   phase condition.
% The orbit is a piecewise polynomial on Ntst mesh intervals, degree Ncol per
% interval (Lagrange basis through Ncol+1 equally spaced nodes), collocated at the
% Ncol Gauss-Legendre points of each interval. Unlike shooting there is no per-step
% ODE integration and no monodromy: the whole orbit is solved at once, and the
% Jacobian stays well conditioned through a fold of cycles (where shooting's M - I
% is singular).
%
% fld carries the vector field as handles of (x, p) with p an nparam-vector:
% fld.n, fld.nparam, fld.F, fld.Fu (n x n), fld.Fp (n x nparam). csys.make_problem
% (Xref) returns a problem (residual + analytic Jacobian) on z = [X(:); T; p], with
% X n-by-Np and Xref the reference profile fixing the integral phase condition.
% csys.seed_profile integrates one period to sample a starting profile.
function csys = orbit_collocation_system(fld, Ntst, Ncol)
    if nargin < 2 || isempty(Ntst); Ntst = 20; end
    if nargin < 3 || isempty(Ncol); Ncol = 4;  end
    n = fld.n;

    [gp, gw] = gauss_legendre(Ncol);          % collocation points/weights on [0,1]
    srep = linspace(0, 1, Ncol+1);            % representation nodes on [0,1]
    [A, B] = lagrange_matrices(srep, gp);     % values / derivatives at gp: Ncol x (Ncol+1)
    h = 1/Ntst;                               % uniform interval width in tau

    Np = Ntst*Ncol + 1;                       % global representation nodes
    csys.n = n;
    csys.nparam = fld.nparam;
    csys.Ntst = Ntst;
    csys.Ncol = Ncol;
    csys.Np = Np;
    csys.tau_nodes = linspace(0, 1, Np);      % global node positions in [0,1]
    % collocation internals exposed for the LPC extended system (second-derivative
    % blocks reuse the same basis matrices and quadrature)
    csys.A = A; csys.B = B; csys.gw = gw; csys.h = h;
    csys.make_problem = @(Xref) make_problem(Xref, fld, n, Ntst, Ncol, Np, A, B, gw, h);
    csys.seed_profile = @(x0, T, p) seed_profile(fld, x0, T, p, csys.tau_nodes);
end

function prob = make_problem(Xref, fld, n, Ntst, Ncol, Np, A, B, gw, h)
    np = fld.nparam;
    nX = n*Np;
    XrefDot = cell(1, Ntst);                  % reference dx/dtau at Gauss points
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        XrefDot{i} = (Xref(:, cols)*B.')/h;   % n x Ncol
    end
    prob.residual = @(z) residual(z, fld, n, np, Ntst, Ncol, Np, A, B, gw, h, Xref, XrefDot);
    prob.jacobian = @(z) jacobian(z, fld, n, np, Ntst, Ncol, Np, nX, A, B, gw, h, XrefDot);
end

function R = residual(z, fld, n, np, Ntst, Ncol, Np, A, B, gw, h, Xref, XrefDot)
    nX = n*Np;
    X = reshape(z(1:nX), n, Np);
    T = z(nX+1);
    p = z(nX+2:nX+1+np);
    R = zeros(Ntst*Ncol*n + n + 1, 1);
    row = 0;
    phase = 0;
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        Xi = X(:, cols);
        Xri = Xref(:, cols);
        for c = 1:Ncol
            xc = Xi*A(c, :).';
            xdot = (Xi*B(c, :).')/h;
            R(row+(1:n)) = xdot - T*fld.F(xc, p);
            row = row + n;
            phase = phase + gw(c)*h*(XrefDot{i}(:, c).'*(xc - Xri*A(c, :).'));
        end
    end
    R(row+(1:n)) = X(:, 1) - X(:, Np);        % periodicity
    R(end) = phase;                           % phase condition
end

function J = jacobian(z, fld, n, np, Ntst, Ncol, Np, nX, A, B, gw, h, XrefDot)
    X = reshape(z(1:nX), n, Np);
    T = z(nX+1);
    p = z(nX+2:nX+1+np);
    m = Ntst*Ncol*n + n + 1;
    J = zeros(m, nX+1+np);
    row = 0;
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        Xi = X(:, cols);
        for c = 1:Ncol
            xc = Xi*A(c, :).';
            Fc = fld.F(xc, p);
            Fu = fld.Fu(xc, p);
            Fp = fld.Fp(xc, p);               % n x np
            rr = row+(1:n);
            for j = 1:Ncol+1
                k = cols(j);
                cc = (k-1)*n + (1:n);
                J(rr, cc) = (B(c, j)/h)*eye(n) - T*Fu*A(c, j);
            end
            J(rr, nX+1) = -Fc;                % d/dT
            J(rr, nX+1+(1:np)) = -T*Fp;       % d/dp
            row = row + n;
        end
    end
    pr = row+(1:n);                           % periodicity
    J(pr, 1:n) = eye(n);
    J(pr, (Np-1)*n + (1:n)) = -eye(n);
    for i = 1:Ntst                            % phase (independent of T, p)
        cols = (i-1)*Ncol + (1:Ncol+1);
        for c = 1:Ncol
            for j = 1:Ncol+1
                k = cols(j);
                cc = (k-1)*n + (1:n);
                J(end, cc) = J(end, cc) + gw(c)*h*A(c, j)*XrefDot{i}(:, c).';
            end
        end
    end
end

% Sample a starting profile by integrating one period from a shooting-style seed.
function X = seed_profile(fld, x0, T, p, tau_nodes)
    o = odeset('RelTol', 1e-9, 'AbsTol', 1e-11);
    [~, Y] = ode45(@(t, y) fld.F(y, p), tau_nodes*T, x0(:), o);
    X = Y.';
end

% Gauss-Legendre nodes and weights on [0,1] (Golub-Welsch).
function [x, w] = gauss_legendre(m)
    if m == 1
        x = 0.5; w = 1; return
    end
    b = 0.5./sqrt(1 - (2*(1:m-1)).^(-2));
    Tm = diag(b, 1) + diag(b, -1);
    [V, D] = eig(Tm);
    [x, ix] = sort(diag(D));
    w = 2*(V(1, ix).^2);
    x = (x(:) + 1)/2;
    w = w(:)/2;
end

% Lagrange basis values A and derivatives B at points seval, for the basis through
% nodes srep. A(i,j) = L_j(seval_i), B(i,j) = L_j'(seval_i).
function [A, B] = lagrange_matrices(srep, seval)
    nn = numel(srep);
    ne = numel(seval);
    A = zeros(ne, nn);
    B = zeros(ne, nn);
    for j = 1:nn
        other = srep([1:j-1, j+1:nn]);
        denom = prod(srep(j) - other);
        for i = 1:ne
            x = seval(i);
            A(i, j) = prod(x - other)/denom;
            d = 0;
            for mm = 1:numel(other)
                rest = other([1:mm-1, mm+1:end]);
                d = d + prod(x - rest);
            end
            B(i, j) = d/denom;
        end
    end
end
