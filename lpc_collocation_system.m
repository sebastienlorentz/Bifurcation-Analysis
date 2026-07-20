% Fold-of-cycles (LPC) continuation by orthogonal collocation, the robust
% replacement for the single-shooting lpc_continuation_system. A fold of cycles is
% a periodic orbit at which the cycle Jacobian G_W (of the collocation BVP w.r.t
% the orbit unknowns W = [X; T]) is singular. Rather than the shooting condition
% det(M) = 1, which makes the periodicity block M - I singular and wrecks the
% corrector, so we continue the extended system with an explicit null vector V:
%       R(W, p)      = 0     cycle collocation BVP            (nW)
%       G_W(W, p) V  = 0     V spans ker G_W (the fold)       (nW)
%       V'V - 1      = 0     normalisation                    (1)
% in the two parameters p = (p1, p2). The unknown is Z = [W; p; V]. Everything is
% algebraic (no ODE integration) and the Jacobian is fully analytic. The second-
% derivative blocks d(G_W V)/dW and d(G_W V)/dp come from the field tensors in d
% (bifurcation_derivatives: Fuu, Fup1, Fup2). make_problem(Xref) returns the
% pseudo_arclength problem; prob.orbit is the underlying cycle problem, used to land
% the seed and to size things (nW = n*Np + 1).
function lsys = lpc_collocation_system(F_sym, u_sym, p_sym, Ntst, Ncol)
    if nargin < 4 || isempty(Ntst); Ntst = 20; end
    if nargin < 5 || isempty(Ncol); Ncol = 4;  end
    d = bifurcation_derivatives(F_sym, u_sym, p_sym);
    fld = struct('n', d.n, 'nparam', 2, 'F', d.F, 'Fu', d.Fu, 'Fp', d.Fp);
    csys = orbit_collocation_system(fld, Ntst, Ncol);
    lsys.n = d.n;
    lsys.csys = csys;
    lsys.Np = csys.Np;
    lsys.nW = d.n*csys.Np + 1;
    lsys.seed_profile = csys.seed_profile;
    lsys.make_problem = @(Xref) make_problem(Xref, d, csys);
end

function prob = make_problem(Xref, d, csys)
    n = d.n; Np = csys.Np; nX = n*Np; nW = nX + 1;
    A = csys.A; B = csys.B; gw = csys.gw; h = csys.h; Ntst = csys.Ntst; Ncol = csys.Ncol;
    orbit = csys.make_problem(Xref);
    XrefDot = cell(1, Ntst);                  % reference dx/dtau at Gauss points
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        XrefDot{i} = (Xref(:, cols)*B.')/h;
    end
    prob.residual = @(Z) res(Z, orbit, d, n, Np, nX, nW, Ntst, Ncol, A, B, gw, h, XrefDot);
    prob.jacobian = @(Z) jac(Z, orbit, d, n, Np, nX, nW, Ntst, Ncol, A);
    prob.orbit = orbit;
    prob.nW = nW;
end

function R = res(Z, orbit, d, n, Np, nX, nW, Ntst, Ncol, A, B, gw, h, XrefDot)
    W = Z(1:nW); V = Z(nW+1:2*nW); p = Z(2*nW+1:2*nW+2);
    Rc = orbit.residual([W; p]);
    b2 = GwV(W, p, V, d, n, Np, nX, Ntst, Ncol, A, B, gw, h, XrefDot);
    R = [Rc; b2; V.'*V - 1];
end

% Directly evaluate G_W(W,p) * V (the linearised cycle residual applied to V),
% matching the analytic Jacobian's first nW columns, cheaper than assembling G_W.
function b2 = GwV(W, p, V, d, n, Np, nX, Ntst, Ncol, A, B, gw, h, XrefDot)
    X = reshape(W(1:nX), n, Np); T = W(nX+1);
    Phi = reshape(V(1:nX), n, Np); phiT = V(nX+1);
    b2 = zeros(Ntst*Ncol*n + n + 1, 1);
    row = 0; ph = 0;
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        Xi = X(:, cols); Phii = Phi(:, cols);
        for c = 1:Ncol
            xc = Xi*A(c, :).';
            Phic = Phii*A(c, :).';
            Phidot = (Phii*B(c, :).')/h;
            b2(row+(1:n)) = Phidot - T*d.Fu(xc, p)*Phic - d.F(xc, p)*phiT;
            row = row + n;
            ph = ph + gw(c)*h*(XrefDot{i}(:, c).'*Phic);
        end
    end
    b2(row+(1:n)) = Phi(:, 1) - Phi(:, Np);   % periodicity
    b2(end) = ph;                             % phase
end

function J = jac(Z, orbit, d, n, Np, nX, nW, Ntst, Ncol, A)
    W = Z(1:nW); V = Z(nW+1:2*nW); p = Z(2*nW+1:2*nW+2);
    Jc = orbit.jacobian([W; p]);              % nW x (nW+2)
    Gw = Jc(:, 1:nW);
    Rp = Jc(:, nW+1:nW+2);
    [GwwV, GwpV] = second_blocks(W, p, V, d, n, Np, nX, nW, Ntst, Ncol, A);
    Z0 = zeros(nW, nW);
    % columns [ W | V | p ];  rows [ cycle BVP | G_W V | normalisation ]
    J = [ Gw,          Z0,       Rp;
          GwwV,        Gw,       GwpV;
          zeros(1, nW), 2*V.',   zeros(1, 2) ];
end

% Analytic second-derivative blocks: d(G_W V)/dW and d(G_W V)/dp. Only the
% collocation rows are nonzero (periodicity and phase rows of G_W are constant in
% W and p). Uses the field tensors Fuu, Fup1, Fup2.
function [GwwV, GwpV] = second_blocks(W, p, V, d, n, Np, nX, nW, Ntst, Ncol, A)
    X = reshape(W(1:nX), n, Np); T = W(nX+1);
    Phi = reshape(V(1:nX), n, Np); phiT = V(nX+1);
    m = Ntst*Ncol*n + n + 1;
    GwwV = zeros(m, nW);
    GwpV = zeros(m, 2);
    row = 0;
    for i = 1:Ntst
        cols = (i-1)*Ncol + (1:Ncol+1);
        Xi = X(:, cols); Phii = Phi(:, cols);
        for c = 1:Ncol
            xc = Xi*A(c, :).';
            Phic = Phii*A(c, :).';
            Fu = d.Fu(xc, p);
            Fuu = d.Fuu(xc, p);
            Dv = zeros(n, n);                 % Dv(:,b) = (dFu/dx_b) Phic
            for b = 1:n
                Dv(:, b) = Fuu(:, :, b)*Phic;
            end
            blk = -(T*Dv + phiT*Fu);          % node multiplier for d/dX
            rr = row+(1:n);
            for jj = 1:Ncol+1
                k = cols(jj);
                ccX = (k-1)*n + (1:n);
                GwwV(rr, ccX) = A(c, jj)*blk;
            end
            GwwV(rr, nW) = -Fu*Phic;          % d/dT
            Fp = d.Fp(xc, p);
            GwpV(rr, 1) = -T*d.Fup1(xc, p)*Phic - phiT*Fp(:, 1);
            GwpV(rr, 2) = -T*d.Fup2(xc, p)*Phic - phiT*Fp(:, 2);
            row = row + n;
        end
    end
end
