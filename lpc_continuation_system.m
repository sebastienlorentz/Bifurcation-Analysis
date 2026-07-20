% Build the fold-of-cycles (LPC) continuation problem for two parameters. A fold
% of cycles is a periodic orbit whose nontrivial Floquet multiplier is +1.
% Tracking it as the parameter pair p = (p1, p2) varies traces an LPC curve in
% the (p1, p2) plane. The unknown is z = [x0; T; p1; p2] and the defining
% equations are
%       phi_T(x0) - x0 = 0        periodicity            (n)
%       x0(2) - I_sec  = 0        Poincare section       (1)
%       det(M(T)) - 1  = 0        nontrivial multiplier +1 (1)
% giving an (n+2) x (n+3) system; pseudo_arclength adds the arclength row. In 2-D
% the Floquet multipliers are {1 (trivial, the flow direction), det(M)}, so
% det(M)=1 pins the nontrivial multiplier WITHOUT having to separate it from the
% trivial one (an eigenvector condition cannot, since the two coincide at the
% LPC). M(T) and dx(T)/d(x0,p) come from the shared monodromy_integrate; the
% periodicity and section rows are analytic, the det row by finite differences.
% Two-parameter analogue of hopf_continuation_system.
function lcsys = lpc_continuation_system(F_sym, u_sym, p_sym)
    d = bifurcation_derivatives(F_sym, u_sym, p_sym);
    lcsys.n = d.n;
    lcsys.F = d.F;
    lcsys.Fu = d.Fu;
    lcsys.make_problem = @(I_sec) make_problem(I_sec, d);
end

% Wrap the residual and Jacobian (sharing one memoized monodromy integration)
% into the form pseudo_arclength expects, given the section level I_sec.
function prob = make_problem(I_sec, d)
    n = d.n;
    rhs  = @(x, p) d.F(x, p);
    jac  = @(x, p) d.Fu(x, p);
    dpar = @(x, p) d.Fp(x, p);     % n x 2

    zc = []; xT = []; M = []; P = []; % cached integration: phi_T, monodromy, dphi/dp
    prob.residual = @residual;
    prob.jacobian = @jacobian;

    function integrate(z)
        if isequal(z, zc)
            return
        end
        [xT, M, P] = monodromy_integrate(rhs, jac, dpar, n, 2, z(1:n), z(n+1), z(n+2:n+3));
        zc = z;
    end

    function R = residual(z)
        integrate(z);
        x0 = z(1:n);
        R = [xT - x0; x0(2) - I_sec; det(M) - 1];
    end

    % columns [x0 (n) | T | p1,p2]. Periodicity rows [M-I, F(x(T)), P] and the
    % constant section row are analytic; det(M)-1 has no cheap closed-form
    % derivative, so its row is forward finite differences against the memoized
    % base det(M), one fresh monodromy integration per column, not two. The
    % Jacobian only steers Newton, so the O(h) forward error costs at most a few
    % extra iterations; the residual stays exact and the final accuracy is unchanged.
    function J = jacobian(z)
        integrate(z);
        p = z(n+2:n+3);
        ephase = zeros(1, n); ephase(2) = 1;
        Jtop = [M - eye(n), rhs(xT, p), P;
                ephase,     0,          zeros(1, 2)];

        h = 1e-6;
        d0 = det(M);
        gd = zeros(1, n+3);
        for k = 1:n+3
            zp = z; zp(k) = zp(k) + h;
            [~, Mp] = monodromy_integrate(rhs, jac, dpar, n, 2, zp(1:n), zp(n+1), zp(n+2:n+3));
            gd(k) = (det(Mp) - d0)/h;
        end
        J = [Jtop; gd];
    end
end
