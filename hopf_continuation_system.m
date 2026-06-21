% Build the Hopf continuation problem for two parameters. A Hopf point
% is where Fu has a pure-imaginary eigenvalue pair +/- i*omega, with eigenvector
% q = qR + i*qI; tracking it as the 2-vector of parameters p varies traces a
% Hopf curve in the (p1,p2) plane. The unknown is z = [u; qR; qI; omega; p]
% and the defining equations are
%       F(u,p) = 0
%       Fu*qR + omega*qI = 0
%       Fu*qI - omega*qR = 0
%       c'*qR - 1 = 0
%       c'*qI = 0
function hsys = hopf_continuation_system(F_sym, u_sym, p_sym)
    d = bifurcation_derivatives(F_sym, u_sym, p_sym);
    hsys.n = d.n;
    hsys.F = d.F;
    hsys.Fu = d.Fu;
    hsys.make_problem = @(c) make_problem(c, d);
end

% Wrap the residual and Jacobian into the form pseudo_arclength expects.
function prob = make_problem(c, d)
    prob.residual = @(z) hopf_residual(z, d, c);
    prob.jacobian = @(z) hopf_jacobian(z, d, c);
end

% The defining equations of a Hopf point, stacked.
function H = hopf_residual(z, d, c)
    n = d.n;
    u = z(1:n);
    qR = z(n+1:2*n);
    qI = z(2*n+1:3*n);
    om = z(3*n+1);
    p = z(3*n+2:3*n+3);

    Fu = d.Fu(u, p);
    H = [d.F(u, p);
         Fu*qR + om*qI;
         Fu*qI - om*qR;
         c'*qR - 1;
         c'*qI];
end

% Jacobian of the Hopf equations w.r.t. z = [u; qR; qI; omega; p].
function J = hopf_jacobian(z, d, c)
    n = d.n;
    u = z(1:n);
    qR = z(n+1:2*n);
    qI = z(2*n+1:3*n);
    om = z(3*n+1);
    p = z(3*n+2:3*n+3);

    Fu = d.Fu(u, p);

    % d(Fu*qR)/du and d(Fu*qI)/du
    Fuu = d.Fuu(u, p);
    BR = zeros(n, n);
    BI = zeros(n, n);
    for l = 1:n
        BR(:, l) = Fuu(:, :, l)*qR;
        BI(:, l) = Fuu(:, :, l)*qI;
    end

    Fp1 = d.Fup1(u, p);
    Fp2 = d.Fup2(u, p);
    In = eye(n);
    Zn = zeros(n, n);

    J = [Fu, Zn, Zn, zeros(n, 1), d.Fp(u, p);
         BR, Fu, om*In, qI, [Fp1*qR, Fp2*qR];
         BI, -om*In, Fu, -qR, [Fp1*qI, Fp2*qI];
         zeros(1, n), c', zeros(1, n), 0, zeros(1, 2);
         zeros(1, n), zeros(1, n), c', 0, zeros(1, 2)];
end
