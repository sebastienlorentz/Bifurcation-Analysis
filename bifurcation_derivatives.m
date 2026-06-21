% Function handles shared by the two-parameter fold and Hopf systems: F, Jacobian
% Fu, parameter derivative Fp, second derivatives Fuu and d(Fu)/dp, all as handles
% of (u, p) with p a 2-vector, plus the state dimension n.
function d = bifurcation_derivatives(F_sym, u_sym, p_sym)
    u_sym = u_sym(:);
    p_sym = p_sym(:);
    n = numel(u_sym);

    Fu_sym = jacobian(F_sym, u_sym);
    Fp_sym = jacobian(F_sym, p_sym); % n x 2: dF/dp
    Fuu_flat_sym = jacobian(Fu_sym(:), u_sym); % (n*n) x n: d(Fu)/du
    Fup1_sym = diff(Fu_sym, p_sym(1)); % n x n: d(Fu)/dp1
    Fup2_sym = diff(Fu_sym, p_sym(2)); % n x n: d(Fu)/dp2

    d.n = n;
    d.F = matlabFunction(F_sym, 'Vars', {u_sym, p_sym});
    d.Fu = matlabFunction(Fu_sym, 'Vars', {u_sym, p_sym});
    d.Fp = matlabFunction(Fp_sym, 'Vars', {u_sym, p_sym});
    Fuu_flat = matlabFunction(Fuu_flat_sym, 'Vars', {u_sym, p_sym});
    d.Fuu = @(u, p) reshape(Fuu_flat(u, p), [n, n, n]);
    d.Fup1 = matlabFunction(Fup1_sym, 'Vars', {u_sym, p_sym});
    d.Fup2 = matlabFunction(Fup2_sym, 'Vars', {u_sym, p_sym});
end
