% Build the state-derivative handles needed for the first Lyapunov coefficient:
% the Jacobian Fu and the second/third derivative tensors Fuu and Fuuu, as handles of (u, p).
function lsys = lyapunov_system(F_sym, u_sym, p_sym)
    u_sym = u_sym(:);
    p_sym = p_sym(:);
    n = numel(u_sym);

    Fu_sym = jacobian(F_sym, u_sym);
    Fuu_flat_sym = jacobian(Fu_sym(:), u_sym);
    Fuuu_flat_sym = jacobian(Fuu_flat_sym(:), u_sym);

    lsys.n = n;
    fu = matlabFunction(Fu_sym, 'Vars', {u_sym, p_sym});
    fuu = matlabFunction(Fuu_flat_sym, 'Vars', {u_sym, p_sym});
    fuuu = matlabFunction(Fuuu_flat_sym, 'Vars', {u_sym, p_sym});
    lsys.Fu = @(u, p) fu(u, p);
    lsys.Fuu = @(u, p) reshape(fuu(u, p), [n, n, n]);
    lsys.Fuuu = @(u, p) reshape(fuuu(u, p), [n, n, n, n]);
end
