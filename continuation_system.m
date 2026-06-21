% Build the function handles and derivatives the library needs from a symbolic
% 2-D ODE system, its state variables, and the continuation parameter.
function sys = continuation_system(F_sym, u_sym, lambda_sym)
    u_sym = u_sym(:);
    n = numel(u_sym);

    Fu_sym = jacobian(F_sym, u_sym);
    Fl_sym = diff(F_sym, lambda_sym);

    % second derivatives
    Fuu_flat_sym = jacobian(Fu_sym(:), u_sym); % (n*n) x n
    Ful_sym = diff(Fu_sym, lambda_sym);
    Fll_sym = diff(Fl_sym, lambda_sym);

    % third u-derivative, so sys also serves as a lyapunov_system (Fu/Fuu/Fuuu/n)
    Fuuu_flat_sym = jacobian(Fuu_flat_sym(:), u_sym); % (n*n*n) x n

    sys.n = n;
    sys.F = matlabFunction(F_sym, 'Vars', {u_sym, lambda_sym});
    sys.Fu = matlabFunction(Fu_sym, 'Vars', {u_sym, lambda_sym});
    sys.Flambda = matlabFunction(Fl_sym, 'Vars', {u_sym, lambda_sym});

    Fuu_flat_fn = matlabFunction(Fuu_flat_sym, 'Vars', {u_sym, lambda_sym});
    sys.Fuu = @(u, l) reshape(Fuu_flat_fn(u, l), [n, n, n]);
    Fuuu_flat_fn = matlabFunction(Fuuu_flat_sym, 'Vars', {u_sym, lambda_sym});
    sys.Fuuu = @(u, l) reshape(Fuuu_flat_fn(u, l), [n, n, n, n]);
    sys.Fulambda = matlabFunction(Ful_sym, 'Vars', {u_sym, lambda_sym});
    sys.Flambdalambda = matlabFunction(Fll_sym, 'Vars', {u_sym, lambda_sym});
end
