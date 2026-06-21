% Build the transcritical continuation problem for two parameters. 
% We use the invariant structure directly: the curve is where the transverse 
% eigenvalue Fu(inv_idx,inv_idx) crosses zero on the subspace. The unknown is
% z = [u; p] and the defining equations are:
%   F(u,p) for the non-invariant components = 0 
%   u(inv_idx) = 0
%   Fu(inv_idx,inv_idx) = 0
function tsys = transcritical_continuation_system(F_sym, u_sym, p_sym, inv_idx)
    u_sym = u_sym(:);
    p_sym = p_sym(:);
    n = numel(u_sym);
    z_sym = [u_sym; p_sym];

    Fu_sym = jacobian(F_sym, u_sym);
    keep = setdiff(1:n, inv_idx); % equilibrium equations off the invariant subspace
    res_sym = [F_sym(keep); u_sym(inv_idx); Fu_sym(inv_idx, inv_idx)];
    jac_sym = jacobian(res_sym, z_sym);

    rf = matlabFunction(res_sym, 'Vars', {z_sym});
    jf = matlabFunction(jac_sym, 'Vars', {z_sym});

    tsys.n = n;
    tsys.F = matlabFunction(F_sym, 'Vars', {u_sym, p_sym});
    tsys.Fu = matlabFunction(Fu_sym, 'Vars', {u_sym, p_sym});
    tsys.make_problem = @() struct('residual', @(z) rf(z), 'jacobian', @(z) jf(z));
end
