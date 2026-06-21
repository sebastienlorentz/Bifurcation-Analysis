% Build the augmented fold (saddle-node) continuation problem for two parameters.
% A fold is a point where F(u,p)=0 and the Jacobian Fu is singular
% The augmented unknown is z = [u; v; p], where v spans the null space of Fu,
% and the defining equations are
%       F(u,p) = 0, Fu(u,p) v = 0, c' v - 1 = 0
function fsys = fold_continuation_system(F_sym, u_sym, p_sym)
    d = bifurcation_derivatives(F_sym, u_sym, p_sym);
    fsys.n = d.n;
    fsys.F = d.F;
    fsys.Fu = d.Fu;
    fsys.make_problem = @(c) make_problem(c, d);
end

% Wrap the residual and Jacobian into the form pseudo_arclength expects.
function prob = make_problem(c, d)
    prob.residual = @(z) fold_residual(z, d, c);
    prob.jacobian = @(z) fold_jacobian(z, d, c);
end

% The defining equations of a fold, stacked.
function H = fold_residual(z, d, c)
    n = d.n;
    u = z(1:n);
    v = z(n+1:2*n);
    p = z(2*n+1:2*n+2);
    H = [d.F(u, p);
         d.Fu(u, p)*v;
         c'*v-1];
end

% Jacobian of the fold equations w.r.t. z = [u; v; p]
function J = fold_jacobian(z, d, c)
    n = d.n;
    u = z(1:n);
    v = z(n+1:2*n);
    p = z(2*n+1:2*n+2);

    Fu = d.Fu(u, p);

    % d(Fu*v)/du: column l is (d^2F/du du_l)*v
    Fuu = d.Fuu(u, p);
    B = zeros(n, n);
    for l = 1:n
        B(:, l) = Fuu(:, :, l)*v;
    end

    % d(Fu*v)/dp: one column per parameter
    Fvp = [d.Fup1(u, p)*v, d.Fup2(u, p)*v];
    J = [Fu, zeros(n, n), d.Fp(u, p);
         B, Fu, Fvp;
         zeros(1, n), c', zeros(1, 2)];
end
