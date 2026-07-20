% First Lyapunov coefficient at a Hopf point, via Kuznetsov's projection formula
% Its sign fixes the criticality of the Hopf:
%   l1 < 0 supercritical: a stable limit cycle is born
%   l1 > 0 subcritical: an unstable cycle surrounds the equilibrium
%   l1 = 0 degenerate: a generalized Hopf point
function l1 = lyapunov_coefficient(lsys, u, p)
    u = u(:);
    n = lsys.n;
    A = lsys.Fu(u, p);
    Fuu = lsys.Fuu(u, p);
    Fuuu = lsys.Fuuu(u, p);

    % critical eigenpair: the eigenvalue nearest the imaginary axis with
    % positive frequency; q is its right eigenvector, A q = i*omega q
    [V, D] = eig(A);
    ev = diag(D);
    cand = find(imag(ev) > 0);
    if isempty(cand) % no imaginary pair, not a Hopf point
        l1 = NaN;
        return
    end
    [~, idx] = min(abs(real(ev(cand))));
    k = cand(idx);
    omega = imag(ev(k));
    q = V(:, k);
    q = q / norm(q);

    % left eigenvector p_l: A' p_l = -i*omega p_l, normalised so p_l'* q = 1
    [W, DT] = eig(A.');
    evt = diag(DT);
    [~, j] = min(abs(evt + 1i*omega));
    pl = W(:, j);
    pl = pl / conj(pl'*q);

    % multilinear forms (non-conjugate contractions of the derivative tensors)
    qb = conj(q);
    h11 = A \ bilin(Fuu, q, qb);
    h20 = (2i*omega*eye(n)-A) \ bilin(Fuu, q, q);
    g21 = pl'*trilin(Fuuu, q, q, qb) - 2*(pl'*bilin(Fuu, q, h11)) + pl'*bilin(Fuu, qb, h20);
    l1 = real(g21) / (2*omega);
end

% Trilinear form C(a,b,c)_i = sum_{j,k,l} Fuuu(i,j,k,l) a_j b_k c_l.
function w = trilin(T, a, b, c)
    n = size(T, 1);
    w = zeros(n, 1);
    for i = 1:n
        Ti = reshape(T(i, :, :, :), [n, n, n]);
        s = zeros(n, 1);
        for l = 1:n
            s(l) = a.'*Ti(:, :, l)*b;
        end
        w(i) = s.'*c;
    end
end
