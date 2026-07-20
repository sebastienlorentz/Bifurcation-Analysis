% Contract the first index of a 3-D tensor's slices with two vectors:
% w(i) = a.' * T(i,:,:) * b = sum_{j,k} T(i,j,k) a_j b_k (no conjugation).
% Shared by the branch-point tangents (trace_branches) and the first Lyapunov
% coefficient (lyapunov_coefficient), which both need this bilinear form.
function w = bilin(T, a, b)
    n = size(T, 1);
    w = zeros(n, 1);
    for i = 1:n
        w(i) = a.' * reshape(T(i, :, :), n, n) * b;
    end
end
