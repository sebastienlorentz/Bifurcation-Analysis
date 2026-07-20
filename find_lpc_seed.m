% Find a fold-of-cycles (LPC) seed for trace_lpc_curve. At each fixed p2 slice it
% builds the one-parameter system, traces the equilibria and the limit cycles
% born at their Hopfs (the same chain as example_one_parameter_R0), and returns
% the first orbit sitting at a genuine fold of cycles, where the nontrivial
% Floquet multiplier crosses +1. Returns the section point x0, period T0 and
% parameter pair p0 = [p1; p2] of that orbit, the seed the two-parameter LPC
% continuation needs. Sweeps the slices in p2_vals order, so pass them ordered
% with the robust ones first: the fold is cleanest well inside the subcritical
% region and fragile right at the generalized Hopf point, where the
% cycle amplitude vanishes. Continuation then carries the curve back to the Bautin.
function [x0, T0, p0] = find_lpc_seed(F_sym, u_sym, p_sym, p2_vals, u_start, p1_range, state_box)
    for p2 = p2_vals
        sys = continuation_system(subs(F_sym, p_sym(2), p2), u_sym, p_sym(1));

        % tracing the equilibria can fail outright on a slice (e.g. no physical
        % branch); skip it and try the next rather than aborting the search
        try
            br = trace_branches(sys, u_start, p1_range(2), p1_range, 1e-4, 1e6, true, state_box);
        catch e
            fprintf('  p2 = %.3g: branch tracing failed (%s)\n', p2, e.message);
            continue
        end

        % continue each Hopf's cycle on its own, inside its own try, so one
        % fragile branch (a near-homoclinic blow-up, a shooting break-down) is
        % skipped without discarding the whole slice and any fold it may hold
        H  = [br.hopfs{:}];    % [u; p1] per Hopf, flattened across branches
        l1 = [br.hopf_l1{:}];
        for k = 1:size(H, 2)
            try
                % stop_on_fold: end each arc at the first +1 Floquet crossing, so a
                % slice that has a fold returns in a handful of points instead of
                % tracing the whole long-period branch both ways first.
                lc = trace_limit_cycle(sys, H(1:end-1, k), H(end, k), l1(k), ...
                                       [], [], [], [], [], true);
            catch
                continue
            end
            j = fold_index(lc);
            if ~isempty(j)
                z = lc.z(:, j); % [x0; T; p1]
                x0 = z(1:end-2);
                T0 = z(end-1);
                p0 = [z(end); p2];
                fprintf('  p2 = %.3g: fold of cycles at p1 = %.4g (nontrivial multiplier = 1)\n', ...
                        p2, p0(1));
                return
            end
        end
        fprintf('  p2 = %.3g: %d Hopf(s), no genuine fold of cycles on this slice\n', ...
                p2, size(H, 2));
    end
    error('find_lpc_seed:none', 'no fold of cycles found in the supplied p2 slices');
end

% Index of a genuine fold of cycles: where the nontrivial Floquet multiplier
% crosses +1 in the branch interior. A fold of cycles is DEFINED by this crossing
% (a saddle-node of periodic orbits). The plain p1-turning test is not enough: it
% is fooled by near-homoclinic shooting break-down, which reverses p1 without the
% multiplier reaching +1, so it seeds on a non-fold and the LPC continuation can't
% start. In 2-D one multiplier is trivially 1 (the flow direction); the other is
% the nontrivial one whose crossing of +1 marks the fold.
function j = fold_index(lc)
    M = lc.multipliers;
    nt = zeros(1, size(M, 2));
    for c = 1:numel(nt)
        m = M(:, c);
        [~, it] = min(abs(m - 1));        % the trivial multiplier (nearest 1)
        r = m([1:it-1, it+1:end]);
        nt(c) = real(r(1));               % the nontrivial multiplier
    end
    g = nt - 1;
    j = find(g(1:end-1).*g(2:end) < 0, 1); % first sign change of (nt - 1) = fold
    if ~isempty(j); j = j + 1; end
end
