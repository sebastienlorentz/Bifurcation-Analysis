% Plot S and I equilibrium curves vs the parameter: stable solid, unstable dashed,
% folds/Hopfs/branch points marked. out is the trace_branches struct.
function plot_branches(out, param_label, param_range)
    figure;
    hold on;

    have_l1 = isfield(out, 'hopf_l1'); % style Hopfs by criticality if available
    hopf_pts = [];
    hopf_l1 = [];

    first = true;
    for i = 1:numel(out.branches)
        pts = out.branches{i};
        if size(pts, 2) < 2 % a single point can't be drawn as a line
            continue
        end
        st = out.stable{i};

        % decimate very long curves so the renderer can cope; a few thousand
        % points is plenty to draw a smooth line
        max_plot_pts = 5000;
        if size(pts, 2) > max_plot_pts
            idx = round(linspace(1, size(pts, 2), max_plot_pts));
            pts = pts(:, idx);
            st = st(idx);
        end

        p = pts(end, :); % the parameter values along the branch

        % split each curve into stable (solid) and unstable (dashed) parts by
        % masking the other part out with NaN. The stable mask is extended by one
        % point each side so the segment at a stability change (Hopf/fold) is
        % still drawn, rather than being dropped from both styles and leaving a gap
        stx = st | [st(2:end), false] | [false, st(1:end-1)];
        Is = pts(2, :);
        Iu = pts(2, :);
        Ss = pts(1, :);
        Su = pts(1, :);
        Is(~stx) = NaN;
        Iu(st) = NaN;
        Ss(~stx) = NaN;
        Su(st) = NaN;

        if first
            plot(p, Is, 'r', 'LineWidth', 1.5, 'DisplayName', 'I (Infected)');
            plot(p, Ss, 'b', 'LineWidth', 1.5, 'DisplayName', 'S (Susceptible)');
            first = false;
        else
            plot(p, Is, 'r', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            plot(p, Ss, 'b', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        plot(p, Iu, 'r--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        plot(p, Su, 'b--', 'LineWidth', 1.5, 'HandleVisibility', 'off');

        mark_points(out.folds{i}, 'ko', 'k', 'Fold');
        if have_l1
            hopf_pts = [hopf_pts, out.hopfs{i}]; 
            hopf_l1 = [hopf_l1, out.hopf_l1{i}]; 
        else
            mark_points(out.hopfs{i}, 'ks', 'g', 'Hopf');
        end
    end
    mark_hopf_criticality(hopf_pts, hopf_l1);
    mark_points(out.branch_points, 'kd', 'm', 'Branch point');

    xlabel(param_label);
    legend();
    grid on;
    title(['Equilibrium curves of S and I vs ', param_label]);
    xlim(param_range);

    % pad the y-axis just past the S/I extremes
    vals = [out.branches{:}];
    vals = vals(1:2, :);
    pad = 0.05*(max(vals(:)) - min(vals(:)));
    ylim([min(vals(:)) - pad, max(vals(:)) + pad]);
end

% Mark Hopf points on both curves, filled for supercritical (l1 < 0, stable cycle
% born) and open for subcritical (l1 > 0), with one legend entry each.
function mark_hopf_criticality(hopf_pts, hopf_l1)
    if isempty(hopf_pts)
        return
    end
    p = hopf_pts(end, :);
    sup = hopf_l1 < 0;
    if any(sup)
        plot(p(sup), hopf_pts(2, sup), 'ks', 'MarkerFaceColor', 'g', 'DisplayName', 'Hopf (supercritical)');
        plot(p(sup), hopf_pts(1, sup), 'ks', 'MarkerFaceColor', 'g', 'HandleVisibility', 'off');
    end
    if any(~sup)
        plot(p(~sup), hopf_pts(2, ~sup), 'ks', 'MarkerFaceColor', 'none', 'DisplayName', 'Hopf (subcritical)');
        plot(p(~sup), hopf_pts(1, ~sup), 'ks', 'MarkerFaceColor', 'none', 'HandleVisibility', 'off');
    end
end

% Mark bifurcation points on both curves with a single legend entry.
function mark_points(P, style, face_color, name)
    if isempty(P)
        return
    end
    plot(P(end, :), P(2, :), style, 'MarkerFaceColor', face_color, 'DisplayName', name);
    plot(P(end, :), P(1, :), style, 'MarkerFaceColor', face_color, 'HandleVisibility', 'off');
end
