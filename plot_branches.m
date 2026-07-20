% Plot equilibrium curves vs the parameter for the state components in idx:
% stable solid, unstable dashed, folds/Hopfs/branch points marked. out is the
% trace_branches struct. idx (optional) is the vector of state indices to draw
% (default all); state_labels (optional) names them in idx order and defaults to
% x_i. Call again with a different idx to draw other components on a fresh figure.
function plot_branches(out, param_label, param_range, idx, state_labels)
    n = size(out.branches{1}, 1) - 1; % number of state components
    if nargin < 4 || isempty(idx)
        idx = 1:n;
    end
    if nargin < 5 || isempty(state_labels)
        state_labels = arrayfun(@(k) sprintf('x_%d', k), idx, 'UniformOutput', false);
    end
    colors = lines(numel(idx)); % one colour per drawn component

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
        p = pts(end, :); % the parameter values along the branch

        % split each curve into stable (solid) and unstable (dashed) parts
        for j = 1:numel(idx)
            [xs, xu] = split_stable(pts(idx(j), :), st);
            if first
                plot(p, xs, 'Color', colors(j, :), 'LineWidth', 1.5, 'DisplayName', state_labels{j});
            else
                plot(p, xs, 'Color', colors(j, :), 'LineWidth', 1.5, 'HandleVisibility', 'off');
            end
            plot(p, xu, 'Color', colors(j, :), 'LineStyle', '--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
        end
        first = false;

        mark_points(out.folds{i}, 'ko', 'k', 'Fold', idx);
        if have_l1
            hopf_pts = [hopf_pts, out.hopfs{i}];
            hopf_l1 = [hopf_l1, out.hopf_l1{i}];
        else
            mark_points(out.hopfs{i}, 'ks', 'g', 'Hopf', idx);
        end
    end
    mark_hopf_criticality(hopf_pts, hopf_l1, idx);
    mark_points(out.branch_points, 'kd', 'm', 'Branch point', idx);

    xlabel(param_label);
    legend();
    grid on;
    title(['Equilibrium curves of ', strjoin(state_labels, ', '), ' vs ', param_label]);
    xlim(param_range);

    % pad the y-axis just past the drawn components' extremes
    vals = [out.branches{:}];
    vals = vals(idx, :);
    pad = 0.05*(max(vals(:)) - min(vals(:)));
    ylim([min(vals(:)) - pad, max(vals(:)) + pad]);
end

% Mark Hopf points on every drawn component, filled for supercritical (l1 < 0,
% stable cycle born) and open for subcritical (l1 > 0), with one legend entry each.
function mark_hopf_criticality(hopf_pts, hopf_l1, idx)
    if isempty(hopf_pts)
        return
    end
    sup = hopf_l1 < 0; % filled marker supercritical, open subcritical
    if any(sup)
        mark_points(hopf_pts(:, sup), 'ks', 'g', 'Hopf (supercritical)', idx);
    end
    if any(~sup)
        mark_points(hopf_pts(:, ~sup), 'ks', 'none', 'Hopf (subcritical)', idx);
    end
end

% Mark bifurcation points on every drawn component with a single legend entry.
function mark_points(P, style, face_color, name, idx)
    if isempty(P)
        return
    end
    h = gobjects(1, numel(idx));
    for j = 1:numel(idx)
        h(j) = plot(P(end, :), P(idx(j), :), style, 'MarkerFaceColor', face_color);
    end
    legend_once(h, name);
end
