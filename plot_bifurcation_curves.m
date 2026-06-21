% Plot all traced fold and Hopf curves in the parameter plane.
function plot_bifurcation_curves(out, p_labels, p_box)
    figure;
    hold on;

    legend_handles = gobjects(0);
    legend_text = {};

    hf = plot_curve_set(out.fold_curves, 'k');
    if ~isempty(hf)
        legend_handles(end+1) = hf;
        legend_text{end+1} = 'Fold';
    end
    [hsup, hsub] = plot_hopf_set(out.hopf_curves, 'r');
    if ~isempty(hsup)
        legend_handles(end+1) = hsup;
        legend_text{end+1} = 'Hopf (supercritical)';
    end
    if ~isempty(hsub)
        legend_handles(end+1) = hsub;
        legend_text{end+1} = 'Hopf (subcritical)';
    end
    if isfield(out, 'transcritical_curves')
        ht = plot_curve_set(out.transcritical_curves, 'b');
        if ~isempty(ht)
            legend_handles(end+1) = ht;
            legend_text{end+1} = 'Transcritical';
        end
    end

    if isfield(out, 'bt_points') && ~isempty(out.bt_points)
        hb = plot(out.bt_points(1, :), out.bt_points(2, :), 'ks', ...
                  'MarkerFaceColor', 'm', 'MarkerSize', 7, 'LineStyle', 'none');
        legend_handles(end+1) = hb;
        legend_text{end+1} = 'Bogdanov-Takens';
    end
    if isfield(out, 'gh_points') && ~isempty(out.gh_points)
        hg = plot(out.gh_points(1, :), out.gh_points(2, :), 'o', ...
                  'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'c', ...
                  'MarkerSize', 7, 'LineStyle', 'none');
        legend_handles(end+1) = hg;
        legend_text{end+1} = 'Generalized Hopf';
    end

    if ~isempty(legend_handles)
        legend(legend_handles, legend_text, 'Location', 'best');
    end
    xlabel(p_labels{1});
    ylabel(p_labels{2});
    title('Bifurcation curves');
    xlim(p_box(1, :));
    ylim(p_box(2, :));
    grid on;
end

% Plot every curve in a set in one colour; return one handle for the legend.
function h = plot_curve_set(curves, color)
    h = gobjects(0);
    for i = 1:numel(curves)
        hi = plot(curves{i}.p(1, :), curves{i}.p(2, :), color, 'LineWidth', 1.5);
        if i == 1
            h = hi;
        else
            set(hi, 'HandleVisibility', 'off');
        end
    end
end

% Plot Hopf curves solid where supercritical (l1 < 0) and dashed where
% subcritical (l1 > 0); return one legend handle for each style. Falls back to a
% plain solid line if no l1 data is attached.
function [hsup, hsub] = plot_hopf_set(curves, color)
    hsup = gobjects(0);
    hsub = gobjects(0);
    for i = 1:numel(curves)
        c = curves{i};
        if ~isfield(c, 'l1')
            h = plot(c.p(1, :), c.p(2, :), color, 'LineWidth', 1.5);
            if isempty(hsup), hsup = h; else, set(h, 'HandleVisibility', 'off'); end
            continue
        end
        sup = c.l1 < 0;
        % extend the supercritical mask by one point each side, so the segment
        % spanning a sign change is still drawn - otherwise that
        % segment is dropped from both styles, leaving a gap at the transition
        supx = sup | [sup(2:end), false] | [false, sup(1:end-1)];
        Xs = c.p(1, :);
        Ys = c.p(2, :);
        Xs(~supx) = NaN;
        Ys(~supx) = NaN;
        Xu = c.p(1, :);
        Yu = c.p(2, :);
        Xu(sup) = NaN;
        Yu(sup) = NaN;
        h1 = plot(Xs, Ys, color, 'LineWidth', 1.5);
        h2 = plot(Xu, Yu, [color '--'], 'LineWidth', 1.5);
        if any(sup) && isempty(hsup), hsup = h1; else, set(h1, 'HandleVisibility', 'off'); end
        if any(~sup) && isempty(hsub), hsub = h2; else, set(h2, 'HandleVisibility', 'off'); end
    end
end
