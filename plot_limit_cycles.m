% Overlay one or more limit-cycle branches on the current equilibrium diagram:
% each cycle's max/min envelope vs the parameter, for the state components in idx,
% stable solid / unstable dashed, with the fold of cycles (LPC) marked. Call after
% plot_branches (with the SAME idx) so the cycles stack onto the same axes; every
% envelope is drawn in one cycle colour (green), distinct from the per-state
% equilibrium colours. lcs is a trace_limit_cycle struct or a cell array of
% them (e.g. from trace_limit_cycles). idx (optional) is the vector of state
% indices to draw (default all). state_labels (optional) is accepted so the call
% can mirror plot_branches, but is unused: all envelopes share one legend entry.
% Also opens a single companion period-vs-parameter figure, where T -> inf marks
% a homoclinic approach.
function plot_limit_cycles(lcs, param_label, idx, state_labels) %#ok<INUSD>
    if ~iscell(lcs)
        lcs = {lcs};
    end
    if isempty(lcs)
        return
    end
    n = size(lcs{1}.umax, 1);
    if nargin < 3 || isempty(idx)
        idx = 1:n;
    end
    green = [0 0.6 0]; % single cycle colour, distinct from the per-state equilibrium colours

    % overlay every envelope on the current (equilibrium) axes; only the first
    % envelope drawn contributes the single "limit cycle" legend entry
    ax = gca;
    hold(ax, 'on');
    for j = 1:numel(lcs)
        draw_envelope(ax, lcs{j}, idx, green, j == 1);
    end
    legend(ax);

    % one shared period figure for all the cycles
    figure;
    hold on;
    for j = 1:numel(lcs)
        draw_period(lcs{j}, green, j == 1);
    end
    xlabel(param_label);
    ylabel('period T');
    legend();
    grid on;
    title(['Limit-cycle period vs ', param_label]);
end

% Max/min envelope of the drawn components for one cycle, plus its LPC marker,
% all in the single cycle colour. label adds the one "limit cycle" legend entry,
% carried by the first envelope drawn.
function draw_envelope(ax, lc, idx, col, label)
    for j = 1:numel(idx)
        name = '';
        if label && j == 1
            name = 'limit cycle';
        end
        cycle_pair(ax, lc.lambda, lc.umax(idx(j), :), lc.umin(idx(j), :), lc.stable, col, name);
    end
    mark_lpc(ax, lc, idx, col, label);
end

% Period of one cycle on the current axes, stable solid / unstable dashed.
function draw_period(lc, col, label)
    p = lc.lambda;
    [Ts, Tu] = split_stable(lc.T, lc.stable);
    if label
        plot(p, Ts, 'Color', col, 'LineStyle', '-',  'LineWidth', 1.2, 'DisplayName', 'stable');
        plot(p, Tu, 'Color', col, 'LineStyle', '--', 'LineWidth', 1.2, 'DisplayName', 'unstable');
    else
        plot(p, Ts, 'Color', col, 'LineStyle', '-',  'LineWidth', 1.2, 'HandleVisibility', 'off');
        plot(p, Tu, 'Color', col, 'LineStyle', '--', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    end
    k = lpc_index(lc.lambda);
    if ~isempty(k)
        mark_pentagram(gca, lc.lambda(k), lc.T(k), col, label);
    end
end

% Draw a max/min envelope pair in one colour, solid where stable and dashed
% where unstable, sharing a single legend entry (empty name = no legend).
function cycle_pair(ax, p, hi, lo, st, col, name)
    [his, hiu] = split_stable(hi, st);
    [lis, lou] = split_stable(lo, st);
    if isempty(name)
        plot(ax, p, his, 'Color', col, 'LineStyle', '-', 'LineWidth', 1.1, 'HandleVisibility', 'off');
    else
        plot(ax, p, his, 'Color', col, 'LineStyle', '-', 'LineWidth', 1.1, 'DisplayName', name);
    end
    plot(ax, p, lis, 'Color', col, 'LineStyle', '-',  'LineWidth', 1.1, 'HandleVisibility', 'off');
    plot(ax, p, hiu, 'Color', col, 'LineStyle', '--', 'LineWidth', 1.1, 'HandleVisibility', 'off');
    plot(ax, p, lou, 'Color', col, 'LineStyle', '--', 'LineWidth', 1.1, 'HandleVisibility', 'off');
end

% Mark the fold of cycles (LPC) on every drawn component's envelope, as a
% pentagram edged in black in the cycle colour. Located as the parameter turning
% point along the branch; label adds the single legend entry (first marker only).
function mark_lpc(ax, lc, idx, col, label)
    k = lpc_index(lc.lambda);
    if isempty(k)
        return
    end
    first = label;
    for j = 1:numel(idx)
        for y = [lc.umax(idx(j), k), lc.umin(idx(j), k)]
            mark_pentagram(ax, lc.lambda(k), y, col, first);
            first = false;
        end
    end
end

% One LPC pentagram; label true adds the (single) legend entry.
function mark_pentagram(ax, x, y, col, label)
    if label
        plot(ax, x, y, 'p', 'MarkerSize', 7, 'MarkerFaceColor', col, ...
             'MarkerEdgeColor', 'k', 'DisplayName', 'LPC (fold of cycles)');
    else
        plot(ax, x, y, 'p', 'MarkerSize', 7, 'MarkerFaceColor', col, ...
             'MarkerEdgeColor', 'k', 'HandleVisibility', 'off');
    end
end

% Index of the parameter turning point (fold) along the branch, if any: where the
% step in the parameter changes sign.
function k = lpc_index(p)
    dp = diff(p);
    k = find(dp(1:end-1).*dp(2:end) < 0, 1) + 1;
end
