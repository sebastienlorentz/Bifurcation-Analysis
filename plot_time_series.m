% Integrate an n-D system and plot the state components in idx against time, one
% figure per initial condition. rhs is the vector field as a handle @(u)
% returning an n-vector; u0 is a 1-by-n row or an N-by-n matrix (one IC per row);
% u_sym labels the states. idx (optional) is the vector of state indices to draw
% (default all). Works in any dimension.
function [t, Y] = plot_time_series(rhs, u0, tspan, u_sym, idx)
    if nargin < 5 || isempty(idx)
        idx = 1:numel(u_sym);
    end
    odefun = @(t, y) rhs(y);
    % tight tolerance so the long near-saddle crawl of a near-homoclinic orbit
    % (where the period diverges) stays accurate
    o = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);

    labels = arrayfun(@(k) char(u_sym(k)), idx, 'UniformOutput', false);
    styles = {'-', '--', ':', '-.'}; % cycled so the states stay distinguishable

    for k = 1:size(u0, 1)
        [t, Y] = ode45(odefun, tspan, u0(k, :)', o);

        figure;
        hold on;
        for j = 1:numel(idx)
            plot(t, Y(:, idx(j)), 'LineWidth', 1);
        end
        xlabel('t');
        ylabel('state');
        legend(labels, 'Location', 'best');
        title(['Time series from u_0 = (', strjoin(compose('%g', u0(k, :)), ', '), ')']);
    end
end
