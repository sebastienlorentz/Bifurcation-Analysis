lambda = 1.5;
mu = 0.1;
gamma = 1;
R0 = 1.205;

syms S I
f = 1/(1+I^2);
f0 = subs(f, I, 0); % treatment term at the disease-free equilibrium (I = 0)
beta = R0*mu*(mu + gamma + f0)/lambda;
F = [lambda - beta*S*I - mu*S + gamma*I;
     beta*S*I - (mu + gamma + f)*I];

u = [S; I];
xrange = [0 20];
yrange = [0 20];
f1 = matlabFunction(F(1), 'Vars', {S, I});
f2 = matlabFunction(F(2), 'Vars', {S, I});

plot_phase_portrait(f1, f2, F, u, xrange, yrange);
plot_phase_trajectories(f1, f2, F, u, xrange, yrange);

function plot_phase_portrait(f1, f2, F_sym, u_sym, xrange, yrange)
    [X, Y, U, V, speed] = unit_field(f1, f2, xrange, yrange, 150);

    figure;
    hold on;

    % speed heatmap: log scale so the slow regions near equilibria read as
    % distinct dark patches instead of being crushed by the fast corners
    logspeed = log10(speed + eps);
    pcolor(X, Y, logspeed);
    colormap(turbo)
    shading interp;
    cb = colorbar;
    cb.Label.String = 'log_{10} speed';
    % clamp to a few decades below the fastest speed, so a single near-zero
    % grid point (an equilibrium landing on a node) can't blow out the scale
    hi = max(logspeed(:));
    clim([hi - 4, hi]);

    % unit-length arrows (direction only); the heatmap carries the magnitude
    arrow_skip = 5; % draw every arrow_skip-th arrow so they don't swamp the map
    s = 1:arrow_skip:size(X, 1);
    quiver(X(s, s), Y(s, s), U(s, s), V(s, s), 0.2, 'w');

    plot_equilibria(F_sym, u_sym, xrange, yrange, 4);
    finalize_axes(u_sym, xrange, yrange, 'Phase portrait');
end

function plot_phase_trajectories(f1, f2, F_sym, u_sym, xrange, yrange)
    odefun = @(t, y) [f1(y(1), y(2)); f2(y(1), y(2))];

    figure;
    hold on;

    % seed a grid of initial conditions and integrate both ways in time;
    % trajectories flow into attractors and away from repellers
    nseed = 9;
    [Sx, Sy] = meshgrid(linspace(xrange(1), xrange(2), nseed), linspace(yrange(1), yrange(2), nseed));
    seeds = [Sx(:), Sy(:)];

    tspan = [0 500];
    opts = odeset('RelTol', 1e-6, ...
                  'Events', @(t, y) leave_box(t, y, xrange, yrange));
    for k = 1:size(seeds, 1)
        [~, Yf] = ode45(odefun, tspan, seeds(k, :)', opts);
        [~, Yb] = ode45(odefun, -tspan, seeds(k, :)', opts);
        plot(Yf(:, 1), Yf(:, 2), 'Color', [0.3 0.4 0.8]);
        plot(Yb(:, 1), Yb(:, 2), 'Color', [0.3 0.4 0.8]);
    end

    % direction field: normalized arrows on a coarse grid show flow direction
    [Ax, Ay, U, V] = unit_field(f1, f2, xrange, yrange, 10);
    quiver(Ax, Ay, U, V, 0.3, 'Color', [0.6 0.6 0.6], 'MaxHeadSize', 0.2);

    plot_equilibria(F_sym, u_sym, xrange, yrange, 5);
    finalize_axes(u_sym, xrange, yrange, 'Phase portrait (trajectories)');
end

function [X, Y, U, V, speed] = unit_field(f1, f2, xrange, yrange, n)
    % grid of unit-length flow vectors (NaN at equilibria) plus raw speed
    [X, Y] = meshgrid(linspace(xrange(1), xrange(2), n), ...
                      linspace(yrange(1), yrange(2), n));
    U = f1(X, Y).*ones(size(X));
    V = f2(X, Y).*ones(size(X));
    speed = hypot(U, V);
    U = U ./ speed;
    V = V ./ speed;
end

function plot_equilibria(F_sym, u_sym, xrange, yrange, markersize)
    % equilibria: solve F = 0, keep real in-window points, classify by Jacobian
    Jfun = matlabFunction(jacobian(F_sym, u_sym), 'Vars', num2cell(u_sym).');
    sol = solve(F_sym == 0, u_sym, 'Real', true);
    eqpts = double([sol.(char(u_sym(1))), sol.(char(u_sym(2)))]);
    in = imag(eqpts(:, 1)) == 0 & imag(eqpts(:, 2)) == 0 ...
       & eqpts(:, 1) >= xrange(1) & eqpts(:, 1) <= xrange(2) ...
       & eqpts(:, 2) >= yrange(1) & eqpts(:, 2) <= yrange(2);
    eqpts = real(eqpts(in, :));

    % a near-imaginary eigenvalue pair is a weak focus the linear test can't
    % resolve, so classify it by the sign of the first Lyapunov coefficient
    lsys = lyapunov_system(F_sym, u_sym, []);
    weak_tol = 0.02;

    cls = zeros(size(eqpts, 1), 1); % 1 stable, 2 unstable, 3 saddle, 4/5 weak focus
    for k = 1:size(eqpts, 1)
        ev = eig(Jfun(eqpts(k, 1), eqpts(k, 2)));
        re = real(ev);
        if ~isreal(ev) && max(abs(re)) < weak_tol
            l1 = lyapunov_coefficient(lsys, eqpts(k, :).', []);
            if l1 < 0
                cls(k) = 4;
            else
                cls(k) = 5;
            end
        elseif all(re < 0)
            cls(k) = 1;
        elseif all(re > 0)
            cls(k) = 2;
        else
            cls(k) = 3;
        end
    end

    % filled circle = stable, open circle = unstable, cross = saddle,
    % filled/open pentagram = stable/unstable weak focus
    specs = {1, 'o', 'k', 'stable'; 2, 'o', 'none', 'unstable'; 3, 'x', 'none', 'saddle'; ...
             4, 'p', 'k', 'weak focus (stable)'; 5, 'p', 'none', 'weak focus (unstable)'};
    hLeg = gobjects(0); labels = {};
    for c = 1:size(specs, 1)
        m = cls == specs{c, 1};
        if any(m)
            hLeg(end+1) = plot(eqpts(m, 1), eqpts(m, 2), specs{c, 2}, ...
                'MarkerFaceColor', specs{c, 3}, 'MarkerEdgeColor', 'k', ...
                'Color', 'k', 'MarkerSize', markersize, 'LineWidth', 1);
            labels{end+1} = specs{c, 4};
        end
    end
    if ~isempty(hLeg)
        legend(hLeg, labels, 'Location', 'best');
    end
end

function finalize_axes(u_sym, xrange, yrange, titleStr)
    xlabel(char(u_sym(1)));
    ylabel(char(u_sym(2)));
    title(titleStr);
    xlim(xrange);
    ylim(yrange);
    daspect([1 1 1]); % equal scale on both axes
end

function [val, term, dir] = leave_box(~, y, xr, yr)
% stop integrating once the trajectory leaves the (slightly padded) window
    pad = 0.05*[diff(xr), diff(yr)];
    inside = y(1) >= xr(1) - pad(1) && y(1) <= xr(2) + pad(1) && ...
             y(2) >= yr(1) - pad(2) && y(2) <= yr(2) + pad(2);
    val = double(inside) - 0.5; % crosses zero when the trajectory exits
    term = 1;
    dir = -1;
end
