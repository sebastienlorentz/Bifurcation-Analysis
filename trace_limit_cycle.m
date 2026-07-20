% Continue the limit cycle born at a Hopf point in one parameter, by single
% shooting. Seeds the first orbit by integrating onto it (forward for a stable
% cycle, reverse for an unstable one), refines with Newton, then follows it with
% pseudo_arclength on the shooting_system. l1 is the first Lyapunov coefficient
% at the Hopf (its sign sets which side the cycle lives and the cycle's
% stability); dlambda is the initial parameter offset. Returns the max/min
% envelope of each state, the period, and Floquet multipliers along the branch,
% ready to overlay on the equilibrium diagram from trace_branches.
function out = trace_limit_cycle(sys, u_H, lambda_H, l1, dlambda, step_size, max_steps, lambda_range, T_cap, stop_on_fold)
    n = sys.n;
    u_H = u_H(:);
    if nargin < 5 || isempty(dlambda);      dlambda = 5e-3;            end
    if nargin < 6 || isempty(step_size);    step_size = 1e-2;          end
    if nargin < 7 || isempty(max_steps);    max_steps = 1e4;           end
    if nargin < 8 || isempty(lambda_range); lambda_range = [-inf inf]; end
    % period-growth cap as a multiple of the seed period T0: stop once the cycle
    % has stretched this far toward the homoclinic. The default 50 traces the full
    % branch (for display); a seed search that only needs the fold passes a small
    % value to skip the slow long-period tail.
    if nargin < 9 || isempty(T_cap);        T_cap = 50;                end
    % seed search only: stop the arc the instant the nontrivial Floquet multiplier
    % crosses +1 (a fold of cycles), instead of tracing the whole ~800-point branch
    % and both directions before scanning. Cuts a successful seed from minutes to
    % the handful of points up to the first fold. Off by default (full branch).
    if nargin < 10 || isempty(stop_on_fold); stop_on_fold = false;     end

    % side where the cycle lives: amplitude^2 ~ -sigma/l1 > 0, so the Hopf
    % eigenpair's real part sigma has sign opposite to l1 (see side_sign)
    s = side_sign(sys, u_H, lambda_H, l1);
    opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true, ...
                        'FunctionTolerance', 1e-12, 'StepTolerance', 1e-14, 'MaxIterations', 200);

    % Seed a little off the Hopf, then Newton-refine to a true cycle. Near a weakly
    % nonlinear (near-Bautin) Hopf the amplitude grows slowly with the offset, so a
    % small offset gives a tiny orbit that the corrector collapses back onto the
    % trivial equilibrium (a shooting solution for any T). Grow the offset until the
    % landed orbit keeps most of the seed's amplitude, so the branch starts on a
    % genuine cycle rather than the degenerate equilibrium line.
    dl = abs(dlambda);
    for attempt = 1:8
        lam0 = lambda_H + s*dl;
        ue = fsolve(@(x) sys.F(x, lam0), u_H, optimset('Display', 'off'));

        % Poincare section at the equilibrium's I-level: the cycle surrounds the
        % equilibrium, so this stays crossed transversally as the cycle grows
        I_sec = ue(2);
        prob = shooting_system(sys, I_sec);

        % seed onto the real orbit by integration (l1 < 0 stable cycle -> forward;
        % l1 > 0 unstable cycle -> reverse time)
        [x0, T0] = seed_cycle(sys, ue, lam0, I_sec, -sign(l1));
        seed_amp = orbit_amp(sys, [x0; T0; lam0], n);
        land = @(z) deal([prob.residual(z); z(end) - lam0], ...
                         [prob.jacobian(z); zeros(1, n), 0, 1]);
        z0 = fsolve(land, [x0; T0; lam0], opts);
        if orbit_amp(sys, z0, n) > 0.5*seed_amp
            break
        end
        dl = 2*dl; % corrector fell back to the equilibrium: seed further out
    end

    % initial tangent = null direction of the shooting Jacobian, oriented to move
    % away from the Hopf (increasing |lambda - lambda_H|)
    t0 = null(prob.jacobian(z0));
    t0 = t0(:, 1)/norm(t0(:, 1));
    if sign(t0(end))*sign(lam0 - lambda_H) < 0
        t0 = -t0;
    end

    % amplitude floor: the cycle vanishes at the Hopf, where it degenerates to
    % the equilibrium (a trivial shooting solution for any T). Stopping just
    % above zero amplitude terminates the branch cleanly at the Hopf instead of
    % crossing that trivial line and retracing a spurious mirror.
    amp_floor = 0.05*orbit_amp(sys, z0, n);

    % also stop when lambda leaves its range, the period blows up, or shooting
    % loses accuracy. T -> inf is the homoclinic-approach signature where single
    % shooting gets fragile; the trivial Floquet multiplier (theoretically 1)
    % drifting from 1 is the built-in accuracy monitor and the hand-off cue to
    % the (future) global solver. monodromy reuses the cached integration.
    % follow the branch both ways from the seed: one way grows the cycle toward
    % the homoclinic, the other shrinks it back to the Hopf (small amplitude),
    % capturing any fold of cycles (LPC) where the stability flips in between.
    % stop() reuses the one monodromy it integrates for both the accuracy monitor
    % and (when seeding) the fold detector; prev_g carries the previous point's
    % (nontrivial multiplier - 1) so a sign change between points ends the arc.
    prev_g = []; % reset before fwd
    fwd = pseudo_arclength(prob, z0,  t0, step_size, max_steps, @stop);
    prev_g = []; % reset before bwd (fresh, independent arc)
    bwd = pseudo_arclength(prob, z0, -t0, step_size, max_steps, @stop);
    pts = [fliplr(bwd(:, 2:end)), fwd]; % drop the shared seed point in bwd

    out.z = pts; % raw continuation points [x0; T; lambda], e.g. to seed an LPC curve
    out.lambda = pts(n+2, :);
    out.T = pts(n+1, :);
    out = envelope_and_floquet(out, prob, sys, pts, n);

    % Nested so it can update prev_g between points. One monodromy integration
    % (reused from the cached shooting solve at z) serves both the accuracy/param/
    % period stops and, when seeding, the fold detector.
    function tf = stop(z)
        M = prob.monodromy(z);
        tf = z(end) < lambda_range(1) || z(end) > lambda_range(2) ...
           || z(n+1) > T_cap*T0 || min(abs(eig(M) - 1)) > 1e-2 ...
           || orbit_amp(sys, z, n) < amp_floor;
        if stop_on_fold && ~tf
            mu = eig(M);
            [~, it] = min(abs(mu - 1));           % trivial multiplier (nearest 1)
            nt = mu([1:it-1, it+1:end]);
            g = real(nt(1)) - 1;                  % nontrivial multiplier - 1
            if ~isempty(prev_g) && isfinite(g) && prev_g*g < 0
                tf = true;                        % crossed +1 between points: fold
            end
            prev_g = g;
        end
    end
end

% Sign of the lambda offset toward the cycle: the cycle lives where amplitude^2 ~
% -sigma/l1 > 0, i.e. where the Hopf eigenpair's real part sigma has sign opposite
% to l1. sigma is read from the critical complex eigenpair (the one nearest the
% imaginary axis), not trace(Fu): trace equals 2*sigma only in 2-D, and in 3-D+ the
% extra real eigenvalue corrupts the sign and picks the wrong side of the Hopf.
% Reads the eigenvalues on the corrected equilibrium so it is not biased by u_H.
function s = side_sign(sys, u_H, lambda_H, l1)
    h = 1e-4;
    u = fsolve(@(x) sys.F(x, lambda_H + h), u_H, optimset('Display', 'off'));
    ev = eig(sys.Fu(u, lambda_H + h));
    cpx = ev(abs(imag(ev)) > 1e-9);
    if isempty(cpx)
        sigma = trace(sys.Fu(u, lambda_H + h)); % fallback: no complex pair resolved
    else
        [~, i] = min(abs(real(cpx))); % the eigenpair nearest the imaginary axis
        sigma = real(cpx(i));
    end
    s = -sign(sigma)*sign(l1);
    if s == 0
        s = 1;
    end
end

% Integrate from just off the equilibrium ue onto the nearby limit cycle and
% read a section-crossing point and the period from the settled orbit. dir = +1
% integrates forward (a stable cycle attracts), dir = -1 reverse (an unstable
% cycle attracts in reverse time). Successive upward crossings of I = I_sec give
% the period; the last one is the most settled and seeds the shooting solve.
function [x0, T] = seed_cycle(sys, ue, lambda, I_sec, dir)
    o = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, ...
               'Events', @(t, y) settle_events(t, y, I_sec, ue));
    [~, ~, te, ye, ie] = ode45(@(t, y) dir*sys.F(y, lambda), [0 5000], ue + 1e-3, o);
    up = ie == 1; % section crossings; event 2 is the blow-up guard
    te = te(up);
    ye = ye(up, :);
    if numel(te) < 3
        error('trace_limit_cycle:seed', 'could not settle onto a cycle to seed from');
    end
    x0 = ye(end-1, :).';
    T = te(end) - te(end-1);
end

% Peak-to-peak amplitude of the I-coordinate over one period of the orbit at z.
function a = orbit_amp(sys, z, n)
    o = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    [~, Y] = ode45(@(t, y) sys.F(y, z(n+2)), [0 z(n+1)], z(1:n), o);
    a = max(Y(:, 2)) - min(Y(:, 2));
end

% Upward section crossings (event 1, non-terminal) to read the period, plus a
% blow-up guard (event 2, terminal): reverse-time integration of a near-
% homoclinic cycle runs away, so stop the instant the orbit leaves a large
% neighbourhood of the equilibrium instead of letting ode45 grind to the horizon.
function [v, isterm, dir] = settle_events(~, y, I_sec, ue)
    v      = [y(2) - I_sec; norm(y - ue) - 1e3];
    isterm = [0;            1];
    dir    = [1;            1];
end

% Per stored cycle: integrate one period for the state max/min envelope, and
% read Floquet multipliers from the cached monodromy. In 2-D one multiplier is
% trivially ~1 (flow direction); the other governs stability (|mu| < 1 stable).
function out = envelope_and_floquet(out, prob, sys, pts, n)
    N = size(pts, 2);
    out.umax = zeros(n, N);
    out.umin = zeros(n, N);
    out.multipliers = zeros(n, N);
    out.stable = false(1, N);
    o = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
    for j = 1:N
        [~, Y] = ode45(@(t, y) sys.F(y, pts(n+2, j)), [0 pts(n+1, j)], pts(1:n, j), o);
        out.umax(:, j) = max(Y, [], 1)';
        out.umin(:, j) = min(Y, [], 1)';

        mu = eig(prob.monodromy(pts(:, j)));
        out.multipliers(:, j) = mu;
        [~, it] = min(abs(mu - 1)); % the trivial multiplier
        nontriv = mu([1:it-1, it+1:end]);
        out.stable(j) = all(abs(nontriv) < 1);
    end
end
