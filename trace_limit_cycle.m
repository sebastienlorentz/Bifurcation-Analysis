% Continue the limit cycle born at a Hopf point in one parameter, by orthogonal
% collocation. Seeds the first orbit by integrating onto it (forward for a stable
% cycle, reverse for an unstable one), samples it onto the collocation mesh, then
% follows it with pseudo_arclength on the one-parameter orbit_collocation_system.
% l1 is the first Lyapunov coefficient at the Hopf (its sign sets which side the
% cycle lives and the cycle's stability); dlambda is the initial parameter offset.
% Returns the max/min envelope of each state, the period, and Floquet multipliers
% along the branch, ready to overlay on the diagram from trace_branches. Floquet
% multipliers come from monodromy_integrate, one variational integration per stored
% point, since collocation carries no monodromy of its own.
function out = trace_limit_cycle(sys, u_H, lambda_H, l1, dlambda, step_size, max_steps, lambda_range, T_cap, stop_on_fold, Ntst, Ncol)
    n = sys.n;
    u_H = u_H(:);
    if nargin < 11 || isempty(Ntst); Ntst = 20; end
    if nargin < 12 || isempty(Ncol); Ncol = 4;  end
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
    % accuracy gate: how far the trivial Floquet multiplier may drift from 1 before the
    % branch is stopped (see the stop function below)
    drift_tol = 1e-2;

    % side where the cycle lives: amplitude^2 ~ -sigma/l1 > 0, so the Hopf
    % eigenpair's real part sigma has sign opposite to l1 (see side_sign)
    s = side_sign(sys, u_H, lambda_H, l1);
    opts = optimoptions('fsolve', 'Display', 'off', 'SpecifyObjectiveGradient', true, ...
                        'FunctionTolerance', 1e-12, 'StepTolerance', 1e-14, 'MaxIterations', 200);

    % the one-parameter collocation problem: z = [X(:); T; lambda]
    fld = struct('n', n, 'nparam', 1, 'F', @(x, l) sys.F(x, l), ...
                 'Fu', @(x, l) sys.Fu(x, l), 'Fp', @(x, l) sys.Flambda(x, l));
    csys = orbit_collocation_system(fld, Ntst, Ncol);
    Np = csys.Np;
    nX = n*Np;
    mono = @(z) monodromy_of(sys, z, n, nX);

    % Seed a little off the Hopf, then Newton-refine to a true cycle. Near a weakly
    % nonlinear (near-Bautin) Hopf the amplitude grows slowly with the offset, so a
    % small offset gives a tiny orbit that the corrector collapses back onto the
    % trivial equilibrium (a solution of the BVP for any T). Grow the offset until the
    % landed orbit keeps most of the seed's amplitude, so the branch starts on a
    % genuine cycle rather than the degenerate equilibrium line.
    dl = abs(dlambda);
    for attempt = 1:8
        lam0 = lambda_H + s*dl;
        ue = fsolve(@(x) sys.F(x, lam0), u_H, optimset('Display', 'off'));

        % a coordinate section at the equilibrium's I-level still gives the seeding
        % integration a period to read off; the collocation problem itself uses the
        % integral phase condition, so no section is carried into the continuation
        I_sec = ue(2);

        % seed onto the real orbit by integration (l1 < 0 stable cycle -> forward;
        % l1 > 0 unstable cycle -> reverse time), then sample it onto the mesh
        [x0, T0] = seed_cycle(sys, ue, lam0, I_sec, -sign(l1));
        Xref = csys.seed_profile(x0, T0, lam0);
        seed_amp = profile_amp(Xref);
        prob = csys.make_problem(Xref);

        % land on a true cycle with lambda held: the orbit residual is one row short
        % of square, so pinning lambda closes it
        pin = [zeros(1, nX + 1), 1];
        land = @(z) deal([prob.residual(z); z(end) - lam0], [prob.jacobian(z); pin]);
        z0 = fsolve(land, [Xref(:); T0; lam0], opts);
        if profile_amp(reshape(z0(1:nX), n, Np)) > 0.5*seed_amp
            break
        end
        dl = 2*dl; % corrector fell back to the equilibrium: seed further out
    end

    % initial tangent = null direction of the collocation Jacobian, oriented to move
    % away from the Hopf (increasing |lambda - lambda_H|)
    t0 = null(prob.jacobian(z0));
    t0 = t0(:, 1)/norm(t0(:, 1));
    if sign(t0(end))*sign(lam0 - lambda_H) < 0
        t0 = -t0;
    end

    % amplitude floor: the cycle vanishes at the Hopf, where it degenerates to
    % the equilibrium (a trivial solution for any T). Stopping just above zero
    % amplitude terminates the branch cleanly at the Hopf instead of crossing that
    % trivial line and retracing a spurious mirror.
    amp_floor = 0.05*profile_amp(reshape(z0(1:nX), n, Np));

    % follow the branch both ways from the seed: one way grows the cycle toward
    % the homoclinic, the other shrinks it back to the Hopf (small amplitude),
    % capturing any fold of cycles (LPC) where the stability flips in between.
    % prev_g carries the previous point's (nontrivial multiplier - 1) so a sign
    % change between points ends the arc when seeding.
    prev_g = []; % reset before fwd
    fwd = pseudo_arclength(prob, z0,  t0, step_size, max_steps, @stop);
    prev_g = []; % reset before bwd (fresh, independent arc)
    bwd = pseudo_arclength(prob, z0, -t0, step_size, max_steps, @stop);
    pts = [fliplr(bwd(:, 2:end)), fwd]; % drop the shared seed point in bwd

    out.z = pts;            % raw continuation points [X(:); T; lambda]
    out.X = pts(1:nX, :);   % full orbit profile per point
    out.x0 = pts(1:n, :);   % first mesh node, the section point an LPC seed wants
    out.lambda = pts(nX+2, :);
    out.T = pts(nX+1, :);
    out = envelope_and_floquet(out, mono, pts, n, nX, Np);

    % Nested so it can update prev_g between points. One monodromy integration per
    % point serves both the accuracy monitor and (when seeding) the fold detector.
    %
    % The trivial multiplier drifting from 1 is that accuracy monitor, and it is
    % still needed on collocation. Collocation removes shooting's corrector
    % fragility, but not two other limits: a fixed mesh cannot represent an orbit
    % whose period has grown far beyond the seed's (the profile starts oscillating
    % and dips unphysically negative), and the monodromy is still an ODE integration
    % that loses accuracy over a long period. Without this stop the branch runs on to
    % T_cap producing profiles and multipliers that are both meaningless.
    function tf = stop(z)
        M = mono(z);
        tf = z(end) < lambda_range(1) || z(end) > lambda_range(2) ...
           || z(nX+1) > T_cap*T0 || min(abs(eig(M) - 1)) > drift_tol ...
           || profile_amp(reshape(z(1:nX), n, Np)) < amp_floor;
        if stop_on_fold && ~tf
            nt = nontrivial_multipliers(eig(M));
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
% the period; the last one is the most settled and seeds the collocation solve.
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

% Peak-to-peak amplitude of the I-coordinate, read straight off the mesh profile.
function a = profile_amp(X)
    a = max(X(2, :)) - min(X(2, :));
end

% Monodromy at a collocation point, by one variational integration from the first
% mesh node. The collocation problem carries no monodromy of its own.
function M = monodromy_of(sys, z, n, nX)
    [~, M] = monodromy_integrate(sys.F, sys.Fu, n, z(1:n), z(nX+1), z(nX+2));
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

% Per stored cycle: the state max/min envelope straight off the mesh profile (no
% integration needed, the whole orbit is already solved for), and Floquet
% multipliers from one variational integration. In 2-D one multiplier is trivially
% ~1 (flow direction); the other governs stability (|mu| < 1 stable).
function out = envelope_and_floquet(out, mono, pts, n, nX, Np)
    N = size(pts, 2);
    out.umax = zeros(n, N);
    out.umin = zeros(n, N);
    out.multipliers = zeros(n, N);
    out.stable = false(1, N);
    for j = 1:N
        X = reshape(pts(1:nX, j), n, Np);
        out.umax(:, j) = max(X, [], 2);
        out.umin(:, j) = min(X, [], 2);

        mu = eig(mono(pts(:, j)));
        out.multipliers(:, j) = mu;
        out.stable(j) = all(abs(nontrivial_multipliers(mu)) < 1);
    end
end
