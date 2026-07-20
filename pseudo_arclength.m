% Trace a 1-D solution curve of prob.residual(z)=0 by adaptive pseudo-arclength continuation.
% Fourth output `reason` reports why the arc ended: 'stop' (left the box),
% 'nonfinite' (orbit couldn't be integrated, the near-homoclinic terminus),
% 'stuck' (corrector couldn't converge even at ds_min), or 'max_steps'.
function [points, tangents, arclen, reason] = pseudo_arclength(prob, z0, t0, step_size, max_steps, stop)
    m = numel(z0);
    points = zeros(m, max_steps);
    tangents = zeros(m, max_steps);
    arclen = zeros(1, max_steps); % cumulative arclength at each stored point
    points(:, 1) = z0;
    tangents(:, 1) = t0;
    z = z0;
    t = t0;
    npts = 1; % number of points stored so far
    reason = 'max_steps'; % overwritten below if the arc ends early

    newton_iters = 50;
    tol = 1e-8;

    % adaptive step control: step_size is the starting/target arclength step.
    % It grows when the corrector solves easily and shrinks (with retry) when a
    % step fails or rounds a corner too fast, staying within [ds_min, ds_max].
    ds = step_size;
    ds_min = step_size*1e-4;
    ds_max = step_size*100; % cap kept modest so a step can't jump over a fold
    newton_target = 4; % Newton iterations we aim for per accepted step
    min_cos_turn = cos(0.3); % reject a step that turns the tangent > 0.3 radians

    for iter = 1:max_steps
        if stop(z)
            reason = 'stop';
            break
        end

        % attempt one step at the current ds; on a bad step shrink and retry
        % until it is accepted or ds bottoms out at ds_min
        accepted = false;
        t_new = t;
        while true
            z_predict = z+ds*t;

            % corrector: solve residual = 0 with the arclength constraint that
            % keeps us ds along the tangent from the previous point
            G = @(w) [prob.residual(w); t'*(w-z)-ds];
            w = z_predict;
            converged = false;
            nonfinite = false;
            for N = 1:newton_iters
                g = G(w);
                if ~all(isfinite(g))
                    nonfinite = true; % integration failed (e.g. near-homoclinic timeout)
                    break
                end
                ng = norm(g);
                if ng < tol
                    converged = true;
                    break
                end

                % backtracking line search along the Newton direction: an
                % ill-conditioned Jacobian (e.g. the doubly-singular monodromy at
                % a fold of cycles) makes a full step overshoot and diverge, so
                % accept only the fraction that actually decreases the residual
                dw = -[prob.jacobian(w); t'] \ g;
                a = 1;
                while a > 1e-6
                    gn = G(w + a*dw);
                    if all(isfinite(gn)) && norm(gn) < ng
                        break
                    end
                    a = a/2;
                end
                if a <= 1e-6
                    break % no downhill step at this ds; leave unconverged and shrink
                end
                w = w + a*dw;
            end

            % a non-finite residual means the orbit itself couldn't be integrated;
            % shrinking ds won't recover it, so end this arc now (keeping the
            % points already traced, including any fold found earlier on it)
            if nonfinite
                reason = 'nonfinite';
                break
            end

            % a merely unconverged step (ill-conditioned corner near a fold, e.g.
            % where a steep nonlinearity keeps Fu near-singular over a range) is
            % recoverable: fall through to the shrink-and-retry path below
            bad = ~converged || norm(w - z_predict) > 2*ds;

            % tangent for the next step (null direction of the Jacobian); also
            % used here to reject steps that round a corner too sharply
            if ~bad
                t_new = [prob.jacobian(w); t'] \ [zeros(m-1, 1); 1];
                t_new = t_new / norm(t_new);
                if t'*t_new < min_cos_turn
                    bad = true;
                end
            end

            if ~bad
                accepted = true;
                break
            end
            if ds <= ds_min*(1+1e-9)
                break % cannot shrink any further; give up on this arc
            end
            ds = max(ds/2, ds_min);
        end

        if ~accepted
            reason = 'stuck';
            break
        end

        z = w;
        t = t_new; % last row of the bordered system pins the sign
        npts = iter+1;
        points(:, npts) = z;
        tangents(:, npts) = t;
        arclen(npts) = arclen(npts-1)+ds;

        % size the next step from how hard the corrector worked: grow when
        % Newton was easy (N small), shrink when N large, but not by more
        % than 2x either way, and never past ds_max
        growth = min(max(newton_target/N, 0.5), 2);
        ds = min(ds*growth, ds_max);
    end

    points = points(:, 1:npts);
    tangents = tangents(:, 1:npts);
    arclen = arclen(1:npts);
end
