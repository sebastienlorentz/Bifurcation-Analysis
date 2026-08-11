% Integrate the augmented variational system over [0, T] in one ode45 call:
%   X' = rhs(X, p)             X(0) = X0     state      (dim)
%   M' = jac(X, p) * M         M(0) = I      monodromy  (dim x dim)
% and return the endpoint values X(T) and M(T) = dX(T)/dX0. rhs and jac are handles of
% (X, p). Sole caller is trace_limit_cycle, which uses MT for the Floquet multipliers:
% collocation solves the orbit but carries no monodromy of its own.
function [XT, MT] = monodromy_integrate(rhs, jac, dim, X0, T, p)
    t_start = tic;
    % tol below the pseudo_arclength corrector target (1e-8): a residual that depends on
    % this integration can't be driven under its noise floor, so the floor must sit
    % beneath the tolerance the corrector has to reach
    o = odeset('RelTol', 1e-9, 'AbsTol', 1e-11, ...
               'Events', @(t, y) blowup_event(y, dim), ...
               'OutputFcn', @(t, y, flag) toc(t_start) > 20); % wall-clock cap (s)
    nM = dim*dim;
    y0 = [X0; reshape(eye(dim), [], 1)];
    [tt, Y] = ode45(@(t, y) rhs_aug(y, rhs, jac, dim, nM, p), [0 T], y0, o);

    % A run that stopped before T was truncated by the blow-up guard or the wall-clock
    % cap: near a homoclinic the period diverges and the orbit crawls past the saddle, so
    % a single ode45 call can take unboundedly many steps. Return NaN so the caller's
    % residual is non-finite and pseudo_arclength ends the arc at once, instead of Newton
    % re-integrating a hopeless orbit.
    if abs(tt(end) - T) > 1e-6*max(1, abs(T))
        XT = nan(dim, 1);
        MT = nan(dim, dim);
        return
    end
    yT = Y(end, :)';
    XT = yT(1:dim);
    MT = reshape(yT(dim+1:dim+nM), dim, dim);
end

function dy = rhs_aug(y, rhs, jac, dim, nM, p)
    X = y(1:dim);
    M = reshape(y(dim+1:dim+nM), dim, dim);
    dy = [rhs(X, p);
          reshape(jac(X, p)*M, [], 1)];
end

% Terminal guard: if the state part diverges the orbit has run away (a near-homoclinic
% blow-up), so stop rather than let ode45 grind to tiny steps. Truncation before T is
% turned into a NaN endpoint above, ending the arc.
function [v, isterm, dir] = blowup_event(y, dim)
    v = norm(y(1:dim)) - 1e6;
    isterm = 1;
    dir = 1;
end
