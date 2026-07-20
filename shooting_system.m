% Build a single-shooting limit-cycle problem for pseudo_arclength. The unknown
% is z = [x0; T; lambda] (length n+2); the residual enforces periodicity
% phi_T(x0) - x0 = 0 (n rows) plus a coordinate Poincare section x0(2) - I_sec
% (1 row). sys supplies F, Fu, Flambda (from continuation_system). Residual and
% jacobian share one augmented ODE integration, memoized on z so Newton's
% separate residual/jacobian calls integrate only once per iteration.
function prob = shooting_system(sys, I_sec)
    n = sys.n;

    zc  = [];   % cached z for the integration below
    xT  = [];   % phi_T(x0)
    Mon = [];   % monodromy M(T) = d phi_T / d x0
    pT  = [];   % parameter sensitivity p(T) = d phi_T / d lambda

    prob.residual  = @residual;
    prob.jacobian  = @jacobian;
    prob.monodromy = @monodromy; % M(T) at a point, for Floquet multipliers

    % integrate x' = F, M' = Fu*M, p' = Fu*p + Flambda over [0, T] once per z
    function integrate(z)
        if isequal(z, zc)
            return
        end
        x0 = z(1:n);
        T = z(n+1);
        lambda = z(n+2);
        [xT, Mon, pT] = monodromy_integrate(@(x, l) sys.F(x, l), @(x, l) sys.Fu(x, l), ...
                                            @(x, l) sys.Flambda(x, l), n, 1, x0, T, lambda);
        zc = z;
    end

    function R = residual(z)
        integrate(z);
        x0 = z(1:n);
        R = [xT - x0; x0(2) - I_sec];
    end

    % columns [x0 (n) | T | lambda]; the section row is the constant [0 1 | 0 | 0]
    function Jz = jacobian(z)
        integrate(z);
        lambda = z(n+2);
        dphi_dT = sys.F(xT, lambda); % endpoint vector field
        e2 = zeros(1, n);
        e2(2) = 1;
        Jz = [Mon - eye(n), dphi_dT, pT;
              e2,           0,       0];
    end

    function M = monodromy(z)
        integrate(z);
        M = Mon;
    end
end
