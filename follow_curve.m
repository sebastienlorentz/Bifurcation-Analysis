% Follow a two-parameter bifurcation curve both ways from a landed start point z0,
% stopping when either parameter leaves p_box. This is the shared tail of the
% trace_*_curve functions: initial tangent from the augmented Jacobian's null
% space, the box-exit stop, the forward/backward arcs, and their assembly into one
% point set (the shared seed dropped once from the backward arm). Each caller keeps
% only its own landing (z0 construction) and output unpacking. Also returns why
% each arm ended and its length, for callers that report the terminus.
function [pts, fwd_reason, bwd_reason, fwd_n, bwd_n] = follow_curve(prob, z0, p_box, step_size, max_steps, extra_stop)
    t0 = null(prob.jacobian(z0));
    t0 = t0(:, 1) / norm(t0(:, 1));

    box = @(z) z(end-1) < p_box(1, 1) || z(end-1) > p_box(1, 2) || z(end)   < p_box(2, 1) || z(end)   > p_box(2, 2);
    % optional extra terminator (e.g. cycle amplitude vanishing at a Bautin point)
    if nargin >= 6 && ~isempty(extra_stop)
        stop = @(z) box(z) || extra_stop(z);
    else
        stop = box;
    end

    [fwd, ~, ~, fwd_reason] = pseudo_arclength(prob, z0,  t0, step_size, max_steps, stop);
    [bwd, ~, ~, bwd_reason] = pseudo_arclength(prob, z0, -t0, step_size, max_steps, stop);
    fwd_n = size(fwd, 2);
    bwd_n = size(bwd, 2);
    pts = [fliplr(bwd(:, 2:end)), fwd]; % drop the shared seed point in bwd
end
