% Continue the limit cycle born at every Hopf point found by trace_branches.
% branches_out is the trace_branches struct; its hopfs / hopf_l1 fields supply
% the seeds. Returns a cell array of trace_limit_cycle structs, one per Hopf,
% ready for plot_limit_cycles. Any trailing arguments (offset, step, max steps,
% range) are forwarded to trace_limit_cycle. Two-parameter analogue: this is to
% trace_limit_cycle what trace_branches is to a single branch.
function lcs = trace_limit_cycles(sys, branches_out, varargin)
    lcs = {};
    for i = 1:numel(branches_out.hopfs)
        H = branches_out.hopfs{i};
        l1 = branches_out.hopf_l1{i};
        for k = 1:size(H, 2)
            % one cycle branch per Hopf, seeded from that Hopf's point + l1
            lcs{end+1} = trace_limit_cycle(sys, H(1:end-1, k), H(end, k), l1(k), varargin{:}); %#ok<AGROW>
        end
    end
end
