# Bifurcation continuation for ODE systems

A small MATLAB library for continuation and bifurcation analysis of
n-dimensional ODE systems, built on an adaptive pseudo-arclength engine. It traces equilibrium
branches in one parameter and fold, Hopf and transcritical curves in two,
classifies stability, and finds codimension-2 points (Bogdanov-Takens,
generalized Hopf). It also continues the limit cycles born at Hopf points (by
single shooting or orthogonal collocation), with Floquet stability and
fold-of-cycles detection, and continues the fold-of-cycles (LPC) curve itself in
two parameters from the generalized-Hopf point.

The example driver scripts apply it to an SIS epidemic model with saturated
treatment,

```
S' = lambda - beta*S*I - mu*S + gamma*I
I' = beta*S*I - (mu + gamma + f(I))*I,   f(I) = 1/(1 + I^2)
```

with state `u = [S; I]` and `R0 = beta*lambda / (mu*(mu + gamma + f(0)))`, but
the core library isn't specific to it; pass any symbolic `F(u, p)`.

## Requirements

MATLAB with the Symbolic Math Toolbox and
the Optimization Toolbox.

## Quick start

Put all the scripts on the MATLAB path and run one of the drivers:

- `example_one_parameter_R0`: equilibrium branches vs `R0`, with stability,
  folds, Hopfs and branch points marked, plus the limit cycles born at each Hopf
  continued by single shooting (with the period vs `R0` companion figure).
- `example_sirs_3d`: the same one-parameter analysis for a 3-D SIRS model,
  exercising the dimension-general Hopf test, limit cycles and phase portrait
  (each drawn for the state components named in an index vector).
- `example_bifurcation_curves`: two-parameter `(R0, gamma)` bifurcation diagram
  (fold, Hopf and transcritical curves plus BT / generalized-Hopf points).
- `example_lpc_curve`: the same two-parameter diagram plus the fold-of-cycles
  (LPC) curve, seeded from a fold of cycles on a subcritical slice and continued
  back up to the generalized-Hopf point, overlaid in green.
- `example_phase_portrait`: phase portrait and trajectories at fixed parameters
  (projected onto any 2 states for higher-dimensional systems).
- `example_time_series`: the states against time from given initial conditions;
  the period (and its divergence near a global bifurcation) is visible directly.

## Library overview

Engine
- `pseudo_arclength`: adaptive pseudo-arclength corrector/predictor; traces a
  curve of `residual(z) = 0` given a point and a tangent.
- `follow_curve`: shared tail of the two-parameter tracers; follows a landed
  curve both ways within a parameter box and assembles the two arms.

Symbolic setup (turns a symbolic `F` into the handles the solvers need)
- `continuation_system`: one continuation parameter; F, Jacobian and the
  derivatives used for branch switching and Lyapunov coefficients.
- `bifurcation_derivatives`: shared derivative handles for the two-parameter
  fold and Hopf systems.
- `lyapunov_system` / `lyapunov_coefficient`: first Lyapunov coefficient at a
  Hopf point (sign gives super/subcriticality).

Augmented continuation systems (two parameters)
- `fold_continuation_system`
- `hopf_continuation_system`
- `transcritical_continuation_system`

Limit cycles (single shooting)
- `shooting_system`: the periodicity + phase-condition problem `z = [x0; T; p]`,
  with the monodromy and parameter sensitivity from one variational ODE.
- `monodromy_integrate`: the shared augmented variational integration returning
  `phi_T`, the monodromy `M(T)` and the parameter sensitivity.

Limit cycles (orthogonal collocation)
- `orbit_collocation_system`: the periodic-orbit BVP on a mesh of Ntst intervals
  with degree-Ncol Lagrange polynomials collocated at Gauss points; solves the
  whole orbit at once, staying well conditioned through a fold of cycles.

Fold of cycles (LPC) in two parameters
- `lpc_continuation_system` / `trace_lpc_curve`: the LPC curve by single
  shooting, using the `det(M) = 1` fold condition.
- `lpc_collocation_system` / `trace_lpc_collocation`: the same curve by
  collocation with an explicit null vector, robust through the near-homoclinic
  tail and the faster production path.
- `find_lpc_seed`: sweeps parameter slices to find a fold of cycles (a +1
  Floquet crossing) to seed the LPC continuation.

Tracers
- `trace_branches`: all equilibrium branches in one parameter, switching at
  branch points; flags folds, Hopfs and stability.
- `trace_fold_curve` / `trace_hopf_curve` / `trace_transcritical_curve`: a
  single two-parameter curve from a starting point.
- `trace_bifurcation_curves`: finds and traces every fold/Hopf/transcritical
  curve in a parameter box, plus Bogdanov-Takens and generalized-Hopf points.
- `trace_limit_cycle` / `trace_limit_cycles`: the limit cycle from a single Hopf
  point in one parameter (envelope, period, Floquet multipliers, stability,
  stopping at the homoclinic approach), and the wrapper that continues the
  cycles from every Hopf in a `trace_branches` result.

Plotting
- `plot_branches`: one-parameter equilibrium diagram.
- `plot_limit_cycles`: overlays the cycles' max/min envelopes on the equilibrium
  diagram, plus a shared period vs parameter figure.
- `plot_bifurcation_curves`: two-parameter bifurcation diagram.
- `plot_time_series`: integrate the system and plot the states in `idx` against
  time, one figure per initial condition.
- `phase_portrait`: speed heatmap, direction field and trajectories at fixed
  parameters, projected onto the 2 states in `idx`.
