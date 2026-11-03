# Bifurcation continuation for ODE systems

A small MATLAB library for continuation and bifurcation analysis of n-dimensional ODE
systems, built on an adaptive pseudo-arclength engine. It traces equilibrium branches in
one parameter and fold, Hopf and transcritical curves in two, classifies stability, and
locates codimension-2 points (Bogdanov–Takens, generalized Hopf). On top of that it
continues the limit cycles born at Hopf points by orthogonal collocation, with Floquet
stability and fold-of-cycles detection, continues the fold-of-cycles (LPC) curve itself
in two parameters, and continues the homoclinic orbits that terminate the cycle family.

The example drivers apply it to an SIS epidemic model with saturated treatment,

```
S' = mu - beta*S*I - mu*S + gamma*I
I' = beta*S*I - (mu + gamma + f(I))*I,   f(I) = 1/(1 + N0^2*I^2)
```

written in normalized form: the states are fractions of the disease-free population
`N0 = lambda/mu`, so the disease-free state sits at `S = 1`, the transmission is
`beta = R0*(mu + gamma + f(0))`, and the treatment saturation carries the scale `N0^2`.
This is a pure change of state variables, so it leaves every bifurcation unchanged. The
core library is not specific to the model; pass any symbolic `F(u, p)`.

## Requirements

MATLAB with the Symbolic Math Toolbox and the Optimization Toolbox.

## Quick start

Put all the scripts on the MATLAB path and run one of the drivers:

- `example_one_parameter_R0`: equilibrium branches vs `R0`, with stability, folds, Hopfs
  and branch points marked, plus the limit cycles born at each Hopf continued by
  orthogonal collocation, and the period vs `R0` companion figure.
- `example_bifurcation_curves`: two-parameter `(R0, gamma)` diagram of the equilibrium
  bifurcations (fold, Hopf and transcritical curves plus the BT and generalized-Hopf
  points).
- `example_lpc_curve`: the same diagram plus the fold-of-cycles (LPC) curve in green,
  seeded from a fold of cycles on a subcritical slice and continued back up to the
  generalized-Hopf point.
- `example_homoclinic_curve`: the same diagram plus the homoclinic curve in purple,
  seeded from a moderate-period orbit on the LPC arc. All five curve families at once.
- `example_homoclinic_orbit`: one homoclinic orbit in the phase plane, with the saddle,
  the truncation ball its two ends are pinned to, and the equilibrium it encircles.
- `example_truncation_failure`: what goes wrong when the truncation distance is held
  constant rather than pinned to the loop amplitude, run before and after the fix on the
  same seed.
- `example_sirs_3d`: the one-parameter analysis for a 3-D SIRS model, exercising the
  dimension-general Hopf test, limit cycles and phase portrait (each drawn for the state
  components named in an index vector).
- `example_phase_portrait`: phase portrait and trajectories at fixed parameters
  (projected onto any 2 states for higher-dimensional systems).
- `example_time_series`: the states against time from given initial conditions.

The drivers that produce the report's figures write vector PDFs to `docs/figures/`
through `export_figure`.

## Library overview

Engine
- `pseudo_arclength`: adaptive pseudo-arclength predictor/corrector; traces a curve of
  `residual(z) = 0` from a point and a tangent.
- `follow_curve`: shared tail of the two-parameter tracers; follows a landed curve both
  ways within a parameter box and assembles the two arms.

Symbolic setup (turns a symbolic `F` into the handles the solvers need)
- `continuation_system`: one continuation parameter; F, Jacobian and the derivatives
  used for branch switching and Lyapunov coefficients.
- `bifurcation_derivatives`: shared derivative handles for the two-parameter systems.
- `lyapunov_system` / `lyapunov_coefficient`: first Lyapunov coefficient at a Hopf point
  (its sign gives super/subcriticality).

Augmented continuation systems (two parameters)
- `fold_continuation_system`
- `hopf_continuation_system`
- `transcritical_continuation_system`

Limit cycles (orthogonal collocation)
- `orbit_collocation_system`: the periodic-orbit BVP on a mesh of Ntst intervals with
  degree-Ncol Lagrange polynomials collocated at Gauss points; solves the whole orbit at
  once, staying well conditioned through a fold of cycles. Used for both the one- and
  two-parameter cycle continuations.
- `monodromy_integrate`: the augmented variational integration returning the monodromy
  `M(T)`. Collocation carries no monodromy of its own, so Floquet multipliers still come
  from here, one integration per stored point.

Fold of cycles (LPC) in two parameters
- `lpc_collocation_system` / `trace_lpc_collocation`: the LPC curve by collocation with
  an explicit null vector, robust through the near-homoclinic tail.
- `find_lpc_seed`: sweeps parameter slices for a fold of cycles (a +1 Floquet crossing)
  to seed the LPC continuation.

Homoclinic orbits in two parameters
- `homoclinic_collocation_system`: the truncated homoclinic BVP on the same collocation
  mesh, with the saddle's eigenvectors carried as unknowns so the Jacobian stays
  analytic, and the truncation distance pinned to a fixed fraction of the loop
  amplitude. See `docs/homoclinic_design_note.md` for the alternatives rejected.
- `trace_homoclinic_curve`: cuts a large-period cycle at its saddle passage, lands it on
  the BVP with one parameter pinned, and follows the curve both ways.

Tracers
- `trace_branches`: all equilibrium branches in one parameter, switching at branch
  points; flags folds, Hopfs and stability.
- `trace_fold_curve` / `trace_hopf_curve` / `trace_transcritical_curve`: a single
  two-parameter curve from a starting point.
- `trace_bifurcation_curves`: finds and traces every fold, Hopf and transcritical curve
  in a parameter box, plus the Bogdanov–Takens and generalized-Hopf points.
- `trace_limit_cycle` / `trace_limit_cycles`: the limit cycle from a single Hopf point in
  one parameter by collocation (envelope, period, Floquet multipliers, stability,
  stopping at the homoclinic approach), and the wrapper that continues the cycles from
  every Hopf in a `trace_branches` result. A branch stops once the trivial Floquet
  multiplier drifts from 1, which bounds it to where the fixed mesh and the variational
  integration are both still trustworthy.

Plotting
- `plot_branches`: one-parameter equilibrium diagram.
- `plot_limit_cycles`: overlays the cycles' max/min envelopes on the equilibrium
  diagram, plus a shared period vs parameter figure.
- `plot_bifurcation_curves`: two-parameter bifurcation diagram.
- `plot_time_series`: integrate the system and plot the states in `idx` against time,
  one figure per initial condition.
- `phase_portrait`: speed heatmap, direction field and trajectories at fixed parameters,
  projected onto the 2 states in `idx`.
- `export_figure`: write the current figure to `docs/figures/` as a vector PDF.

Shared helpers
- `bilin`: contracts a 3-D derivative tensor with two vectors, the bilinear form the
  branch-point tangents and the Lyapunov coefficient both need.
- `nontrivial_multipliers`: the Floquet multipliers with the trivial one at +1 removed.
- `split_stable`: masks a data series into its stable (solid) and unstable (dashed)
  parts for plotting.
- `legend_once`: gives a group of line handles a single legend entry.
