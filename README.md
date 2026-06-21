# Bifurcation continuation for a 2D SIS model

A small MATLAB library for continuation and bifurcation analysis of 2D ODE
systems, built on an adaptive pseudo-arclength engine. It traces equilibrium
branches in one parameter and fold, Hopf and transcritical curves in two,
classifies stability, and finds codimension-2 points (Bogdanov-Takens,
generalized Hopf).

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

- `branch_switching_R0`: equilibrium branches vs `R0`, with stability, folds,
  Hopfs and branch points marked.
- `branch_switching_gamma`: equilibrium branches vs the treatment rate `gamma`.
- `example_bifurcation_curves`: two-parameter `(R0, gamma)` bifurcation diagram
  (fold, Hopf and transcritical curves plus BT / generalized-Hopf points).
- `phase_portrait_combined`: phase portrait and trajectories at fixed parameters.

## Library overview

Engine
- `pseudo_arclength`: adaptive pseudo-arclength corrector/predictor; traces a
  curve of `residual(z) = 0` given a point and a tangent.

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

Tracers
- `trace_branches`: all equilibrium branches in one parameter, switching at
  branch points; flags folds, Hopfs and stability.
- `trace_fold_curve` / `trace_hopf_curve` / `trace_transcritical_curve`: a
  single two-parameter curve from a starting point.
- `trace_bifurcation_curves`: finds and traces every fold/Hopf/transcritical
  curve in a parameter box, plus Bogdanov-Takens and generalized-Hopf points.

Plotting
- `plot_branches`: one-parameter equilibrium diagram.
- `plot_bifurcation_curves`: two-parameter bifurcation diagram.
