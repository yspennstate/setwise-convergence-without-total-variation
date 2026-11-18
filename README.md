# Setwise convergence of shifted laws without total variation convergence

Thorisson asked whether setwise convergence of the shifted laws of a discrete-time process to the law of a stationary
process forces convergence in total variation, that is, agreement on the tail sigma-field. It does not.
If the law of the process has a density with respect to a stationary mixing law, the shifted laws converge setwise,
and their total variation distance to the limit is governed by the conditional expectation of the density given the
tail sigma-field. Orbits of an invertible mixing transformation started from a non-invariant density are therefore
counterexamples, and zero-entropy mixing processes give counterexamples with two states. Processes with independent
coordinates give none. For Thorisson's companion question on couplings, a single coupling under which the shifted
processes asymptotically agree on every event is sufficient but not necessary: for the horocycle flow, no coupling
brings the processes together in probability.

- `paper/` the LaTeX source and PDF.
- `code/` `sanity.py` checks several finite statements of the paper in exact rational arithmetic; its output is
  `sanity_output.txt`. No proof in the paper depends on it.

Requirements: Python 3.8 or later, standard library only. Run `cd code && python sanity.py`.
