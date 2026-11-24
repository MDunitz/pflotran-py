# Forecasting test: fit early, predict late

Each experiment here was split in time rather than by experiment. The two
salinity inhibition parameters -- the chloride concentration at which the
term centres and the width of its transition -- were fitted using only the
first 2 sampling rounds of that experiment, and then asked to
reproduce the rounds that came after. Nothing else was adjusted.

This is the question a modeller faces mid-experiment: a few weeks of data
exist, the rest of the incubation does not yet. Can the early rounds pin the
parameters well enough to say where the bottles end up?

Day zero is excluded from every window. It is a baseline reading taken
before anything has been produced: the model sits at its trace floor by
construction while the instrument reports its background, and the ratio
between those measures nothing about the model. In an objective built on
log ratios that one point dominates everything else -- including it gave
fit scores near 4.7, a nominal miss by a factor of fifty thousand, and made
the held-out window score better than the fitted one.

Each experiment therefore has four production rounds. The first 2
are fitted on and the rest are predicted.

## Results

Score is the root mean squared error of log10(modelled/measured). Lower is
better; a value of one is a typical miss by a factor of ten.

| Experiment | Chloride threshold | Transition width | Fitted rounds | Held-out rounds |
|---|---|---|---|---|
| Exp003 | 2.0 mol/L | 1.0 decades | 0.89 | 0.84 |
| Exp004 | 2.0 mol/L | 1.0 decades | 0.82 | 1.47 |

## The selection hit the edge of the grid

Both experiments chose the corner of the parameter grid: the highest
chloride threshold and the widest transition available, which together are
the weakest inhibition it can express. A search that stops at its own
boundary has not found a minimum; it has run out of room. The true optimum
on early data lies outside the grid, at weaker inhibition still.

That is itself the finding. Fitted on the whole timecourse the same
objective prefers a threshold of 0.75 mol/L, which is strong inhibition.
Fitted on the first two rounds it prefers 2.0 mol/L or beyond, which is
almost none. The early rounds do not merely constrain the parameters
loosely -- they point in the opposite direction.

The reason is visible in the figures. At days 10 and roughly 55 the salted
bottles have not yet separated from the controls by much, so a model with
little salt inhibition matches them well. The collapse develops later, and
by then the parameters have been chosen.

## What each folder holds

- `methane_forecast.png` -- the timecourse, with the fitting window shaded.
  Circles inside the shading were used to choose the parameters; squares
  outside it were not.
- `fit_grid.csv` -- every parameter combination and its score on the fitting
  window, so the shape of the objective is inspectable rather than summarised
  by its minimum alone.
- `per_timepoint.csv` -- modelled and measured methane at every batch and
  sampling day, labelled by which window it fell in.

## Caveat

These experiments were examined during the work that produced this model, so
this is a pre-registered protocol rather than a blind forecast. The objective
and the parameter grid are fixed in advance and the selection is arithmetic,
but knowledge of how the incubations end cannot be unlearned.

Regenerate with `python -m pflotran_py.comparison.forecast`.
