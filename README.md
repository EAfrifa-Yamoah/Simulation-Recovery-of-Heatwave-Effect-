# Declined, displaced, or just harder to catch?

A simulation framework for attributing marine heatwave impacts on fished populations from catch per unit effort (CPUE). This repository reproduces every number, table and figure in the manuscript from R. The manuscript itself is not distributed here.

## Authors

- Eben Afrifa-Yamoah, ORCID [0000-0003-1741-9249](https://orcid.org/0000-0003-1741-9249)
- Yaw Kwaafo Awuah-Mensah, ORCID [0009-0002-7893-0948](https://orcid.org/0009-0002-7893-0948)

## What is here

```
R/
  01_operating_model.R   configuration, 2^4 scenario design, events, temperature, niche,
                         movement with a dynamic refuge beyond the frame, state transition,
                         non negative removal based fishing mortality, catchability
  02_observation.R       tow level delta lognormal catch (encounter and variance scale with
                         effort), effort allocation with an event driven contraction that
                         enters as an interaction, survey with frequency and stratification
  03_estimators.R        nominal, positive only, delta two way (ridge on separated cells,
                         convergence and imputation accounting, fished cells only variant),
                         delta with spatial interaction, survey (design and model based),
                         centroids, heatwave slope with Newey West standard errors
  04_simulate.R          ONE sequential loop for the factual world; counterfactual with the
                         same effort and noise; two estimands (total and in frame abundance)
  05_identifiability.R   invariance illustration (uniform direction singular) and the
                         informative check: spatial catchability is identified
  06_experiments.R       2^4 factorial + warm edge cells, frame geometry, survey design
                         curve (coverage x frequency x design, CI coverage and power),
                         matched calibration, temporal recovery with nested GCV,
                         global sensitivity, classifier
scripts/
  run_all.R              every experiment, stage cached, parallel over replicates
  make_tables.R          all tables -> results/tables/*.md and *.csv
  make_figures.R         all figures -> figures/*.pdf and 300 dpi *.png (base graphics)
tests/testthat/          40 checks, including the ones that would have caught the
                         bugs removed in version 2 (see Version 2 changes below)
```

## Reproduce

Requires R >= 4.3 with `mgcv` and `testthat`, plus `pdftoppm` (poppler). No other language is used.

```bash
make test        # unit and property tests (about 2 minutes)
make quick       # the whole pipeline at small sizes (about 2 minutes)
make run         # full run: ~12,000 replicates, ~1 hour on two cores
make tables      # results/tables/
make figures     # figures/
```

`make run` caches each stage in `results/parts/`; delete a stage's RDS to rerun it. `MHW_CORES=n` sets the number of forked workers (default: all but one). Seeds are a deterministic function of the replicate index, so scenarios share random numbers and the results are bit reproducible for a given R version.

## Design, in one paragraph

The operating model couples abundance, availability and catchability through a single thermal driver. The fished domain is a 12 x 12 lattice with an open poleward boundary; a refuge one band beyond the edge experiences the same seasonal cycle, warming trend and event anomaly as the frame. Total abundance follows a density dependent transition with an occupancy weighted thermal dose and a removal based fishing mortality. A fleet allocates fixed monthly effort by a multinomial on its lagged belief about cell catch rates, with an event driven contraction that sharpens its targeting. Catch is observed tow by tow through a mean coherent delta lognormal. An optional survey samples the frame by a fixed design with protocol fixed catchability. Every replicate is run twice under common random numbers, with and without the event anomaly, which defines two estimands: the event effect on total abundance and on in frame abundance. Six estimators are benchmarked against both over a 2^4 factorial (decline, redistribution, uniform catchability, fleet contraction) plus warm edge catchability cells.

## Version 2 changes (October 2026)

Fixes to the simulation that followed a code level review: observations are now drawn from the same abundance the truth is defined on; fishing mortality cannot be negative; the fleet contraction, which previously cancelled in the softmax and did nothing, now contracts the footprint; the refuge beyond the frame warms with the frame; the declared estimand is the one used in every comparison; the survey is benchmarked against in frame abundance, for which it is unbiased, and against total abundance, for which it is not when biomass crosses the frame boundary; survey frequency and stratification are part of the design curve; confidence interval coverage and power are reported; the distributed lag smoothing parameter is chosen inside each training block; the integrated estimator, which shared a month effect across streams and absorbed the fishery's catchability change, was removed. The Python port and the Node document pipeline were removed; everything is R.
