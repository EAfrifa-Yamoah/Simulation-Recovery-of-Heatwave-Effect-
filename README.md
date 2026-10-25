# Declined, displaced, or just harder to catch?

A simulation framework for attributing marine heatwave impacts on fished populations from catch per unit effort (CPUE). This repository reproduces every number, table and figure in the manuscript, and builds the manuscript itself, from R.

## Authors

- Yaw Kwaafo Awuah-Mensah, Edith Cowan University, ORCID [0009-0002-7893-0948](https://orcid.org/0009-0002-7893-0948)
- Eben Afrifa-Yamoah, Edith Cowan University, ORCID [0000-0003-1741-9249](https://orcid.org/0000-0003-1741-9249)

## Citation

If you use this code, please cite it as follows. The DOI below is the Zenodo concept DOI and resolves to the latest version; each release also has its own version DOI on Zenodo.

Awuah-Mensah, Y. K., & Afrifa-Yamoah, E. (2026). Code: Declined, displaced, or just harder to catch? A simulation framework for attributing marine heatwave impacts on fished populations from catch per unit effort [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23293721

```bibtex
@software{awuahmensah2026code,
  author    = {Awuah-Mensah, Yaw Kwaafo and Afrifa-Yamoah, Eben},
  title     = {Code: Declined, displaced, or just harder to catch? A simulation framework for attributing marine heatwave impacts on fished populations from catch per unit effort},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.23293721},
  url       = {https://doi.org/10.5281/zenodo.23293721}
}
```

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
  make_tex_tables.R      LaTeX table bodies from results/tables/*.csv -> results/tables/tex/
  assemble_manuscript.R  templates + table bodies -> manuscript/sim-rec_main.tex and sim-rec_supp.tex
  build_manuscript.R     manuscript.md -> manuscript_sim-rec_revised.docx via pandoc,
                         with every quoted number evaluated from results/
  build_latex.R          the same resolved source -> manuscript_sim-rec_revised.tex
                         (numbered equations, no indent after equations, the
                         manuscript's own figure and table numbers) and, when
                         latexmk is present, the PDF
tests/testthat/          62 checks, including the ones that would have caught the
                         bugs removed in version 2 (see RESPONSE_TO_REVIEW.md)
manuscript/
  main_template.tex      main text source (Fish and Fisheries); every quoted number and
                         table body is inserted from results/ by assemble_manuscript.R
  supp_template.tex      Supporting Information source
  references.tex         reference list shared by both
  sim-rec_main.tex, sim-rec_supp.tex   the assembled files to compile (latexmk -pdf)
manuscript.md            source of the earlier Word draft (make doc, make tex); superseded by
                         manuscript/main_template.tex and supp_template.tex
```

## Reproduce

Requires R >= 4.3 with `mgcv` and `testthat`, plus `pandoc` and `pdftoppm` (poppler). No other language is used.

```bash
make test        # unit and property tests (about 2 minutes)
make quick       # the whole pipeline at small sizes (about 2 minutes)
make run         # full run: ~12,000 replicates, ~1 hour on two cores
make tables      # results/tables/
make figures     # figures/
make paper       # manuscript/build/sim-rec_main.pdf and sim-rec_supp.pdf
make doc         # manuscript_sim-rec_revised.docx (earlier Word draft)
make tex         # manuscript_sim-rec_revised.tex (+ PDF via latexmk if installed)
```

`make run` caches each stage in `results/parts/`; delete a stage's RDS to rerun it. `MHW_CORES=n` sets the number of forked workers (default: all but one). Seeds are a deterministic function of the replicate index, so scenarios share random numbers and the results are bit reproducible for a given R version.

## Design, in one paragraph

The operating model couples abundance, availability and catchability through a single thermal driver. The fished domain is a 12 x 12 lattice with an open poleward boundary; a refuge one band beyond the edge experiences the same seasonal cycle, warming trend and event anomaly as the frame. Total abundance follows a density dependent transition with an occupancy weighted thermal dose and a removal based fishing mortality. A fleet allocates fixed monthly effort by a multinomial on its lagged belief about cell catch rates, with an event driven contraction that sharpens its targeting. Catch is observed tow by tow through a mean coherent delta lognormal. An optional survey samples the frame by a fixed design with protocol fixed catchability. Every replicate is run twice under common random numbers, with and without the event anomaly, which defines two estimands: the event effect on total abundance and on in frame abundance. Six estimators are benchmarked against both over a 2^4 factorial (decline, redistribution, uniform catchability, fleet contraction) plus warm edge catchability cells.

## Version 2 changes (October 2026)

Fixes to the simulation that followed a code level review: observations are now drawn from the same abundance the truth is defined on; fishing mortality cannot be negative; the fleet contraction, which previously cancelled in the softmax and did nothing, now contracts the footprint; the refuge beyond the frame warms with the frame; the declared estimand is the one used in every comparison; the survey is benchmarked against in frame abundance, for which it is unbiased, and against total abundance, for which it is not when biomass crosses the frame boundary; survey frequency and stratification are part of the design curve; confidence interval coverage and power are reported; the distributed lag smoothing parameter is chosen inside each training block; the integrated estimator, which shared a month effect across streams and absorbed the fishery's catchability change, was removed. The Python port and the Node document pipeline were removed; everything is R.

## Fish and Fisheries manuscript (October 2026)

The manuscript is now maintained as LaTeX in `manuscript/main_template.tex` and `manuscript/supp_template.tex`, assembled by `scripts/assemble_manuscript.R` (`make paper`). Changes to code and results made with it:

- Global sensitivity: the fleet contraction parameter `beta_J` was sampled over [-1.6, 0], a range that excluded the model's default (+0.8) and, below about -0.5, reversed the fleet's targeting. It is now sampled over [0, 1.6] (`SENS_RANGES` in `R/06_experiments.R`), and a test checks that every range contains its default. The sensitivity stage was rerun: contraction strength now has almost no effect on the nominal bias (SRC 0.004), which is governed mainly by uniform catchability loss (-0.58) and the number of events (-0.31).
- Figures: Figure 4 (sampling distributions) read non-existent columns and drew no violins; fixed. Figures use `cairo_pdf` so degree signs render; titles were removed from inside panels; non-estimator elements use neutral colours so no estimator colour carries a second meaning; the combined sensitivity and identifiability figure is split into `fig11_identifiability` and `fig12_sensitivity`.
- Centre of gravity: shifts are event contrasts of the standardised northing and are now reported in km (173 km per unit; 3.46 latitude bands of 50 km), not as bands.

## Changes after the code review of 23 October 2026

- With redistribution off, the fixed distribution is now the mean suitability over the series without heatwaves, identical in the factual and counterfactual runs (`distribution_shares()`, `run_world()` arguments `Tmat_base`, `phi_out_base`). Previously each run used the mean over its own series, so the factual run's shares included future events and the paired contrast was non-zero before the first heatwave. A test now checks that the contrast is exactly zero before the first anomaly. All stages were rerun; changes to reported means are within Monte Carlo error.
- Month 1 is the initial state and receives no transition or fishing mortality; the loop and the Methods now say so.
- `scripts/run_all.R` keys every cached stage to an md5 of the run sizes, the configuration and the R source; a stage with a different key is recomputed.
- The delta index is described as effort-unadjusted; an effort-adjusted comparison is in the Supporting Information.
- The sensitivity, thermal-stress-beyond-the-frame and effort-adjustment checks referred to in the Supporting Information were run with `scripts/review_checks.R`.
