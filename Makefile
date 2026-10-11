# Everything runs in R. See README.
.PHONY: all run tables figures test quick clean

all: run tables figures

run:
	Rscript scripts/run_all.R

quick:
	MHW_QUICK=1 Rscript scripts/run_all.R

tables:
	Rscript scripts/make_tables.R

figures:
	Rscript scripts/make_figures.R

test:
	NOT_CRAN=true Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'

clean:
	rm -rf results/parts results/*.rds results/*.csv results/tables results/manuscript_build.md results/run.log
