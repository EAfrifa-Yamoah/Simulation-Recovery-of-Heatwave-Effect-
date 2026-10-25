# Everything runs in R. See README.
.PHONY: all run tables figures paper doc tex test quick clean

all: run tables figures paper

run:
	Rscript scripts/run_all.R

quick:
	MHW_QUICK=1 Rscript scripts/run_all.R

tables:
	Rscript scripts/make_tables.R

figures:
	Rscript scripts/make_figures.R

paper: tables figures
	Rscript scripts/make_tex_tables.R
	Rscript scripts/assemble_manuscript.R
	cd manuscript && latexmk -pdf -interaction=nonstopmode -outdir=build sim-rec_main.tex
	cd manuscript && latexmk -pdf -interaction=nonstopmode -outdir=build sim-rec_supp.tex

doc:
	Rscript scripts/build_manuscript.R

tex: doc
	Rscript scripts/build_latex.R

test:
	NOT_CRAN=true Rscript -e 'testthat::test_dir("tests/testthat", reporter = "summary")'

clean:
	rm -rf results/parts results/*.rds results/*.csv results/tables results/manuscript_build.md results/run.log
