# =============================================================================
# scripts/assemble_manuscript.R
# Fill the manuscript templates with table bodies (results/tables/tex/) and with
# every number quoted in the text, computed here from results/, then write
#   manuscript/sim-rec_main.tex   (main text, Fish and Fisheries)
#   manuscript/sim-rec_supp.tex   (Supporting Information)
# Usage:
#   Rscript scripts/make_tex_tables.R && Rscript scripts/assemble_manuscript.R
#   cd manuscript && latexmk -pdf sim-rec_main.tex && latexmk -pdf sim-rec_supp.tex
# =============================================================================
Sys.setlocale("LC_CTYPE", "C.UTF-8")
for (f in list.files("R", full.names = TRUE)) source(f)
tex_dir <- "results/tables/tex"

## ---- numbers quoted in the text ---------------------------------------------
fs  <- read.csv("results/factorial_summary.csv")
fac <- readRDS("results/factorial.rds")
cs  <- read.csv("results/coverage_summary.csv")
fr  <- readRDS("results/frame.rds")
tp  <- readRDS("results/temporal.rds")
cl  <- readRDS("results/classifier.rds")
mm  <- read.csv("results/matched.csv")
sens <- read.csv("results/sensitivity_design.csv")
cfg <- readRDS("results/config.rds"); grid <- make_grid(cfg)
KM  <- cfg$cell_km * sd(grid$lat)

f3 <- function(x) { v <- formatC(round(x, 3), digits = 3, format = "f"); if (abs(round(x, 3)) < 5e-4) v <- "0.000"
  if (substr(v, 1, 1) == "-") paste0("$-$", substring(v, 2)) else v }
f2 <- function(x) { v <- formatC(round(x, 2), digits = 2, format = "f"); if (abs(round(x, 2)) < 5e-3) v <- "0.00"
  if (substr(v, 1, 1) == "-") paste0("$-$", substring(v, 2)) else v }
pc <- function(x) sprintf("%.0f", 100 * x)
sc <- function(s, col) fs[[col]][fs$scenario == s]
mg <- function(k, d, col) mm[[col]][abs(mm$kappa - k) < 1e-9 & abs(mm$delta_dep - d) < 1e-9]
cov <- function(every, design, coverage, col, cv = cfg$cv_surv)
  cs[[col]][cs$every == every & cs$design == design & abs(cs$coverage - coverage) < 1e-9 & abs(cs$cv_surv - cv) < 1e-9][1]
frm <- function(frame, scen, col) mean(fr[[col]][fr$frame == frame & fr$scenario == scen], na.rm = TRUE)
cg  <- function(scen, col) KM * mean(fac[[col]][fac$scenario == scen], na.rm = TRUE)
srcv <- function(resp, par) { z <- src(sens, resp); z$src[z$parameter == par] }
srcr2 <- function(resp) src(sens, resp)$r2[1]
cm <- as.matrix(cl$confusion)
ab <- cm["abundance", ]

redis <- fs$scenario[grepl("^S.1..$", fs$scenario) | fs$scenario %in% c("SPQT0", "SPQT1")]
surv_err_in <- max(abs(fs$beta_survey[fs$scenario %in% redis] - fs$truth_in[fs$scenario %in% redis]))
dev1001 <- max(abs(c(sc("S1001", "beta_nominal"), sc("S1001", "beta_delta"), sc("S1001", "beta_deltaint"),
                     sc("S1001", "beta_survey")) - sc("S1001", "truth_total")))

vals <- list(
  MCSE = signif(max(fs[, grep("^mcse_beta", names(fs))], na.rm = TRUE), 2),
  MG1T = f3(mg(0, 0.07, "beta_truth")),   MG1N = f3(mg(0, 0.07, "beta_nominal")),   MG1S = f3(mg(0, 0.07, "beta_survey")),
  MG2T = f3(mg(0.06, 0.035, "beta_truth")), MG2N = f3(mg(0.06, 0.035, "beta_nominal")), MG2S = f3(mg(0.06, 0.035, "beta_survey")),
  MG3T = f3(mg(0.12, 0, "beta_truth")),   MG3N = f3(mg(0.12, 0, "beta_nominal")),   MG3S = f3(mg(0.12, 0, "beta_survey")),
  T_S1001 = f3(sc("S1001", "truth_total")), DEV_S1001 = f3(ceiling(dev1001 * 1000) / 1000),
  NOM_S0011 = f3(sc("S0011", "beta_nominal")), DEL_S0011 = f3(sc("S0011", "beta_delta")), SUR_S0011 = f3(sc("S0011", "beta_survey")),
  NOM_S1011 = f3(sc("S1011", "beta_nominal")), T_S1011 = f3(sc("S1011", "truth_total")), SUR_S1011 = f3(sc("S1011", "beta_survey")),
  TI_S0101 = f3(sc("S0101", "truth_in")), SUR_S0101 = f3(sc("S0101", "beta_survey")), SURVERR = f3(ceiling(surv_err_in * 1000) / 1000),
  DI_S0101 = f3(sc("S0101", "beta_deltaint")), DI_S1101 = f3(sc("S1101", "beta_deltaint")), TI_S1101 = f3(sc("S1101", "truth_in")),
  NOM_S0101 = f3(sc("S0101", "beta_nominal")), DEL_S0101 = f3(sc("S0101", "beta_delta")),
  DI_S0111 = f3(sc("S0111", "beta_deltaint")), TI_S0111 = f3(sc("S0111", "truth_in")),
  NOM_S0001 = f3(sc("S0001", "beta_nominal")),
  CG_TRUE = sprintf("%.0f", cg("S0101", "dG_true")), CG_FISH = sprintf("%.0f", cg("S0101", "dG_fishery")),
  CG_SURV = sprintf("%.0f", cg("S0101", "dG_survey")), CG_GRID = sprintf("%.0f", cg("S0101", "dG_grid")),
  CG_FISH_WE = sprintf("%.0f", cg("SPQ1", "dG_fishery")), CG_GRID_WE = sprintf("%.0f", cg("SPQ1", "dG_grid")),
  CG_SURV_WE = sprintf("%.1f", cg("SPQ1", "dG_survey")),
  CG_FISH_WER = sprintf("%.0f", cg("SPQT1", "dG_fishery")), CG_SURV_WER = sprintf("%.0f", cg("SPQT1", "dG_survey")),
  NOM_SPQ1 = f3(sc("SPQ1", "beta_nominal")), T_SPQ1 = f3(sc("SPQ1", "truth_total")),
  RMSE_NOM = f3(cov(1, "srs", 0.10, "rmse_nominal")), RMSE_DEL = f3(cov(1, "srs", 0.10, "rmse_delta")),
  COV_NOM = pc(cov(1, "srs", 0.10, "cover_nominal")), COV_DEL = pc(cov(1, "srs", 0.10, "cover_delta")),
  RMSE_M7 = f3(cov(1, "srs", 0.05, "rmse")), RMSE_M14 = f3(cov(1, "srs", 0.10, "rmse")),
  COV_M14 = pc(cov(1, "srs", 0.10, "cover_survey")), EMPTY_M3 = sprintf("%.2f", cov(1, "srs", 0.02, "empty_months")),
  RMSE_Q29 = f3(cov(3, "srs", 0.20, "rmse")), COV_Q29 = pc(cov(3, "srs", 0.20, "cover_survey")),
  COV_H14 = pc(cov(6, "srs", 0.10, "cover_survey")),
  RMSE_M14CV = f3(cov(1, "srs", 0.10, "rmse", cv = 0.90)),
  POW_M14 = pc(cov(1, "srs", 0.10, "power_survey")), POW_M58 = pc(cov(1, "srs", 0.40, "power_survey")),
  POW_NOM = pc(cov(1, "srs", 0.10, "power_nominal")), POW_DEL = pc(cov(1, "srs", 0.10, "power_delta")),
  FR_DYN_TI = f3(frm("open_dynamic", "S0101", "beta_cf_in")), FR_DYN_S = f3(frm("open_dynamic", "S0101", "beta_survey")),
  FR_STAT_TI = f3(frm("open_static", "S0101", "beta_cf_in")), FR_STAT_S = f3(frm("open_static", "S0101", "beta_survey")),
  FR_CL_N = f3(frm("closed", "S0101", "beta_nominal")), FR_DYN_N = f3(frm("open_dynamic", "S0101", "beta_nominal")),
  FR_STAT_N = f3(frm("open_static", "S0101", "beta_nominal")),
  TP_R2A = f2(mean(tp$summary_obs$r2_acute)), TP_R2D = f2(mean(tp$summary_obs$r2_dl)),
  TP_CTA = f2(mean(tp$summary_obs$cor_obs_truth_ac)), TP_CTD = f2(mean(tp$summary_obs$cor_obs_truth_dl)),
  TP_TRUA = f2(mean(tp$summary_tru$cor_acute)), TP_TRUD = f2(mean(tp$summary_tru$cor_dl)),
  TP_OOSA = f2(mean(tp$summary_obs$oos_acute)), TP_OOSD = f2(mean(tp$summary_obs$oos_dl)),
  SR_GAP_DEP = f2(srcv("gap_nominal", "delta_dep")), SR_GAP_NEV = f2(srcv("gap_nominal", "n_events")),
  SR_GAP_KAP = f2(srcv("gap_nominal", "kappa")), SR_GAP_STH = f2(srcv("gap_nominal", "sigma_th")),
  SR_GAP_RHO = f2(srcv("gap_nominal", "rho_move")), SR_GAP_BJ = f3(srcv("gap_nominal", "beta_J")),
  SR_GAP_LR = f2(srcv("gap_nominal", "belief_lr")),
  SR_RMSE_COV = f2(srcv("rmse_survey", "surv_cover")), SR_RMSE_AS = f2(srcv("rmse_survey", "a_surv")),
  SR_RMSE_KAP = f2(srcv("rmse_survey", "kappa")),
  SR_BIAS_RHO = f2(srcv("bias_survey_total", "rho_move")), SR_BIAS_NEV = f2(srcv("bias_survey_total", "n_events")),
  SR_BIAS_WT = f2(srcv("bias_survey_total", "warm_trend")), SR_BIAS_STH = f2(srcv("bias_survey_total", "sigma_th")),
  SR_R2_GAP = pc(srcr2("gap_nominal")), SR_R2_BIAS = pc(srcr2("bias_survey_total")),
  KEY_ACC = pc(cl$accuracy), KEY_IND = pc(cl$indeterminate),
  KEY_WE = pc(cm["warm edge catchability", "warm edge catchability"] / sum(cm["warm edge catchability", ])),
  KEY_RE = pc(cm["redistribution", "redistribution"] / sum(cm["redistribution", ])),
  KEY_AB_OK = ab["abundance"], KEY_AB_N = sum(ab), KEY_AB_NONE = ab["none"], KEY_AB_IND = ab["indeterminate"]
)
fr_txt <- strsplit(readLines(file.path(tex_dir, "design_curve_fishery_rmse.txt")), " ")[[1]]
r2s <- strsplit(readLines(file.path(tex_dir, "sensitivity_r2.txt")), ",\\s*")[[1]]
vals$FISHRMSE <- paste(fr_txt[1], "and", fr_txt[2])
vals$SENSR2 <- paste0(paste(r2s[1:3], collapse = ", "), " and ", r2s[4])

bib <- readLines("manuscript/references.tex", warn = FALSE, encoding = "UTF-8")

fill <- function(template, out_file) {
  tpl <- readLines(template, warn = FALSE, encoding = "UTF-8")
  out <- character(0)
  for (ln in tpl) {
    m <- regmatches(ln, regexec("^%%TAB:([A-Za-z0-9_]+)%%$", ln))[[1]]
    if (length(m)) {
      f <- file.path(tex_dir, paste0(m[2], ".tex"))
      if (!file.exists(f)) stop("missing table body: ", f)
      out <- c(out, readLines(f, warn = FALSE, encoding = "UTF-8"))
    } else if (ln == "%%BIBLIOGRAPHY%%") out <- c(out, "%%BIBLIOGRAPHY%%") else out <- c(out, ln)
  }
  # keep only the references this document cites, in the order of references.tex
  cited <- unique(unlist(lapply(regmatches(out, gregexpr("\\\\cite[a-z]*(\\[[^]]*\\])*\\{([^}]*)\\}", out)),
                                function(x) unlist(strsplit(gsub(".*\\{([^}]*)\\}$", "\\1", x), "[, ]+")))))
  cited <- cited[nzchar(cited)]
  item_key <- sub("^\\\\bibitem\\[[^]]*\\]\\{([^}]*)\\}.*$", "\\1", bib)
  item_key[!grepl("^\\\\bibitem", bib)] <- NA
  grp <- cumsum(!is.na(item_key)); grp[grp == 0] <- NA
  keys <- item_key[!is.na(item_key)]
  keep <- is.na(grp) | keys[grp] %in% cited | grepl("^\\\\end\\{thebibliography\\}", bib)
  missing <- setdiff(cited, item_key[!is.na(item_key)])
  if (length(missing)) stop("cited but not in references.tex: ", paste(missing, collapse = ", "))
  message(basename(template), ": ", length(cited), " references cited of ", sum(!is.na(item_key)))
  i <- which(out == "%%BIBLIOGRAPHY%%"); out <- c(out[seq_len(i - 1)], bib[keep], out[-seq_len(i)])
  for (nm in names(vals)) out <- gsub(paste0("%%", nm, "%%"), as.character(vals[[nm]]), out, fixed = TRUE)
  out[1:9] <- sub("TEMPLATE.*", paste("Generated by scripts/assemble_manuscript.R from", basename(template)), out[1:9])
  left <- unique(unlist(regmatches(out, gregexpr("%%[A-Z][A-Za-z0-9_]*%%", out))))
  if (length(left)) stop("unresolved placeholders in ", template, ": ", paste(left, collapse = " "))
  writeLines(out, out_file, useBytes = TRUE)
  cat("wrote", out_file, "\n")
}
fill("manuscript/main_template.tex", "manuscript/sim-rec_main.tex")
fill("manuscript/supp_template.tex", "manuscript/sim-rec_supp.tex")
