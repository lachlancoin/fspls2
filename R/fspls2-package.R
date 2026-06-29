#' @importFrom grDevices colorRampPalette
#' @importFrom graphics title
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics title
#' @importFrom methods as slotNames
#' @importFrom stats ar coef cor deviance dnorm ecdf glm knots logLik median pchisq plogis qchisq quantile update var
#' @importFrom utils capture.output head read.csv read.delim
#' @importFrom Matrix sparseMatrix
#' @importFrom glmnet glmnet
#' @importFrom tidyr pivot_wider pivot_longer unite separate
#' @importFrom binom binom.confint
#' @importFrom ggrepel geom_text_repel
#' @importFrom jsonlite fromJSON toJSON
#' @importFrom confintr ci_cor ci_var
#' @importFrom pROC ci roc
#' @importFrom MASS polr
#' @importFrom DBI dbConnect dbDisconnect dbWriteTable dbReadTable dbListTables dbGetQuery
#' @importFrom RSQLite SQLite
#' @importFrom R6 R6Class

#' @importFrom ggplot2 ggplot aes geom_point geom_line theme_bw labs geom_ribbon aes_string geom_hline facet_grid facet_wrap scale_color_manual geom_abline ggtitle geom_text ggsave guides element_text scale_shape_manual scale_x_continuous scale_y_continuous scale_y_log10 sec_axis theme
#' @importFrom RColorBrewer brewer.pal
#' @importFrom Matrix nnzero Matrix
NULL

utils::globalVariables(c(
  "angle","pval",
  "spec",
  "Column",
  "counts",
  "counts ",
  "data",
  "dataset",
  "drug",
  "factor",
  "func",
  "gene",
  "knots",
  "label",
  "lens",
  "lineage",
  "model",
  "name",
  "ni1",
  "pheno",
  "pheno_code",
  "Row",
  "sample",
  "sens",
  "sign",
  "subpheno",
  "typ1",
  "type",
  "value",
  "X33",
  "X50",
  " X66",
  "xx",
  "y",
  "y_1",
  "CV", "beam", "beta_new2", "beta_scale", "cohort_measure_pheno_trained", "comb", "cv",
  "cv_full", "dims", "experiment_id" ,"measure", "mid", "nrep", "phens", "roc" ,"submeasure"
))
