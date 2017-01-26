#' Effect size values for Connectivity Map build 02 drugs.
#'
#' Unbiased effect sizes values for all 1309 drugs in the Connectivity Map build 02.
#'
#' @usage data(cmap_es)
#' @return A matrix where columns correspond to drugs and rows to gene symbols.
"cmap_es"


#' Variance values for Connectivity Map build 02 drugs.
#'
#' Variances of unbiased effect sizes values for all 1309 drugs in the
#' Connectivity Map build 02.
#'
#' @usage data(cmap_var)
#' @return A matrix where columns correspond to drugs and rows to gene symbols.
"cmap_var"


#' Effect size values for LINCS l1000 signatures.
#'
#' Unbiased effect sizes values for all 230829 LINCS l1000 signatures.
#'
#' @usage data(l1000_es)
#' @return A matrix where columns correspond to perturbagens and rows to gene symbols.
"l1000_es"


#' HGNC symbols used for NNet predictions.
#'
#' Order is as required for input and produced by output of net1/net2 predictions.
#'
#' @usage data(genes)
#' @return A character vector of 11525 HGNC symbols.
"genes"


#' Neural network model 1 for treatment combinations.
#'
#' Contains weight matrices and bias vectors needed to make predictions.
#'
#' @usage #NA
#' @return List with matrices W1/W2 and vectors b1/b2.
"net1"


#' Neural network model 2 for treatment combinations.
#'
#' Contains weight matrices and bias vectors needed to make predictions.
#'
#' @usage #NA
#' @return List with matrices W1/W2 and vectors b1/b2.
"net2"


#' XGBoost model for treatment combinations.
#'
#' Model stacks predictions from net1 and net2 with effect size values from
#' cmap_es and variance values from cmap_var.
#'
#' @usage #NA
#' @return Object of class xgb.Booster
"xgb_mod"
