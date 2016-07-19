#' Effect size values for Connectivity Map build 02 drugs.
#'
#' Effect size values for all 1309 drugs in the Connectivity Map build 02.
#'
#' @usage data(cmap_es)
#' @return A matrix where columns correspond to drugs and rows to gene symbols.
"cmap_es"



#' Top tables for Connectivity Map build 02 drugs.
#'
#' Adjusted p-values from call to limma \code{topTables} for all 1309
#' Connectivity Map build 02 drugs. Also includes unbiased effect sizes
#' (\code{dprimes}) from results of call to metaMA \code{effectsize}.
#'
#' @usage data(cmap_tables)
#' @return A named list containing data.frames with dprime and adjusted p-values
#'    for each drug in the Connectivity Map build 02.
"cmap_tables"



#' XGBoost model for treatment combinations.
#'
#' XGBoost model used to predict gene expression changes resulting from a
#' combination treatment.
#'
#' Predictions for combinations are made using data from the individual
#' treatments. Model was trained using \code{combo_train} data and is used to
#' predict Connectivity Map build 02 combinations using \code{cmap_tables}.
#'
#' @usage data(combo_model)
#' @format An object of class xgb.Booster.
#' @return xgb.Booster object
"combo_model"
