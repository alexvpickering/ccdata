#' Effect size values for Connectivity Map build 02 drugs.
#'
#' Effect size values for all 1309 drugs in the Connectivity Map build 02.
#'
#' @docType data
#' @keywords datasets
#' @name cmap_es
#' @usage data(cmap_es)
#' @format A matrix where columns correspond to drugs and rows to gene symbols.
NULL



#' Top tables for Connectivity Map build 02 drugs.
#'
#' Adjusted p-values from call to limma \code{topTables} for all 1309 Connectivity
#' Map build 02 drugs. Also includes unbiased effect sizes ('dprimes') from
#' results of call to metaMA \code{effectsize}.
#'
#' @docType data
#' @keywords datasets
#' @name cmap_tables
#' @usage data(cmap_tables)
#' @format A named list containing data.frames with dprime and adjusted p-values
#'    for each drug in CMap build 02.
NULL



#' XGBoost model for treatment combinations.
#'
#' XGBoost model used to predict gene expression changes resulting from a
#' combination treatment.
#'
#' Predictions for combinations are made using data from the individual
#' treatments. Model was trained using 'combo_train' data and is used to predict
#' Connectivity Map build 02 combinations using 'cmap_tables'.
#'
#' @docType data
#' @keywords datasets
#' @name combo_model
#' @usage data(combo_model)
NULL
