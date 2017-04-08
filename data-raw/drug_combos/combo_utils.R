# Setup drug combination training data
#
# User is asked to select drug1, drug2, and combo contrasts for each GSE.
# Selections are then used create training data.frame.
#
# @param diff_exprs Result from call to \code{\link[crossmeta]{diff_expr}}.
# @param prev_selections Used to re-use previous selections.
#
# @return data.frame

setup_combo_data <- function(diff_exprs, prev_selections=NULL) {

    #add moderated effect sizes
    diff_exprs <- add_es(diff_exprs, cols = "dprime")

    #add co-expression data
    #diff_exprs <- add_coxpres(diff_exprs)

    #get selections
    selections <- select_combo_data(diff_exprs, prev_selections)

    #combine
    data <- combine_combo_data(diff_exprs, selections)

    #return result
    combo_data <- list(data=data, selections=selections)
    return(combo_data)
}


add_coxpres <- function(diff_exprs) {

    hs <- readRDS("~/Documents/Batcave/GEO/ccdata/data-raw/coxpresdb/hs_coxpres.rds")
    mm <- readRDS("~/Documents/Batcave/GEO/ccdata/data-raw/coxpresdb/mm_coxpres.rds")

    diff_exprs <- lapply(diff_exprs, function(anal) {

        #determine species
        sp <- levels(pData(anal$eset)$organism_ch1)
        if (sp == "Homo sapiens") {
            coxpres <- hs
        } else {
            coxpres <- mm
        }

        #merge coxpres and anal entrez ids
        eids <- row.names(anal$top_tables[[1]])
        dt <- data.table(gene1 = eids)
        coxpres <- merge(dt, coxpres, by = "gene1")

        #remove coxpres rows where gene2 not in anal entrez ids
        coxpres <- coxpres[coxpres$gene2 %in% eids, ]

        #coxpress gene1 must have >= 10 entries
        three <- table(coxpres$gene1) >= 10
        coxpres <- coxpres[coxpres$gene1 %in% names(three)[three]]

        #remove anal entrez ids not in coexpress gene1
        #order top_tables the same
        eids <- eids[eids %in% coxpres$gene1]
        top_tables <- lapply(anal$top_tables, function(x) {
            x[eids, ]
        })

        #get first 10 entries for each gene1 in coxpres
        coxpres <- coxpres[, .SD[1:10], by = gene1]

        #cast coxpres
        coxpres$rank <- 1:10
        coxpres <- dcast(coxpres, gene1 ~ rank, value.var = c("mrank", "gene2"))
        coxpres <- coxpres[match(eids, gene1), ]
        coxpres$gene1 <- NULL

        # merge coxpres and top tables
        top_tables <- lapply(top_tables, function(x) {
            cbind(x, coxpres)
        })

        #fill in gene2 values with dprimes
        for (i in seq_along(top_tables)) {

            tt      <- top_tables[[i]]
            tt_full <- anal$top_tables[[i]]

            top_tables[[i]]$gene2_1 <- tt_full[tt$gene2_1, "dprime"]
            top_tables[[i]]$gene2_2 <- tt_full[tt$gene2_2, "dprime"]
            top_tables[[i]]$gene2_3 <- tt_full[tt$gene2_3, "dprime"]
            top_tables[[i]]$gene2_4 <- tt_full[tt$gene2_4, "dprime"]
            top_tables[[i]]$gene2_5 <- tt_full[tt$gene2_5, "dprime"]
            top_tables[[i]]$gene2_6 <- tt_full[tt$gene2_6, "dprime"]
            top_tables[[i]]$gene2_7 <- tt_full[tt$gene2_7, "dprime"]
            top_tables[[i]]$gene2_8 <- tt_full[tt$gene2_8, "dprime"]
            top_tables[[i]]$gene2_9 <- tt_full[tt$gene2_9, "dprime"]
            top_tables[[i]]$gene2_10 <- tt_full[tt$gene2_10, "dprime"]
        }

        #store results
        anal$top_tables <- top_tables
        anal
    })
    return(diff_exprs)
}


#---------------


# Select drug1, drug2, and combo contrasts for each GSE
#
# Used by setup_combo_data to ask user to select drug1, drug2, and combo
# contrasts for each GSE.
#
# @inheritParams setup_combo_data
#
# @return list (one per GSE) of lists (one per combination treatment)
#    with contrast names for each combination treatment in each GSE.

select_combo_data <- function(diff_exprs, prev_selections = list()){

    #setup selections
    selections <- vector("list", length(diff_exprs))
    names(selections) <- names(diff_exprs)

    #transfer prev_selections from GSEs also in diff_exprs
    prev_selections <- prev_selections[names(prev_selections) %in% names(selections)]
    selections[names(prev_selections)] <- prev_selections

    #select drug 1, drug 2, and combo from each study
    for (i in seq_along(diff_exprs)) {

        #check if previously selected
        gse_name <- names(diff_exprs)[i]
        if (!is.null(selections[[gse_name]])) next

        #if not, input selection
        choices <- names(diff_exprs[[i]]$top_tables)

        while (TRUE) {
            #select drug 1, 2, and combo contrast
            drug1 <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug 1")
            if (drug1 == "") {break}
            drug2 <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug 2")
            combo <- tcltk::tk_select.list(choices,
                                           title="Contrast for drug combo")

            #store selections
            selection <- c(drug1=drug1, drug2=drug2, combo=combo)
            selections[[gse_name]] <- append(selections[[gse_name]],
                                             list(selection))
        }
    }
    return(selections)
}

#---------------

# Creates final drug combination training data.frame
#
# Used by setup_combo_data to create final drug combination data.frame.
#
# @inheritParams setup_combo_data
#
# @return data.frame

combine_combo_data <- function(cbo_data, selections){

    # setup
    nrows  <- length(unlist(selections, recursive = FALSE))
    dr1_cols <- paste0('drug1_', colnames(cbo_data))
    dr2_cols <- paste0('drug2_', colnames(cbo_data))
    cbo_cols <- paste0('combo_', colnames(cbo_data))

    Xtrain   <- matrix(nrow = nrows,
                       ncol = 2 * ncol(cbo_data),
                       dimnames = list(1:nrows, c(dr1_cols, dr2_cols)))

    ytrain   <- matrix(nrow = nrows,
                       ncol = ncol(cbo_data),
                       dimnames = list(1:nrows, cbo_cols))


    #loop through each GSE
    n = 1
    cbo_data <- t(cbo_data)
    for (i in seq_along(selections)) {

        print(i)
        sels <- selections[[i]]

        #loop through each selection within each GSE
        for (sel in sels) {

            # add drug1 and drug2 values to Xtrain
            Xtrain[n, dr1_cols] <- cbo_data[, make.names(sel['drug1'])]
            Xtrain[n, dr2_cols] <- cbo_data[, make.names(sel['drug2'])]

            # add combo values to ytrain
            ytrain[n, ] <- cbo_data[, make.names(sel['combo'])]

            n <- n + 1
        }
    }
    return(list(X = Xtrain, y = ytrain))
}


#-----------------

# Add metaMA effectsize values to top tables.
#
# Used internally by \code{setup_combo_data} and \code{\link[crossmeta]{es_meta}}
# to add moderated unbiased standardised effect sizes (dprimes) to top tables
# from differential expression analysis.
#
# @param diff_exprs Result from call to \code{\link[crossmeta]{diff_expr}}.
# @param cols Columns from \code{\link[metaMA]{effectsize}} result to add to
#    top tables.
#
# @export
# @seealso \link[crossmeta]{diff_expr}, \link[crossmeta]{es_meta}.
#
# @return diff_exprs with specified columns added to top_tables for each contrast.
#
# @examples
# library(crossmeta)
# library(lydata)
#
# # location of raw data
# data_dir <- system.file("extdata", package = "lydata")
#
# # load previous analysis for eset
# anal <- load_diff("GSE9601", data_dir)
#
# # add dprime and vardprime to top tables
# anal <- add_es(anal)


add_es <- function(diff_exprs, cols = c("dprime", "vardprime")) {

    for (i in seq_along(diff_exprs)) {

        # get study degrees of freedom and group classes
        study <- diff_exprs[[i]]

        df <- study$ebayes_sv$df.residual + study$ebayes_sv$df.prior
        classes <- Biobase::pData(study$eset)$group

        for (con in names(study$top_tables)) {
            # get group names for contrast
            groups <- gsub("GSE.+?_", "", con)
            groups <- strsplit(groups, "-")[[1]]

            # get sample sizes for groups
            ni <- sum(classes == groups[2])
            nj <- sum(classes == groups[1])

            # bind effect size values with top table
            tt <- study$top_tables[[con]]
            es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
            tt <- cbind(tt, es)

            # store results
            study$top_tables[[con]] <- tt
        }
        diff_exprs[[i]] <- study
    }
    return(diff_exprs)
}
