setwd("~/Documents/Batcave/GEO/ccdata/data-raw/")

# get compounds
cmap_instances <- read.table("raw/cmap_instances_02.csv",
                             header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

compounds <- unique(cmap_instances$cmap_name)

# change Prestwick to be like pubchem
compounds <- gsub('^Prestwick-', 'Prestwick_', compounds)

# adapted from webchem::get_cid
get_cid <- function (query, from = "name", verbose = TRUE, arg = NULL, ...) {
  foo <- function(query, from, first, verbose, ...) {
    cat('Working on:', query, '...\n')
    prolog <- "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    input <- paste0("/compound/", from, "/property/CanonicalSMILES")
    output <- "/JSON"
    if (!is.null(arg))
      arg <- paste0("?", arg)
    qurl <- paste0(prolog, input, output, arg)
    if (verbose)
      message(qurl)
    Sys.sleep(0.2)
    cont <- try(httr::content(httr::POST(qurl, body = paste0(from, "=", query), type = "text", encoding = "UTF-8")), silent = TRUE)
    if (inherits(cont, "try-error")) {
      warning("Problem with web service encountered... Returning NA.")
      return(NA)
    }
    out <- unlist(cont$PropertyTable$Properties[[1]])
    names(out) <- NULL
    return(out)
  }
  out <- lapply(query, foo, from = from, first = first, verbose = verbose)
  out <- setNames(out, query)
  if (first)
    out <- unlist(out)
  return(out)
}

# get pubchem cid and canonical smiles
cids <- get_cid(compounds)
pubids <- unlist(lapply(cids, `[`, 1))
smiles <- unlist(lapply(cids, `[`, 2))

annot <- data.frame(cmap_name=names(cids), row.names = names(cids))
annot$cid <- ''
annot$smiles <- ''

annot[names(pubids), 'cid'] <- pubids
annot[annot$cid == 'NULL', 'cid'] <- ''
annot[names(smiles), 'smiles'] <- smiles


# save missing for manual annotation
# googled e.g. 'BCB000038 DSigDB' then clicked 'InChIKey' or 'Links'
write.csv(annot[annot$cid == '', ], 'cmap_es/annot_missing.csv')

# replace missing
miss <- read.csv('cmap_es/annot_missing.csv', row.names = 1, na.strings = '', stringsAsFactors = FALSE)
miss[is.na(miss)] <- ''
annot[row.names(miss), 'cid'] <- miss$cid
annot[row.names(miss), 'smiles'] <- miss$smiles

# change Prestwick back
annot$cmap_name <- gsub('^Prestwick_', 'Prestwick-', annot$cmap_name)
row.names(annot) <- annot$cmap_name

saveRDS(annot, 'cmap_es/annot.rds')
