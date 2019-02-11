# load.R
# @Matthew McNeil
# Loads COSMIC mutation file and adds HGNC gene identifiers to data

# mut_filename is the filepath which contains the large mutation data file.


COSMIC <- function(mut_filename = "./data/CosmicCompleteTargetedScreensMutantExport.tsv"){
  source("./R/INCLUDE.R)
  # Set columns to be used for analysis
  colTypes <- cols_only(
    `Gene name` = col_character(),
    `Gene CDS length` = col_double(),
    ID_sample = col_double(),
    ID_tumour = col_double(),
    `Mutation ID` = col_character(),
    `Mutation CDS` = col_character(),
    `Primary site` = col_character()
  )
  
  # Read in HGNC symbols for cross reference
  hgnc <- readr::read_tsv("./data/hgnc.tsv")
  
  # Use readr to quickly read in the large file
  mut <- readr::read_tsv(mut_filename, col_types = colTypes)
  # Rename columns to avoid confusion
  mut <- dplyr::rename(mut, Gene = `Gene name`, SampleID = ID_sample,
                       Mutation = `Mutation CDS`, Site = `Primary site`)
  
  # Look for unique gene symbols only, avoid wasting time matching duplicates
  HGNC <- unique(mut$Gene)
  # Create data frame for matching listed identifier with real HGNC
  matchTable <- data.frame(original = HGNC, stringsAsFactors = F)
  matchTable$Replace <- NA
  # Set the correct ones
  matchTable$Replace[HGNC %in% hgnc$`Approved symbol`] <- HGNC[HGNC %in% hgnc$`Approved symbol`]
  
  # Now start investigating the missing identifiers
  missing <- which(!(HGNC %in% hgnc$`Approved symbol`))
  # Get the names which are an amalgamation of different identifiers
  Names_ <- sapply(missing, function(x) strsplit(HGNC[x], "_")[[1]][1])
  
  # Select the ones that now match and put them in our match table
  sel <- Names_ %in% hgnc$`Approved symbol`
  matchTable$Replace[missing[sel]] <- Names_[sel]
  
  # Genes that are still missing
  missing <- which(is.na(matchTable$Replace))
  # Look for them in previous symbols and then match them into the data frame
  previous <- sapply(missing, function(x) grep(HGNC[x], hgnc$`Previous symbols`)[1])
  reset <- which(is.na(previous))
  previous[is.na(previous)] <- 1
  matchTable$Replace[missing] <- hgnc$`Approved symbol`[previous]
  matchTable$Replace[reset] <- NA
  
  # # Replace the newly found gene names, and make mutation a binary variable
  mut$newSymbol <- matchTable$Replace[match(mut$Gene, matchTable$original)]
  mut$Mutation <- ifelse(is.na(mut$Mutation), 0, 1)
  return(mut)
}
