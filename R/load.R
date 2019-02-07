# load.R

# Loads COSMIC mutation file and adds HGNC gene identifiers to data

# mut_filename is the filepath which contains the large mutation data file.


COSMIC <- function(mut_filename = "./data/CosmicCompleteTargetedScreensMutantExport.tsv"){
  colTypes <- cols_only(
    `Gene name` = col_character(),
    `Gene CDS length` = col_double(),
    ID_sample = col_double(),
    ID_tumour = col_double(),
    `Mutation ID` = col_character(),
    `Mutation CDS` = col_character(),
    `Primary site` = col_character()
  )
  
  hgnc <- readr::read_tsv("./data/hgnc.tsv")
  mut <- readr::read_tsv(mut_filename, col_types = colTypes)
  mut <- dplyr::rename(mut, Gene = `Gene name`, SampleID = ID_sample,
                       Mutation = `Mutation CDS`, Site = `Primary site`)
  
  
  sum(is.na(mut$Gene)) # 0
  HGNC <- unique(mut$Gene)
  sum(HGNC %in% hgnc$`Approved symbol`)/length(HGNC)
  head(HGNC[!(HGNC %in% hgnc$`Approved symbol`)])
  matchTable <- data.frame(original = HGNC, stringsAsFactors = F)
  matchTable$Replace <- NA
  matchTable$Replace[HGNC %in% hgnc$`Approved symbol`] <- HGNC[HGNC %in% hgnc$`Approved symbol`]
  
  
  missing <- which(!(HGNC %in% hgnc$`Approved symbol`))
  
  Names_ <- sapply(missing, function(x) strsplit(HGNC[x], "_")[[1]][1])
  sel <- Names_ %in% hgnc$`Approved symbol`
  matchTable$Replace[missing[sel]] <- Names_[sel]
  mean(is.na(matchTable$Replace))
  
  missing <- which(is.na(matchTable$Replace))
  previous <- sapply(missing, function(x) grep(HGNC[x], hgnc$`Previous symbols`)[1])
  reset <- which(is.na(previous))
  previous[is.na(previous)] <- 1
  matchTable$Replace[missing] <- hgnc$`Approved symbol`[previous]
  matchTable$Replace[reset] <- NA
  mean(is.na(matchTable$Replace))
  
  mut$newSymbol <- matchTable$Replace[match(mut$Gene, matchTable$original)]
  mut$Mutation <- ifelse(is.na(mut$Mutation), 0, 1)
  mean(mut$Mutation)
}