library(dplyr)
library(readr)
library(ggplot2)

mut_filename <- "./data/CosmicCompleteTargetedScreensMutantExport.tsv"

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

top_genes <- sort(table(mut$newSymbol), decreasing = T)[1:10]
x <- mut$newSymbol[mut$newSymbol %in% names(top_genes)]
ggplot() + geom_bar(aes(x = x)) + xlab("Gene Name")

ggplot(data = mut, aes(x = as.factor(Mutation), y = log(`Gene CDS length`))) + 
  geom_boxplot() + xlab("Mutation") + ylab("log(Gene Length)")

mut_rates <- tapply(mut$Mutation, mut$newSymbol, mean)
gene_length <- tapply(mut$`Gene CDS length`, mut$newSymbol, mean)
counts <- table(mut$newSymbol)
ggplot() + geom_point(aes(x = log(gene_length), y = mut_rates))

which(!(names(counts)== names(gene_length)))
model <- lm(mut_rates ~ log(gene_length), weights = counts^(1/2))
summary(model)

gene_set <- c("AMBRA1", "ATG14", "ATP2A1", "ATP2A2", "ATP2A3", "BECN1", "BECN2", 
              "BIRC6", "BLOC1S1", "BLOC1S2", "BORCS5", "BORCS6", "BORCS7", 
              "BORCS8", "CACNA1A", "CALCOCO2", "CTTN", "DCTN1", "EPG5", "GABARAP", 
              "GABARAPL1", "GABARAPL2", "HDAC6", "HSPB8", "INPP5E", "IRGM", 
              "KXD1", "LAMP1", "LAMP2", "LAMP3", "LAMP5", "MAP1LC3A", "MAP1LC3B", 
              "MAP1LC3C", "MGRN1", "MYO1C", "MYO6", "NAPA", "NSF", "OPTN", 
              "OSBPL1A", "PI4K2A", "PIK3C3", "PLEKHM1", "PSEN1", "RAB20", "RAB21", 
              "RAB29", "RAB34", "RAB39A", "RAB7A", "RAB7B", "RPTOR", "RUBCN", 
              "RUBCNL", "SNAP29", "SNAP47", "SNAPIN", "SPG11", "STX17", "STX6", 
              "SYT7", "TARDBP", "TFEB", "TGM2", "TIFA", "TMEM175", "TOM1", 
              "TPCN1", "TPCN2", "TPPP", "TXNIP", "UVRAG", "VAMP3", "VAMP7", 
              "VAMP8", "VAPA", "VPS11", "VPS16", "VPS18", "VPS33A", "VPS39", 
              "VPS41", "VTI1B", "YKT6")

subset <- mut[mut$Gene %in% gene_set,]
transform <- as.data.frame(tapply(subset$Mutation, list(subset$Gene, subset$Site), mean))
transform$gene <- rownames(transform)

plotable <- tidyr::gather(transform, "Site", "Ratio", 1:39,
                          na.rm = T)

ggplot(data = plotable, aes(x = gene, y = Ratio, color = Site)) +
  geom_point() 

readr::write_tsv(transform, "./data/Example.tsv")


