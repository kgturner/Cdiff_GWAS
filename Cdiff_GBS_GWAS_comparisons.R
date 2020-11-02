#Compare GWAS hit and surrounding to:
#1) genes rapidly evolving in invasive Cdiff (Hodgins et al 2015)
#2) genes differentially expressed between native and invasive Cdiff, in drought and control (Turner et al 2017)
#11/20/2018

# GWAS_blue <- read.delim("BlueLine_chr12.10895280-11723974.txt")
GWAS_1Mb <- read.delim("genes_near_gwas.1Mbp.summary_table_20200602.txt")
Hodgins <- read.delim("Hodgins_2015_AllCdiff.txt") #combined Cdiff related results from Hodgins
DE <- read.delim("DE_LRTs.txt") #Qvals for LRT of genes with at least one sig fixed effect from Turner 2017
Expr <- read.table("lme4_qval_PC1.txt", header=T, sep="\t") #contigs and qvals for all LRTs from Turner 2017

####match by TAIR hits?####
#overlap in TAIR best hits with decimal
GWAS_athal <- subset(GWAS_1Mb, !is.na(Athal_BEST_HIT))
overlap_HodginsAthal <- subset(GWAS_athal, Athal_BEST_HIT %in% Hodgins$AthalBestHit)

#overlap in TAIR best hits, NO DECIMAL
DE_athal1 <- subset(DE, !is.na(TAIR_1))
GWAS_athal <- subset(GWAS_1Mb, !is.na(Athal_BEST_HIT))

overlap_DEathal1 <- subset(GWAS_athal, Athal_BEST_HIT_noDec %in% DE_athal1$TAIR_1)
DE_athal2 <- subset(DE, !is.na(TAIR_2))
overlap_DEathal2 <- subset(GWAS_athal, Athal_BEST_HIT_noDec %in% DE_athal2$TAIR_2)

overlap_DE <- rbind(overlap_DEathal1, overlap_DEathal2)
write.table(overlap_DE, "Cdiff_GWAS_overlap_Turner2017.txt", sep = "\t") #overlap between ref annotations and DE

overlap_LRT1 <- subset(DE, TAIR_1 %in% overlap_DE$Athal_BEST_HIT_noDec)
overlap_LRT2 <- subset(DE, TAIR_2 %in% overlap_DE$Athal_BEST_HIT_noDec)
overlap_LRT <- rbind(overlap_LRT1, overlap_LRT2)

write.table(overlap_LRT, "Cdiff_GWAS_overlap_Turner2017_analysis.txt", sep = "\t")

# subset(DE, TAIR_1 %in% "AT1G29050") #tair hit from DE with no contigs in GWAS table?

####match by contig name####
GWAS_overlap_vector <- c(as.character(GWAS_1Mb$Transcripts.DK_TR001.1L.1),as.character(GWAS_1Mb$Transcripts.DK_US022.31E.1),
                         as.character(GWAS_1Mb$Transcripts.DK_TR001.1L.2),as.character(GWAS_1Mb$Transcripts.DK_US022.31E.2),
                         as.character(GWAS_1Mb$Transcripts.DK_TR001.1L.3),as.character(GWAS_1Mb$Transcripts.DK_US022.31E.3),
                         as.character(GWAS_1Mb$Transcripts.DK_TR001.1L.4),as.character(GWAS_1Mb$Transcripts.DK_US022.31E.4))
GWAS_overlap_vector <-GWAS_overlap_vector[!is.na(GWAS_overlap_vector)]
GWAS_overlap_vector <- unique(GWAS_overlap_vector)
GWAS_overlap_vector

#overlap using contig/transcript id from Turner 2017, whether sig or not, match by Conting no.
overlap_Expr1 <- subset(Expr, Contig %in% GWAS_overlap_vector ) #subset of contigs and qvals for all LRTs from Turner 2017
write.table(overlap_Expr1,"Cdiff_GWAS_overlapAll_Turner2017_analysis.txt", sep = "\t")

overlap_Expr2 <- subset(GWAS_1Mb, Transcripts.DK_TR001.1L.1 %in% overlap_Expr1$Contig |
                          Transcripts.DK_TR001.1L.2 %in% overlap_Expr1$Contig |
                          Transcripts.DK_TR001.1L.3 %in% overlap_Expr1$Contig |
                          Transcripts.DK_TR001.1L.4 %in% overlap_Expr1$Contig |
                          Transcripts.DK_US022.31E.1 %in% overlap_Expr1$Contig |
                          Transcripts.DK_US022.31E.2 %in% overlap_Expr1$Contig |
                          Transcripts.DK_US022.31E.3 %in% overlap_Expr1$Contig |
                          Transcripts.DK_US022.31E.4 %in% overlap_Expr1$Contig)

write.table(overlap_Expr2,"Cdiff_GWAS_overlapAll_Turner2017.txt", sep = "\t")

#check GO ids?E