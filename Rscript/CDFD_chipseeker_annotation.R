library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
library(clusterProfiler)
library(org.Hs.eg.db)


peaks <- readPeakFile("")

peaks <- readPeakFile("CDFD_final/tag_count/Control_h2bub_final_sorted.bed")

peaks <- readPeakFile("CDFD_final/tag_count/shp1_control_final_sorted.bed")

peaks <- readPeakFile("CDFD_final/tag_count/shp1_KO_Pol_II_final_sorted.bed")

peaks <- readPeakFile("MACS_folder/h2bub_chip_new/ctrl2_shp1/diff_peak.bed")

peaks <- readPeakFile("MACS_folder/shp1_chip/shp1_Ctrl_input/shp1_R1_R2_overlapped.bed")

peaks <- readPeakFile("MACS_folder/shp1_chip/shp1_Ctrl_input/shp1_ctrl_input/shp1_R1_peak.bed_flt.bed")

peaks <- readPeakFile("MACS_folder/pol2_chip/Pol_2pausing/pol2_diff_peak.bed")

peaks <- readPeakFile("emborep/bed/SHP-BL_peaks.filtered.bed")


print(peaks)


#peakAnno <- annotatePeak("merge/blast.bed", tssRegion=c(-3000, 3000),
#                        TxDb=txdb, annoDb="org.Hs.eg.db")

peakAnno <- annotatePeak(peaks, tssRegion=c(-2000, 2000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno)

plotAnnoBar(peakAnno)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

#files <- list("merge/HSC.bed","merge/pHSC.bed","merge/LSC.bed","merge/blast.bed")
######################################################################################################################3
files <- list(HSC = "merge/HSC.bed",CMP = "merge/CMP.bed",GMP = "merge/GMP.bed", Mono = "merge/Mono.bed")

files <- list(Control = "CDFD_final/tag_count/Control_h2bub_final_sorted.bed",shp1 = "CDFD_final/tag_count/shp1_control_final_sorted.bed",
              shp1_ko ="CDFD_final/tag_count/shp1_KO_Pol_II_final_sorted.bed")




files <- list(Ctrl = "MACS_folder/h2bub_chip_new/Ctrl_2_shp_ko2/Ctrl_Unique.bed",
              shp1_ko = "MACS_folder/h2bub_chip_new/Ctrl_2_shp_ko2/shp_kd_Unique.bed",
              Common = "MACS_folder/h2bub_chip_new/Ctrl_2_shp_ko2/Ctrl_shp1_kd_common.bed")



files <- list(shp1 = "MACS_folder/shp1_chip/shp1_Ctrl_input/shp1_ctrl_input/shp1_R2_peak.bed_flt.bed",
              shp2 = "MACS_folder/shp1_chip/shp1_Ctrl_input/shp1_ctrl_input/shp1_R2_peak.bed_flt.bed")


#################################################################### EMBOREP revision
files <- list(SHP_AB = "emborep/IgG_bed//SHP-AB_peaks.filtered.bed",
              SHP_BL = "emborep/IgG_bed//SHP-BL_peaks.filtered.bed",
              SHP_BI = "emborep/IgG_bed//SHP-BI_peaks.filtered.bed")


###################################################################



#files <- list(Blast = "merge/LSC.bed")

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, annoDb="org.Hs.eg.db",
                       tssRegion=c(-2000, 2000), verbose=FALSE)

plotAnnoPie(peakAnnoList)
plotAnnoBar(peakAnnoList)

promoter <- getPromoters(TxDb=txdb, upstream=1000, downstream=1000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000),weightCol="V4")
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), facet="row")

#####################################################################################
SHP_AB <- as.data.frame(peakAnnoList[["SHP_AB"]]@anno)
SHP_BL <- as.data.frame(peakAnnoList[["SHP_BL"]]@anno)
SHP_BI <- as.data.frame(peakAnnoList[["SHP_BI"]]@anno)

clust1 = as.data.frame(peakAnno@anno)

write_tsv(SHP_AB,"SHP_AB_IgG_+-2kb_peaks.txt")


#<1000 & > -1000


Prom_shp1_ctrl2 <- clust1 %>% filter( distanceToTSS < 1000  & distanceToTSS > -1000)

write_tsv(Prom_shp1_ctrl2,"Pol2_diff_annotated_promoters.txt")

#######################################################################################

library(ChIPpeakAnno)

HSC <- readPeakFile("merge/HSC.bed")
CMP <- readPeakFile("merge/CMP.bed")
GMP <- readPeakFile("merge/GMP.bed")
Mono <- readPeakFile("merge/Mono.bed")


##########################################################################################################################
promoter <- getPromoters(TxDb=txdb, upstream=2000, downstream=2000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)
##
## to speed up the compilation of this vigenette, we load a precaculated tagMatrixList
#data("tagMatrixList")

pdf("All_tagmatrix_atac.pdf", height = 10, width = 10)

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
dev.off()

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=10, facet="row")

#tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE) 

pdf("feature_distribution_atac.pdf", height = 15, width = 10)

plotAnnoBar(peakAnnoList)
dev.off()

pdf("feature_TSS_atac.pdf", height = 15, width = 10)
plotDistToTSS(peakAnnoList)
dev.off()

pdf("feature_UPSET_atac.pdf", height = 15, width = 10)

upsetplot(peakAnnoList, vennpie=FALSE)


dev.off()

pdf("Indiviual_sample_readcount_freq_atac.pdf", height = 15, width = 10)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=10, facet="row")
dev.off()



genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           #fun           = "enrichGO",
                           fun="enrichGO", OrgDb='org.Hs.eg.db',
                           pvalueCutoff  = 0.005,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "Pathway Enrichment Analysis")

###########################################################
compKEGG <- compareCluster(geneCluster = genes, 
                           fun = "enrichKEGG",
                           organism = "human",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")

dotplot(compKEGG, showCategory = 15, title = "Pathway Enrichment Analysis")




genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

HSC_annot <- as.data.frame(peakAnnoList[["HSC"]]@anno)
pHSC_annot <- as.data.frame(peakAnnoList[["pHSC"]]@anno)
LSC_annot <- as.data.frame(peakAnnoList[["LSC"]]@anno)
Blast_annot <- as.data.frame(peakAnnoList[["Blast"]]@anno)


# Get unique entrez gene Ids
entrezids <- unique(Blast_annot$geneId)

# Get hg19 entrez to gene symbol mappings
entrez2gene <- grch38 %>% 
  filter(entrez %in% entrezids) %>% 
  dplyr::select(entrez, symbol)

# Match to each annotation dataframe
m <- match(nanog_annot$geneId, entrez2gene$entrez)
nanog_annot <- cbind(nanog_annot[,1:14], geneSymbol=entrez2gene$symbol[m], nanog_annot[,15:16])

write.table(nanog_annot, file="results/Nanog_annotation.txt", sep="\t", quote=F, row.names=F)


ego <- enrichGO(gene = entrezids, 
                keyType = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)


cluster_summary <- data.frame(compKEGG)
write.csv(cluster_summary, "clusterProfiler_All_sample_atac_seq.csv")

dotplot(ego, showCategory=20)

