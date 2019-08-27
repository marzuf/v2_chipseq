# Rscript 786-O_enrichment.R

setDir=""

hicds <- "GSE99051_786_O_40kb" 
exprds <- "TCGAkich_norm_kich"


all_histMarks <- c("H3K27ac", "H3K4me1")

plotCex <- 1.2
outFolder <- "786-0_ENRICHMENT"
dir.create(outFolder, recursive = TRUE)

plotType <- "png"
myHeight <- ifelse(plotType=="png", 400, 7)
myWidth <- myHeight
plotCex <- 1.2

source("../../2_Yuanlong_Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../../Cancer_HiC_data_TAD_DA/utils_fct.R")

final_table_file <- file.path("..", "..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA", "CREATE_FINAL_TABLE", "all_result_dt.Rdata")
stopifnot(file.exists(final_table_file))
final_dt <- get(load(final_table_file))
curr_final_dt <- final_dt[final_dt$hicds == hicds & final_dt$exprds == exprds,]
stopifnot(nrow(curr_final_dt) > 0)
stopifnot(!duplicated(curr_final_dt$region))

tadFile <- file.path("..", "..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA", hicds, "genes2tad/all_assigned_regions.txt")
stopifnot(file.exists(tadFile))
tad_DT <- read.delim(tadFile, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
tad_DT <- tad_DT[grep("_TAD", tad_DT$region),]

# GSM2585472_ChIP-Seq_786-O_H3K27ac_1.bw
# GSM2585473_ChIP-Seq_786-O_H3K27ac_2.bw
# GSM2585474_ChIP-Seq_786-O_H3K4me1_1.bw
# GSM2585475_ChIP-Seq_786-O_H3K4me1_2.bw

#1 chrom - Name of the chromosome (or contig, scaffold, etc.).
#2 chromStart - The starting position of the feature in the chromosome or scaffold. The first base in a chromosome is numbered 0.
#3 chromEnd - The ending position of the feature in the chromosome or scaffold. The chromEnd base is not included in the display of the feature. For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.
#4 name - Name given to a region (preferably unique). Use "." if no name is assigned.
#5 score - Indicates how dark the peak will be displayed in the browser (0-1000). If all scores were "'0"' when the data were submitted to the DCC, the DCC assigned scores 1-1000 based on signal value. Ideally the average signalValue per base spread is between 100-1000.
#6 strand - +/- to denote strand or orientation (whenever applicable). Use "." if no orientation is assigned.
#7 signalValue - Measurement of overall (usually, average) enrichment for the region.
#8 pValue - Measurement of statistical significance (-log10). Use -1 if no pValue is assigned.
#9 qValue - Measurement of statistical significance using false discovery rate (-log10). Use -1 if no qValue is assigned.
#10 peak - Point-source called for this peak; 0-based offset from chromStart. Use -1 if no point-source called.

for(histMark in all_histMarks) {
  

  # The output columns are:
  #   name - name field from bed, which should be unique
  # size - size of bed (sum of exon sizes
  #                     covered - # bases within exons covered by bigWig
  #                       sum - sum of values over all bases covered
  #                     mean0 - average over bases with non-covered bases counting as zeroes
  #                     mean - average over just covered bases
  if(histMark == "H3K4me1") {
    file1 <- paste0("GSM2585474_ChIP-Seq_786-O_", histMark, "_1.bed")
    file2 <- paste0("GSM2585475_ChIP-Seq_786-O_", histMark, "_2.bed")
  } else if(histMark == "H3K27ac") {
    file1 <- paste0("GSM2585472_ChIP-Seq_786-O_", histMark, "_1.bed")
    file2 <- paste0("GSM2585473_ChIP-Seq_786-O_", histMark, "_2.bed")
  } else {
    stop("error histmark")
  }
  
  dt1 <- read.delim(file1, header=FALSE, stringsAsFactors = FALSE)
  dt2 <- read.delim(file2, header=FALSE, stringsAsFactors = FALSE)

  peakDT <- rbind(dt1,dt2)
  colnames(peakDT) <- c("region", "size","covered", "sum", "mean0", "mean")
  
  
  
  final_overlap_dt <- merge(curr_final_dt, peakDT, by="region" )
  
  final_overlap_dt <- final_overlap_dt[order(final_overlap_dt$adjPvalComb),]
  final_overlap_dt$signifPval_0.001 <- final_overlap_dt$adjPvalComb <= 0.001
  final_overlap_dt$adjPvalComb_log10 <- -log10(final_overlap_dt$adjPvalComb)
  
  all_box_vars <- c("signifFDR_0.1", "signifFDR_0.2", "signifPval_0.001")
  
  all_x_vars <- c("meanLogFC", "adjPvalComb_log10")
  
  all_y_vars <-  c("sum", "mean0", "mean")
  
  y_var = all_y_vars[1]
  box_var = all_box_vars[1]
  x_var = all_x_vars[1]
  
  for(y_var in all_y_vars) {
    for(box_var in all_box_vars) {
      cat("... start ", y_var, " vs. ", box_var, "\n")
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", histMark, "_", y_var, "_by_", box_var, "_boxplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      boxplot(as.formula(paste0(y_var, " ~  ", box_var)), ylab=y_var, xlab=box_var,data=final_overlap_dt, main=paste0(hicds, " - ", exprds, " - ", histMark))  
      mtext(side=3, text=paste0(box_var, ": ", sum(final_overlap_dt[,paste0(box_var)]), "/", nrow(final_overlap_dt), "\n"))
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
    }
    for(x_var in all_x_vars) {
      cat("... start ", y_var, " vs. ", x_var, "\n")
      
      myx <- final_overlap_dt[,paste0(x_var)]
      myy <- final_overlap_dt[,paste0(y_var)]
      
      outFile <- file.path(outFolder, paste0(hicds, "_", exprds,"_", histMark, "_", y_var, "_vs_", x_var, "_densplot.", plotType))
      do.call(plotType, list(outFile, height=myHeight, width=myWidth))
      densplot(
        x=(myx),
        y=(myy),
        xlab = paste0(x_var, ""),
        ylab = paste0(y_var, ""),
        main = paste0( y_var, " vs. ", x_var, ""),
        # sub = paste0(hicds, " - ", exprds),
        cex.lab=plotCex,
        cex.axis = plotCex
      )
      addCorr(x = myx, y = myy, bty="n", legPos = "topright")
      mtext(side=3, text = paste0(hicds, " - ", exprds, " - ", histMark))
      # chr12_TAD220; logFC > 0
      final_overlap_dt <- final_overlap_dt[order(abs(final_overlap_dt[,paste0(x_var)]), decreasing = TRUE),]
      text(
        x = final_overlap_dt[,paste0(x_var)][1:3],
        y = final_overlap_dt[,paste0(y_var)][1:3],
        labels = final_overlap_dt[,paste0("region")][1:3],
        cex=0.6,
        pos=3
      )
      final_overlap_dt[1:3,]
      final_overlap_dt <- final_overlap_dt[order(final_overlap_dt[,paste0(y_var)], decreasing = TRUE),]
      text(
        x = final_overlap_dt[,paste0(x_var)][1:3],
        y = final_overlap_dt[,paste0(y_var)][1:3],
        labels = final_overlap_dt[,paste0("region")][1:3],
        cex=0.6,
        pos=3
      )
      foo <- dev.off()
      cat(paste0("... written: ", outFile, "\n"))
      final_overlap_dt[1:3,]
    }
  }
}


