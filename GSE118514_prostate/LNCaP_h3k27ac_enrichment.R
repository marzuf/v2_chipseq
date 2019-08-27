setDir=""

require(GenomicRanges)

hicds <- "ENCSR346DCU_LNCaP_40kb" 
exprds <- "TCGAprad_norm_prad"
histMark <- "H3K27ac"

plotCex <- 1.2
outFolder <- "LNCAP_H3K27AC_ENRICHMENT"
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

signifPeakFile <- "GSM3358195_LNCaP_H3K27Ac_signifpeaks.gff"
stopifnot(file.exists(signifPeakFile))

tadFile <- file.path("..", "..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA", hicds, "genes2tad/all_assigned_regions.txt")
stopifnot(file.exists(tadFile))
tad_DT <- read.delim(tadFile, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
tad_DT <- tad_DT[grep("_TAD", tad_DT$region),]

peakDT <- read.delim(signifPeakFile, header=F, stringsAsFactors = FALSE)
peakDT <- peakDT[,c(1,4,5, 6)]
colnames(peakDT) <- c("chromo", "start", "end", "peakValue")  # !!! THIS IS NOT PEAK VALUE; DO NOT USE
peakDT <- peakDT[order(peakDT$chromo, peakDT$start, peakDT$end),]
peakDT$region <- paste0("peak", 1:nrow(peakDT))
stopifnot(is.numeric(peakDT$start))
stopifnot(is.numeric(peakDT$end))
stopifnot(is.numeric(peakDT$peakValue))

tad_GR <- GRanges(seqnames=tad_DT$chromo, ranges=IRanges(start=tad_DT$start, end=tad_DT$end, names=tad_DT$region))
peak_GR <- GRanges(seqnames=peakDT$chromo, ranges=IRanges(start=peakDT$start, end=peakDT$end, names=peakDT$region))

# determine which features from the query overlap which features in the subject
overlap_GR <- findOverlaps(query=peak_GR, subject=tad_GR)
tadID <- names(tad_GR[subjectHits(overlap_GR)])
peakID <- names(peak_GR[queryHits(overlap_GR)])

overlapDT <- data.frame(
  tadID = tadID,
  peakID = peakID,
  overlapBp = width(pintersect(tad_GR[tadID], peak_GR[peakID])),
  stringsAsFactors = FALSE)
stopifnot(overlapDT$tadID %in% tad_DT$region)
stopifnot(overlapDT$peakID %in% peakDT$region)

# peak values
peak_overlap_DT <- overlapDT
colnames(peak_overlap_DT)[colnames(peak_overlap_DT) == "peakID"] <- "region"
tad_peak_values_DT <- merge(peak_overlap_DT, peakDT, by="region")
tad_peakValues_DT <- aggregate(peakValue ~ tadID, data=tad_peak_values_DT, sum)
colnames(tad_peakValues_DT)[colnames(tad_peakValues_DT) == "tadID"] <- "region"

# take the total overlap for each TAD
tad_overlap_DT <- aggregate(overlapBp ~ tadID, data=overlapDT, sum)
colnames(tad_overlap_DT)[colnames(tad_overlap_DT) == "tadID"] <- "region"
tad_overlapDT <- merge(tad_DT, tad_overlap_DT, by="region", all=TRUE)
tad_overlapDT$tad_size <- tad_overlapDT$end - tad_overlapDT$start + 1
tad_overlapDT$overlapBp[is.na(tad_overlapDT$overlapBp)] <- 0
stopifnot(tad_overlapDT$overlapBp <= tad_overlapDT$tad_size)
tad_overlapDT$overlap_ratio <- tad_overlapDT$overlapBp/tad_overlapDT$tad_size

# plot(density(tad_overlapDT$overlap_ratio))
final_overlap_dt <- merge(curr_final_dt, tad_overlapDT[,c("region", "overlap_ratio")], by="region" )
final_overlap_dt <- merge(final_overlap_dt, tad_peakValues_DT[,c("region", "peakValue")], by="region")
final_overlap_dt <- final_overlap_dt[order(final_overlap_dt$adjPvalComb),]
final_overlap_dt$signifPval_0.001 <- final_overlap_dt$adjPvalComb <= 0.001
final_overlap_dt$adjPvalComb_log10 <- -log10(final_overlap_dt$adjPvalComb)

all_box_vars <- c("signifFDR_0.1", "signifFDR_0.2", "signifPval_0.001")
# all_y_vars <- c("overlap_ratio", "peakValue")  # !! this is not the peak value, this is for drawing
all_y_vars <- c("overlap_ratio")  #
all_x_vars <- c("meanLogFC", "adjPvalComb_log10")

y_var = all_y_vars[1]
box_var = all_box_vars[1]
x_var = all_x_vars[1]

for(y_var in all_y_vars) {
  for(box_var in all_box_vars) {
    cat("... start ", y_var, " vs. ", box_var, "\n")
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", y_var, "_by_", box_var, "_boxplot.", plotType))
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
    
    outFile <- file.path(outFolder, paste0(hicds, "_", exprds, "_", y_var, "_vs_", x_var, "_densplot.", plotType))
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








