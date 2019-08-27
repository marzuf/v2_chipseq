options(scipen=100)
# Rscript create_bed_from_tad.R

require(foreach)
require(doMC)
registerDoMC(4)

hicds <- "GSE99051_786_O_40kb" 

exec <- "/mnt/ed4/marie/software/bigWigAverageOverBed"
stopifnot(file.exists(exec))

tadFile <- file.path("..", "..", "v2_Yuanlong_Cancer_HiC_data_TAD_DA", hicds, "genes2tad/all_assigned_regions.txt")
stopifnot(file.exists(tadFile))
tad_DT <- read.delim(tadFile, header=F, stringsAsFactors = FALSE, col.names=c("chromo", "region", "start", "end"))
tad_DT <- tad_DT[grep("_TAD", tad_DT$region),]


bed_tad_DT <- tad_DT[,c("chromo", "start", "end", "region")]
bed_tad_DT$start <- bed_tad_DT$start-1
outFile <- file.path(paste0(hicds, "_all_assigned_regions_bed.txt"))
write.table(bed_tad_DT, file=outFile, sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE, quote=FALSE)

bedfile <- outFile

all_bwfiles <- c("GSM2585472_ChIP-Seq_786-O_H3K27ac_1.bw",
                 "GSM2585473_ChIP-Seq_786-O_H3K27ac_2.bw",
                 "GSM2585474_ChIP-Seq_786-O_H3K4me1_1.bw",
                 "GSM2585475_ChIP-Seq_786-O_H3K4me1_2.bw"
                 )

foo <- foreach(bwfile = all_bwfiles) %dopar% {
  
  
  outfile <- gsub(".bw$", ".bed", basename(bwfile))
  
  stopifnot(file.exists(bwfile))
  stopifnot(file.exists(bedfile))
  
  mycmd <- paste(exec, bwfile, bedfile, outfile)
  cat(mycmd, "\n")
  system(mycmd)
  
}
# http://www.columbia.edu/~ii2135/readme.txt
# To use the utility, you must prepare a bed file, a tab-separated 4-column file.  
# The first column is the chromosome, the second is the zero-based coordinate of the position of interest, the third is that zero-based coordinate plus one, and the fourth is a unique identifier for the position.
# The output columns are:
#   name - name field from bed, which should be unique
# size - size of bed (sum of exon sizes
#                     covered - # bases within exons covered by bigWig
#                       sum - sum of values over all bases covered
#                     mean0 - average over bases with non-covered bases counting as zeroes
#                     mean - average over just covered bases