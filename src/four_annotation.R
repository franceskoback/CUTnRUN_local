# load in required packages
library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(clusterProfiler)
library(ReactomePA)
library(tidyverse)
library(ggupset)
library(ggimage)

# which human genome version to use
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# set the output directory
#setwd("<path>/annotated_stringent_peaks")
setwd("/Users/fkoback/Documents/Projects/CUTnRUN/results_3_calledpeaks/results_3_calledpeaks")


# locations of all input files
#files <- list.files(path = "<path>/03_peak_calling_v1beds", 
#                    pattern = "stringent.auc.threshold.merge.bed",
#                    full.names = TRUE,
#                    include.dirs = FALSE)

files <- list.files(pattern = "stringent.bed$",
                    full.names = TRUE,
                    include.dirs = FALSE)

# Make a list of all the names because.... thats what it wants
#length_files=length(files)

files<- as.list(files)
#files <- as.list(c(files[1], files[2], files[3], files[4]))
# name each of the peak files accordingly for batch processing later
names(files)<-c(substring(sub(".stringent.bed*", "",files),3))

# iterate through and make individual plots for each peaks file
# current_num <- 1
for( current_num in 1:length(files) ) {
  # what is the name for the current sample
  current_name <- names(files)[current_num]
  
  # name of the current peak file
  current_peak_file_name <- basename( files[[current_num]] ) 
  message(current_peak_file_name)
  #  read in a file
  peak <- readPeakFile( files[[current_num]] )
  newStyle <- mapSeqlevels(seqlevels(peak), "UCSC")
  peak <- renameSeqlevels(peak, newStyle)
  
  # TSS peak binding
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  
  tagMatrix <- getTagMatrix(peak, windows=promoter)
  
  # average profile straight from bed
  print("average peaks")
  png(filename = paste0(current_peak_file_name, " - Average Peaks.png"),
      width = 1920, height = 1080,
      units = "px", pointsize = 6, res = 450)
  #plot( plotAvgProf2(files[[current_num]], TxDb=txdb, upstream=3000, downstream=3000,
  #                   xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95)
  #)
  plot(plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000,
                                        xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95))
  
#ggsave(".plots/Average_Peaks.png")
  dev.off()
  
  # peak annotation
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  paste("saving annotation texts")
  write.table(peakAnno, file = paste0(current_peak_file_name, " - Annotated Peaks.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(peakAnno@annoStat, file = paste0(current_peak_file_name, " - Annotated Peaks Distribution.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # visualise
  print("pie chart")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - PieChart.png"),
      width = 1920, height = 1080,
      units = "px", pointsize = 6, res = 450)
  plotAnnoPie(peakAnno)
  dev.off()
  
  print("bar chart")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - BarChart.png"),
      width = 1920, height = 750,
      units = "px", pointsize = 6, res = 450)
  plotAnnoBar(peakAnno)
  dev.off()
  
  print("venn diagram")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - VennDiagram.png"),
      width = 2500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  vennpie(peakAnno)
  dev.off()
  
  print("upset plot")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - UpsetPlot.png"),
      width = 3500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  upsetplot(peakAnno)
  dev.off()
  
  print("upset plot venn")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - UpsetPlot VennDiagram.png"),
      width = 3500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  upsetplot(peakAnno, vennpie=TRUE)
  dev.off()
  
  # TF binding relative to TSS
  print("tss binding")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding.png"),
      width = 2500, height = 750,
      units = "px", pointsize = 8, res = 450)
  plot( plotDistToTSS(peakAnno,
                      title="Distribution of transcription factor-binding loci\nrelative to TSS")
  )
  dev.off()
}

files <- list.files(pattern = "relaxed.bed$",
                    full.names = TRUE,
                    include.dirs = FALSE)

# Make a list of all the names because.... thats what it wants
#length_files=length(files)

files<- as.list(files)
#files <- as.list(c(files[1], files[2], files[3], files[4]))
# name each of the peak files accordingly for batch processing later
names(files)<-c(substring(sub(".relaxed.bed*", "",files),3))

# iterate through and make individual plots for each peaks file
# current_num <- 1
for( current_num in 1:length(files) ) {
  # what is the name for the current sample
  current_name <- names(files)[current_num]
  
  # name of the current peak file
  current_peak_file_name <- basename( files[[current_num]] ) 
  message(current_peak_file_name)
  #  read in a file
  peak <- readPeakFile( files[[current_num]] )
  newStyle <- mapSeqlevels(seqlevels(peak), "UCSC")
  peak <- renameSeqlevels(peak, newStyle)
  
  # TSS peak binding
  promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
  
  tagMatrix <- getTagMatrix(peak, windows=promoter)
  
  # average profile straight from bed
  print("average peaks")
  png(filename = paste0(current_peak_file_name, " - Average Peaks.png"),
      width = 1920, height = 1080,
      units = "px", pointsize = 6, res = 450)
  #plot( plotAvgProf2(files[[current_num]], TxDb=txdb, upstream=3000, downstream=3000,
  #                   xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95)
  #)
  plot(plotAvgProf2(peak, TxDb=txdb, upstream=3000, downstream=3000,
                    xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency", conf = 0.95))
  
  #ggsave(".plots/Average_Peaks.png")
  dev.off()
  
  # peak annotation
  peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                           TxDb=txdb, annoDb="org.Hs.eg.db")
  paste("saving annotation texts")
  write.table(peakAnno, file = paste0(current_peak_file_name, " - Annotated Peaks.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  write.table(peakAnno@annoStat, file = paste0(current_peak_file_name, " - Annotated Peaks Distribution.txt"),
              row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
  
  # visualise
  print("pie chart")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - PieChart.png"),
      width = 1920, height = 1080,
      units = "px", pointsize = 6, res = 450)
  plotAnnoPie(peakAnno)
  dev.off()
  
  print("bar chart")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - BarChart.png"),
      width = 1920, height = 750,
      units = "px", pointsize = 6, res = 450)
  plotAnnoBar(peakAnno)
  dev.off()
  
  print("venn diagram")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - VennDiagram.png"),
      width = 2500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  vennpie(peakAnno)
  dev.off()
  
  print("upset plot")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - UpsetPlot.png"),
      width = 3500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  upsetplot(peakAnno)
  dev.off()
  
  print("upset plot venn")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding - UpsetPlot VennDiagram.png"),
      width = 3500, height = 2500,
      units = "px", pointsize = 8, res = 450)
  upsetplot(peakAnno, vennpie=TRUE)
  dev.off()
  
  # TF binding relative to TSS
  print("tss binding")
  png(filename = paste0(current_peak_file_name, " - TSS Peak Binding.png"),
      width = 2500, height = 750,
      units = "px", pointsize = 8, res = 450)
  plot( plotDistToTSS(peakAnno,
                      title="Distribution of transcription factor-binding loci\nrelative to TSS")
  )
  dev.off()
}


