#BioManager::install(c("sangerseqR","annotate")) 

library(sangerseqR)

folder_path <- "Downloads/11108409869_SCF_SEQ_ABI/"

list_files <- list.files(path = folder_path,
                         pattern = ".ab1",
                         full.names = TRUE)

for(f in list_files) {
  abif <-read.abif(f)
  seq <- sangerseq(abif)
  str(seq)
  # chromatogram(seq, width = 500, 
  #              height = 2, 
  #              showcalls = "both")
  SeqX<-makeBaseCalls(seq)
  ## blast
  SeqXBlastDF<-blastSequences(paste(SeqX@primarySeq),as='data.frame', timeout=99)
  
  print(SeqXBlastDF$Hit_def)
  
}
