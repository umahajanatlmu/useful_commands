library("tidyverse")
library("FDb.InfiniumMethylation.hg19")
library("biomaRt")
library("here")


get_TCGA_ICGC_data <- function(project=c("PAAD-US","PACA-AU", "PACA-CA"),
                               clinical=TRUE,
                               rna_seq =TRUE,
                               methylation= FALSE,
                               mirna = FALSE,
                               cnv=TRUE,
                               somatic_mut=TRUE
) {
  data <- list()
  compiled_dataset <- c()
  # create download directory and set it
  .exdir = 'tmp'
  dir.create(.exdir)
  
  url="https://dcc.icgc.org/api/v1/download?fn=/current/Projects/"
  
  if (isTRUE(cnv) || isTRUE(somatic_mut)) {
    #Set up an gene annotation template to use
    mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
    annot <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol","chromosome_name","start_position","end_position"), mart=mart)
  }
  
  if(isTRUE(clinical)) {
    donor = paste0(url, project, "/", paste("donor", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("donor", project, "tsv.gz", sep="."))
    download.file(donor, .file)
    R.utils::gunzip(.file)
    donor=read_tsv(gsub(".gz$","", .file))
    
    sample = paste0(url, project, "/", paste("sample", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("sample", project, "tsv.gz", sep="."))
    download.file(sample, .file)
    R.utils::gunzip(.file)
    sample=read_tsv(gsub(".gz$","", .file))
    
    specimen = paste0(url, project, "/", paste("specimen", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("specimen", project, "tsv.gz", sep="."))
    download.file(specimen, .file)
    R.utils::gunzip(.file)
    specimen=read_tsv(gsub(".gz$","", .file))
    
    clinical_data <- merge.data.frame(donor, specimen, all = TRUE)
    clinical_data <- merge.data.frame(clinical_data, sample, all = TRUE)
    
    ## save clinical data
    data[["clinical"]] <- clinical_data
    
    compiled_dataset <- append(compiled_dataset, "clinical")
  }
  if(isTRUE(rna_seq)) {
    exp_seq = paste0(url, project, "/", paste("exp_seq", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("exp_seq", project, "tsv.gz", sep="."))
    download.file(exp_seq, .file)
    R.utils::gunzip(.file)
    exp_seq=read_tsv(gsub(".gz$","", .file))
    
    column_to_retain <- c("icgc_donor_id", "gene_id", "normalized_read_count", "raw_read_count")
    exp_seq <- exp_seq[, colnames(exp_seq) %in% column_to_retain]
    
    exp_seq_raw <- exp_seq %>%
      dplyr::select(-normalized_read_count) %>%
      pivot_wider(names_from = "icgc_donor_id", values_from = "raw_read_count", values_fn = {median}) %>%
      column_to_rownames("gene_id")
    
    exp_seq_normalized <- exp_seq %>%
      dplyr::select(-raw_read_count) %>%
      pivot_wider(names_from = "icgc_donor_id", values_from = "normalized_read_count", values_fn = {median}) %>%
      column_to_rownames("gene_id")
    
    ## save transcriptome data
    data[["rna_seq_raw"]] <- exp_seq_raw
    data[["rna_seq_normalized"]] <- exp_seq_normalized
    compiled_dataset <- append(compiled_dataset, "rna_seq")
  }
  if(isTRUE(methylation)) {
    meth_array = paste0(url, project, "/", paste("meth_array", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("meth_array", project, "tsv.gz", sep="."))
    download.file(meth_array, .file)
    R.utils::gunzip(.file)
    meth_array=read_tsv(gsub(".gz$","", .file))
    
    column_to_retain <- c("icgc_donor_id", "probe_id", "methylation_value")
    meth_array <- meth_array[, colnames(meth_array) %in% column_to_retain]
    
    ## convert probe_id to gene_id
    mh450 <- get450k()
    probenames <- unique(meth_array$probe_id)
    probes <- mh450[probenames]
    annot <- getNearestTranscript(probes)
    annot <- annot %>%
      dplyr::select(-queryHits) %>%
      rownames_to_column("probe_id")
    
    meth_array <- meth_array %>%
      pivot_wider(names_from = "icgc_donor_id", values_from = "methylation_value", values_fn = {median}) 
    
    ## merge_data
    meth_array <- annot %>%
      full_join(meth_array, by="probe_id")
    
    ## save transcriptome data
    data[["methylation"]] <- meth_array
    compiled_dataset <- append(compiled_dataset, "methylation")
  }
  
  if(isTRUE(mirna)) {
    mirna_seq = paste0(url, project, "/", paste("mirna_seq", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("mirna_seq", project, "tsv.gz", sep="."))
    download.file(mirna_seq, .file)
    R.utils::gunzip(.file)
    mirna_seq=read_tsv(gsub(".gz$","", .file))
    
    ## save transcriptome data
    data[["mirna_seq"]] <- mirna_seq
    compiled_dataset <- append(compiled_dataset, "mirna")
  }
  if(isTRUE(cnv)) {
    
    cnv_data = paste0(url, project, "/", paste("copy_number_somatic_mutation", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("copy_number_somatic_mutation", project, "tsv.gz", sep="."))
    download.file(cnv_data, .file)
    R.utils::gunzip(.file)
    cnv_data=read_tsv(gsub(".gz$","", .file))
    
    column_to_retain <- c("icgc_donor_id","chromosome","chromosome_start","chromosome_end","copy_number","segment_mean")
    cnv_data <- cnv_data[, colnames(cnv_data) %in% column_to_retain]
    colnames(cnv_data) <- c("icgc_donor_id","copy_number","segment_mean","Chr","Start","End")
    cnv_data <- cnv_data[, c("Chr","Start","End","copy_number","segment_mean", "icgc_donor_id")]
    
    ## annotate for gene_affected
    annot_gr <- annot[, c("hgnc_symbol","chromosome_name","start_position","end_position")]
    colnames(annot_gr) <- c("GeneSymbol","Chr","Start","End")
    annot_genomeRange <- makeGRangesFromDataFrame(annot_gr,keep.extra.columns = TRUE)
    
    cnv_data_genomeRange <- makeGRangesFromDataFrame(cnv_data, keep.extra.columns = TRUE)
    
    hits <- findOverlaps(annot_genomeRange, cnv_data_genomeRange, type="within")
    
    a <- as.data.frame(cnv_data[subjectHits(hits),])
    b <- as.data.frame(annot_gr[queryHits(hits),])
    
    cnv_data1 <- left_join(a,b,
                          by = c("Chr", "Start", "End"))
    
    ## save transcriptome data
    data[["cnv"]] <- cnv_data
    compiled_dataset <- append(compiled_dataset, "cnv")
  }
  if(isTRUE(somatic_mut)) {
    mut_data = paste0(url, project, "/", paste("simple_somatic_mutation.open", project, "tsv.gz", sep="."))
    .file = file.path(.exdir, paste("simple_somatic_mutation.open", project, "tsv.gz", sep="."))
    download.file(mut_data, .file)
    R.utils::gunzip(.file)
    mut_data=read_tsv(gsub(".gz$","", .file))
    
    column_to_retain <- c("icgc_donor_id","chromosome","chromosome_start","chromosome_end","mutation_type","mutated_from_allele","mutated_to_allele","platform","sequencing_strategy","gene_affected", "total_read_count", "mutant_allele_read_count")
    
    mut_data = mut_data[,colnames(mut_data) %in% column_to_retain]
    
    ## Ensembl ID to Gene Symbol
    annot1 <- annot %>%
      dplyr::select(ensembl_gene_id, hgnc_symbol)
    
    mut_data <- mut_data %>%
      left_join(annot1, by = c("gene_affected"="ensembl_gene_id")) %>%
      distinct()
    
    ## save transcriptome data
    data[["somatic_mutation"]] <- mut_data
    compiled_dataset <- append(compiled_dataset, "somatic_mut")
  }
  unlink(.exdir, recursive = TRUE)
  write_rds(data, path=paste0(project, "_", paste(compiled_dataset, collapse="_"), ".rds"))
  
}
