TS_signature_score <- function(norm_matrix, organism = c("hsapiens", "mmusculus")) {
  organism <- match.arg(organism, c("hsapiens", "mmusculus"))
  
  ## download EMT_signature
  TS_signature_genes <- read.csv("https://raw.githubusercontent.com/umahajanatlmu/useful_commands/main/transcriptome/TS_score.csv", sep = ";", na ="")
  
  ## convert to mouse
  if (organism == "mmusculus") {
      TS_signature_genes$TS_genes_mmusculus <- .convertHumanGeneList(TS_signature_genes$TS_genes)
      # Extract the expression of the 61 TS-related genes
      TS_related_expression <- norm_matrix[rownames(norm_matrix) %in% TS_signature_genes$TS_genes_mmusculus,]
  } else
    # Extract the expression of the 61 TS-related genes
    TS_related_expression <- norm_matrix[rownames(norm_matrix) %in% TS_signature_genes$TS_genes, ]
  
  # Perform Principal Component Analysis (PCA) based on TCGA BLCA data
  pca_result <- prcomp(t(TS_related_expression))
  
  # Calculate the TS score using PC1 and PC2
  pca_result <- as.data.frame(pca_result$x)
  
  TS_score <- pca_result %>%
    dplyr::select(PC1, PC2) %>%
    mutate(TS_score = PC2-PC1) %>%
    dplyr::select(TS_score)
  
  return(TS_score)

}

.convertHumanGeneList <- function(gene_list){
  mouse_human_genes =
    read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% dplyr::filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% dplyr::filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      if (length(mouse_genes) != 1) {
        mouse_genes <- mouse_genes[1]
      }
      for(mouse_gene in mouse_genes){
        output = append(output,mouse_gene)
      }
    } else
      output = append(output, NA)
  }
  
  
  return (output)
}
