EMT_signature_score <- function(norm_matrix, organism = c("hsapiens", "mmusculus")) {
  organism <- match.arg(organism, c("hsapiens", "mmusculus"))
  
  ## download EMT_signature
  EMT_signature_genes <- read.csv("https://raw.githubusercontent.com/umahajanatlmu/useful_commands/main/transcriptome/EMT_signatures.csv", sep = ",", na ="")
  
  ## convert to mouse
  if (organism == "mmusculus") {
    for (i in colnames(EMT_signature_genes)) {
      EMT_signature_genes[[i]] <- .convertHumanGeneList(EMT_signature_genes[[i]])
    }
  }
  
  # Groger signature
  groger_EMT_score <- apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Groger_M)), ], 2, median) - 
    apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Groger_E)), ], 2, median)
  # Creighton signature
  creighon_EMT_score <- apply(apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Creighton_M)), ], 1, .z_transformation), 1, sum, na.rm = TRUE) - 
    apply(apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Creighton_E)), ], 1, .z_transformation), 1, sum, na.rm = TRUE)
  # Byers signature
  byers_EMT_score <- apply(apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Byers_M)), ], 1, .z_transformation), 1, sum, na.rm = TRUE)/length(unique(na.omit(EMT_signature_genes$Creighton_M))) - 
    apply(apply(norm_matrix[rownames(norm_matrix) %in% unique(na.omit(EMT_signature_genes$Creighton_E)), ], 1, .z_transformation), 1, sum, na.rm = TRUE)/length(unique(na.omit(EMT_signature_genes$Creighton_E)))
  # Combining two EMT score.
  EMT_score <- data.frame(Groger_EMT_score = groger_EMT_score[colnames(norm_matrix)], Creighton_EMT_score = creighon_EMT_score[colnames(norm_matrix)],
                          Byers_EMT_score = byers_EMT_score[colnames(norm_matrix)], row.names = colnames(norm_matrix))
  
  return(EMT_score)
  
  
}

.z_transformation <- function(x) (x - median(x))/sd(x)

.convertHumanGeneList <- function(gene_list){
  mouse_human_genes =
    read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  output = c()
  
  for(gene in gene_list){
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      mouse_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
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