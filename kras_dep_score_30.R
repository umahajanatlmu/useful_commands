kras_dep_score_30 <- function (se) {
  kds30 <- read.csv("https://raw.githubusercontent.com/umahajanatlmu/useful_commands/main/kds30_pmid_36852277.csv", 
                          sep = ";")
  kds30$KDS30 <- scan(text=kds30$KDS30, dec=",", sep=".")
  kds30_genes <- unlist(kds30$Gene.ID)
  
  mat <- log2(assay(se)+1)
  ## se_to_mat
  mat_sub <- mat[rownames(mat) %in% kds30_genes,]
  
  mat_sub <- heatmaply::percentize(mat_sub)
  
  kds30 <- kds30[kds30$Gene.ID %in% rownames(mat_sub),]
  
  kds30_r <- data.frame()

  for ( i in colnames(mat_sub)) {
    pcc <- cor.test(mat_sub[[i]], kds30$KDS30)
    kds30_r <- rbind(kds30_r, data.frame(samples=i,
                                 kds30=pcc$estimate))
  }
  
  return (kds30_r)
}
