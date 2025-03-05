library("msigdb")

msigdb.hs = getMsigdb(org = 'hs', id = 'SYM', version = '7.4')
listCollections(msigdb.hs)

Hs.H <- list()
temp <- subsetCollection(msigdb.hs, 'h')
for (i in names(temp)) {
  name <- gsub("^HALLMARK_", "H: ", i)
  Hs.H[[name]] <- temp[[i]]@geneIds
}

Hs.c2 <- list()
temp <- subsetCollection(msigdb.hs, 'c2')
for (i in names(temp)) {
  name <- paste0("C2: ", i)
  Hs.c2[[name]] <- temp[[i]]@geneIds
}

Hs.c4 <- list()
temp <- subsetCollection(msigdb.hs, 'c4')
for (i in names(temp)) {
  name <- paste0("C4: ", i)
  Hs.c4[[name]] <- temp[[i]]@geneIds
}

Hs.c5bp <- list()
temp <- subsetCollection(msigdb.hs, collection = 'c5', subcollection = 'GO:BP')

# Filter names that start with "GOBP_" and process them all at once
valid_names <- grep("^GOBP_", names(temp), value = TRUE)
new_names <- gsub("^GOBP_", "BP: ", valid_names)

# Create the list directly without looping
Hs.c5bp <- setNames(
  lapply(valid_names, function(i) temp[[i]]@geneIds),
  new_names
)

Hs.c5mf <- list()
temp <- subsetCollection(msigdb.hs, collection = 'c5', subcollection = 'GO:MF')

# Filter names that start with "GOBP_" and process them all at once
valid_names <- grep("^GOMF_", names(temp), value = TRUE)
new_names <- gsub("^GOMF_", "MF: ", valid_names)

# Create the list directly without looping
Hs.c5mf <- setNames(
  lapply(valid_names, function(i) temp[[i]]@geneIds),
  new_names
)

Hs.c6 <- list()
temp <- subsetCollection(msigdb.hs, 'c6')
for (i in names(temp)) {
  name <- paste0("C6: ", i)
  Hs.c6[[name]] <- temp[[i]]@geneIds
}

Hs.c7 <- list()
temp <- subsetCollection(msigdb.hs, 'c7')
for (i in names(temp)) {
  name <- paste0("C7: ", i)
  Hs.c7[[name]] <- temp[[i]]@geneIds
}

Hs.c8 <- list()
temp <- subsetCollection(msigdb.hs, 'c8')
for (i in names(temp)) {
  name <- paste0("C8: ", i)
  Hs.c8[[name]] <- temp[[i]]@geneIds
}

save(Hs.H, Hs.c2, Hs.c4, Hs.c5bp, Hs.c5mf,Hs.c6, Hs.c7, Hs.c8, file = "data/Human_msigdb.RData")
