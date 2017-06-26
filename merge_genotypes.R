# Helper functions

read.genotypes <- function(pheno.dict, sep.char = ",") {
  geno.data <- read.table(pheno.dict, sep = sep.char, header = F, stringsAsFactors = F)
  genotypes <- geno.data[[1]]
  full.genotypes.list <- lapply(rownames(geno.data), function(i) sapply(strsplit(geno.data[i,2], split = ',')[[1]], trimws))
  full.genotypes.list <- lapply(full.genotypes.list, function(x) make.names(x))
  names(full.genotypes.list) <- genotypes
  return(full.genotypes.list)
}


write.genotypes <- function(genotypes.list, out.file) {
  geno.data <- sapply(genotypes.list, function(x) paste('\"', paste(x, collapse = ", "), '\"', sep = ""))
  geno.data <- sapply(names(geno.data), function(x) paste(x, geno.data[[x]], sep = ","))
  fileConn = file(out.file)
  writeLines(geno.data, fileConn)
  close(fileConn)
}


# Command line arguments
options <- commandArgs(trailingOnly = T)
out.file <- options[[1]]

combined.genotypes.list <- list()
for (i in 2:length(options)) {
  j <- i - 1
  genotypes.list <- read.genotypes(options[[j]])
  genotypes.list <- lapply(genotypes.list, function(cells) sapply(cells, function(x) paste(x, j, sep = "-")))
  for (g in names(genotypes.list)) {
    if (g %in% names(combined.genotypes.list)) {
      combined.genotypes.list[[g]] <- c(combined.genotypes.list[[g]], genotypes.list[[g]])
    }
    else {
      combined.genotypes.list[[g]] <- genotypes.list[[g]]
    }
  }
}

print(sapply(combined.genotypes.list, length))
print(paste("Total cells", sum(sapply(combined.genotypes.list, length))))
write.genotypes(combined.genotypes.list, out.file)
