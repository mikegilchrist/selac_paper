combine_all_files <- function() {
  file.names <- system("ls -1  ../SELAC_DATA_FILES/DATA/*fasta", intern=TRUE)
  command.string <- c("phyutility -concat -in ")
  for (i in seq_along(file.names)) {
    command.string <- paste(command.string, file.names[i])
  }
  command.string <- paste(command.string, "-out combined.nex")
  cat(command.string)
  system(command.string, intern=TRUE)
  return("combined.nex")
}

convert_to_fasta <- function(combined) {
  seqs <- ape::read.nexus.data(combined)
  ape::write.dna(seqs, format="fasta", file="combined.fasta")
  return("combined.fasta")
}

run_IQtree <- function(fasta) {
  iqrun <- system(paste0("iqtree -s ", fasta, " -st CODON -nt AUTO"), intern=TRUE)
  return(iqrun)
}
