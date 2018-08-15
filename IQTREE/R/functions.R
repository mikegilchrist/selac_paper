combine_all_files <- function() {
  file.names <- system("ls -1  ../SELAC_DATA_FILES/DATA/*fasta", intern=TRUE)
  command.string <- c("phyutility -concat -in ")
  for (i in seq_along(file.names)) {
    command.string <- paste(command.string, file.names[i])
  }
  command.string <- paste(command.string, "-out combined.nex")
  cat(command.string)
  return(system(command.string, intern=TRUE))
}

convert_to_fasta <- function() {
  seqs <- ape::read.nexus.data("combined.nex")
  return(ape::write.FASTA(seqs, file="combined.fasta"))
}
