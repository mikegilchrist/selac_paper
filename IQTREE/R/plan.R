testmodels <- drake_plan(
  combined_worked = combine_all_files(),
  fasta_conversion = convert_to_fasta(),
  test_run = system("iqtree -s combined.fasta -st CODON")
)
