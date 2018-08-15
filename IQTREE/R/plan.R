testmodels <- drake_plan(
  combined_file = combine_all_files(),
  fasta_conversion = convert_to_fasta(combined_file),
  test_run = run_IQtree(fasta_conversion),
  save(test_run, file_out("iqresult.rda"))
)
