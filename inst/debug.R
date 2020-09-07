# This is a script to aid debugging.
# Copy the config from the stdout of running MotifBinner2.R below.

library(MotifBinner2)
stop_at_op <- 'n031'

config <- list(operation_list = list(n001 = list(name = "fwd_loadData", 
    op = "loadData", data_source = "/tmp/left.fastq", max_seq = 0, 
    cache_data = TRUE), n002 = list(name = "fwd_basicQC", op = "basicQC", 
    data_source = "n001", cache_data = FALSE), n003 = list(name = "fwd_ambigSeqs", 
    op = "ambigSeqs", data_source = "n001", threshold = 0.02, 
    cache_data = TRUE), n004 = list(name = "fwd_primerDimer", 
    op = "primerDimer", data_source = "n003", threshold = 80, 
    cache_data = TRUE), n005 = list(name = "fwd_seqLength", op = "seqLength", 
    data_source = "n004", threshold = 295, cache_data = TRUE), 
    n006 = list(name = "fwd_qualTrim", op = "qualTrim", data_source = "n005", 
        avg_qual = 20, bad_base_threshold = 10, max_bad_bases = 0.15, 
        cache_data = TRUE), n007 = list(name = "fwd_trimAffixes", 
        op = "trimAffixes", data_source = "n006", primer_seq = "GAGGAGATATGNNNNNNNNAGGGACAATTG", 
        primer_lens = c(11, 8, 11), min_score = 25, primer_location = "front", 
        front_gaps_allowed = 0, cache_data = TRUE), n008 = list(
        name = "rev_loadData", op = "loadData", data_source = "/tmp/right.fastq", 
        max_seq = 0, cache_data = TRUE), n009 = list(name = "rev_basicQC", 
        op = "basicQC", data_source = "n008", cache_data = FALSE), 
    n010 = list(name = "rev_ambigSeqs", op = "ambigSeqs", data_source = "n008", 
        threshold = 0.02, cache_data = TRUE), n011 = list(name = "rev_primerDimer", 
        op = "primerDimer", data_source = "n010", threshold = 80, 
        cache_data = TRUE), n012 = list(name = "rev_seqLength", 
        op = "seqLength", data_source = "n011", threshold = 295, 
        cache_data = TRUE), n013 = list(name = "rev_qualTrim", 
        op = "qualTrim", data_source = "n012", avg_qual = 20, 
        bad_base_threshold = 10, max_bad_bases = 0.15, cache_data = TRUE), 
    n014 = list(name = "rev_trimAffixes", op = "trimAffixes", 
        data_source = "n013", primer_seq = "AAGTGAAGAGGTTAACAGGGA", 
        primer_lens = 21, min_score = 21, primer_location = "front", 
        front_gaps_allowed = 0, cache_data = TRUE), n015 = list(
        name = "fwd_extractPIDs", op = "extractPIDs", data_source = "n007", 
        pid_in_which_fragment = 2, pattern_to_chop_from_names = " [0-9]:N:[0-9]*:[0-9]*$", 
        pid_gaps_allowed = 0, cache_data = TRUE), n016 = list(
        name = "rev_extractPIDs", op = "extractPIDs", data_source = "n014", 
        pid_in_which_fragment = NULL, pattern_to_chop_from_names = " [0-9]:N:[0-9]*:[0-9]*$", 
        pid_gaps_allowed = 0, cache_data = TRUE), n017 = list(
        name = "matchPairs", op = "matchPairs", data_source = c(fwd = "n015", 
        rev = "n016"), cache_data = TRUE), n018 = list(name = "processBadPIDs", 
        op = "processBadPIDs", data_source = "n017", cache_data = TRUE), 
    n019 = list(name = "fwd_extractReads", op = "extractData", 
        data_source = "n018", extract_levels = c("seq_dat", "fwd"
        ), cache_data = TRUE), n020 = list(name = "fwd_alignBinsMSA", 
        op = "alignBinsMSA", bins_to_process = Inf, data_source = "n019", 
        cache_data = TRUE), n021 = list(name = "fwd_buildConsensus", 
        op = "buildConsensus", data_source = "n020", cache_data = TRUE), 
    n022 = list(name = "fwd_removeGaps", op = "removeChars", 
        data_source = "n021", char_to_remove = "-", cache_data = TRUE), 
    n023 = list(name = "rev_extractReads", op = "extractData", 
        data_source = "n018", extract_levels = c("seq_dat", "rev"
        ), cache_data = TRUE), n024 = list(name = "rev_alignBinsMSA", 
        op = "alignBinsMSA", bins_to_process = Inf, data_source = "n023", 
        cache_data = TRUE), n025 = list(name = "rev_buildConsensus", 
        op = "buildConsensus", data_source = "n024", cache_data = TRUE), 
    n026 = list(name = "rev_removeGaps", op = "removeChars", 
        data_source = "n025", char_to_remove = "-", cache_data = TRUE), 
    n030 = list(name = "primerSeqErr", op = "primerSeqErr", data_source = c(fwd = "n007", 
    rev = "n014"), cache_data = FALSE), n031 = list(name = "binSeqErr", 
        op = "binSeqErr", data_source = c(bin_msa_fwd = "n020", 
        bin_msa_rev = "n024", cons_fwd = "n021", cons_rev = "n025", 
        primer_err = "n030"), cache_data = FALSE), n100 = list(
        name = "dataTracing", op = "dataTracing", data_source = c(fwdReads.01 = "n001", 
        fwdReads.02 = "n003", fwdReads.03 = "n004", fwdReads.04 = "n005", 
        fwdReads.05 = "n006", fwdReads.06 = "n007", fwdReads.07 = "n017", 
        fwdReads.08 = "n018", fwdReads.09 = "n020", fwdReads.10 = "n021", 
        revReads.01 = "n008", revReads.02 = "n010", revReads.03 = "n011", 
        revReads.04 = "n012", revReads.05 = "n013", revReads.06 = "n014", 
        revReads.07 = "n017", revReads.08 = "n018", revReads.09 = "n024", 
        revReads.10 = "n025"), cache_data = FALSE), n101 = list(
        name = "prepConfig", op = "prepConfig", cache_data = FALSE)), 
    output_dir = "/tmp/mb2_test/", base_for_names = "mb2_test", 
    header_format = "SRA", intermediate_reports = TRUE, verbosity = 3, 
    erase_history = FALSE, report_type = "html", ncpu = "2")

unlink(file.path(config$output_dir, config$base_for_names), recursive=T)
all_results <- list()
class(all_results) <- 'allResults'

for (op_num in names(config$operation_list)){
  if (op_num == stop_at_op){
    stop('Reached stopping point')
  }
  all_results <-  applyOperation(all_results, config, op_number = op_num)
}

print(names(all_results))
print(config$operation_list[[stop_at_op]])
config$current_op_number <- stop_at_op

print('now go to the code for the following function and step for it:')
print(config$operation_list[[stop_at_op]]$op)


#file.copy(from = file.path(all_results$n027_cons_ambigSeqs$config$op_dir, 
#                           paste(config$base_for_names, "kept", "cons_ambigSeqs.fastq", sep = "_")),
#          to = file.path(config$output_dir, config$base_for_names))
#
#x <- genReport(all_results, config)
#
