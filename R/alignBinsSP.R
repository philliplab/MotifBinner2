#' Removes sequences shorter than a given cutoff
#' @inheritParams applyOperation
#' @export

alignBinsSP <- function(all_results, config)
{
  op_number <- config$current_op_number
  op_args <- config$operation_list[[op_number]]
  op_full_name <- paste(op_number, op_args$name, sep = '_')

  op_dir <- file.path(config$output_dir, config$base_for_names, op_full_name)
  dir.create(op_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(op_dir, 'bins'), showWarnings = FALSE, recursive = TRUE)
  
  which_pair <- op_args$which_pair
  data_source_indx <- grep(op_args$data_source, names(all_results))

  stopifnot(length(data_source_indx) == 1)
  stopifnot(all(sort(names(all_results[[data_source_indx]]$seq_dat)) == c("fwd", "rev")))
  stopifnot(which_pair %in% c('fwd', 'rev'))
  seq_dat <- all_results[[data_source_indx]]$seq_dat[[which_pair]]

  bins_to_process <- op_args$bins_to_process

  per_read_metrics <- data.frame(read_name = as.character(seq_dat@id),
                                 stringsAsFactors = F)
  per_read_metrics$pid <- gsub("^.*_PID:" , "", per_read_metrics$read_name)
  per_read_metrics$clean_pid <- gsub("_" , "", per_read_metrics$pid)

  profile_seqs <- readDNAStringSet(op_args$profile_file)

  registerDoMC(cores = config$ncpu)
  uniq_pids <- unique(per_read_metrics$clean_pid)
  if (is.null(bins_to_process)){
    bins_to_process <- length(uniq_pids)
  } else {
    bins_to_process <- min(length(uniq_pids), bins_to_process)
  }
  pid <- uniq_pids[1]
  tmp_x <- foreach(pid = uniq_pids[1:bins_to_process], .combine = "c") %dopar% {
    cur_seq_indxs <- which(per_read_metrics$clean_pid == pid)
    cur_seq_names <- per_read_metrics$read_name[cur_seq_indxs]
    stopifnot(all(cur_seq_names == as.character(seq_dat@id)[cur_seq_indxs]))

    cur_seqs <- seq_dat[cur_seq_indxs]

    aligned_with_qual <- alignBinsSP_internal(cur_seqs, profile_seqs, op_dir, pid, which_pair)
    aligned_with_qual
  }
  all_bins_aligned_with_qual <- shortReadQ_forced_append(tmp_x)
  rm(tmp_x)

  per_read_metrics <- data.frame('read_exists' = rep(1, length(all_bins_aligned_with_qual)))
  trim_steps <- list(step1 = list(name = 'read_exists',
                                  threshold = 1,
                                  breaks = c(1)))

  result <- list(trim_steps = trim_steps,
                 metrics = list(per_read_metrics = per_read_metrics))
  class(result) <- 'alignBinsSP'
  if (op_args$cache){
    result$seq_dat <- all_bins_aligned_with_qual
  }
  result$input_dat <- all_bins_aligned_with_qual
  result$config <- list(op_number = op_number,
                        op_args = op_args,
                        op_full_name = op_full_name,
                        op_dir = op_dir)
  return(result)
}

#' aligns fwd and rev sequences to a profile
#'
#' a custom guide tree is constructed to ensure that the alignment of the fwd
#' reads does not interfere with the alignment of the rev reads and vice versa.
#' @export

alignBinsSP_internal <- function(cur_seqs, profile_seqs, working_dir, pid, which_pair)
{
  interleaving_vector <- 1:length(cur_seqs)

  if (which_pair == 'fwd'){
    interleaved_seqs <- cur_seqs@sread
  } else if (which_pair == 'rev'){
    interleaved_seqs <- reverseComplement(cur_seqs@sread)
  } else {
    stop('which_pair must be either fwd or rev')
  }
  names(interleaved_seqs) <- paste(as.character(cur_seqs@id), which_pair, sep = '_')
  interleaved_seqs <- c(profile_seqs, interleaved_seqs)
  mafft_guide_tree <- data.frame(r1 = 1,
                                 r2 = 2:length(interleaved_seqs),
                                 r3 = 0.01,
                                 r4 = 0.01)
  gt_file_name <- file.path(working_dir, 'bins', paste(pid, '_guide_tree.txt', sep = ''))
  write.table(mafft_guide_tree,
              gt_file_name,
              col.names = FALSE, row.names = FALSE, sep = '\t')
  interleaved_file_name <- file.path(working_dir, 'bins', paste(pid, '_interleaved', '.fasta', sep = ''))
  writeXStringSet(interleaved_seqs,
                  interleaved_file_name,
                  width=20000)
  aligned_file_name <- file.path(working_dir, 'bins', paste(pid, '_aligned', '.fasta', sep = ''))
  system(paste('mafft --quiet --retree 1 --treein ', gt_file_name, ' ', interleaved_file_name, 
               ' > ', aligned_file_name, sep = ''))
  stopifnot(file.exists(aligned_file_name))
  aligned_seqs <- readDNAStringSet(aligned_file_name)
  aligned_seqs <- aligned_seqs[!(names(aligned_seqs) %in% names(profile_seqs))]
  
  if (which_pair == 'fwd'){
    interleaved_quals <- cur_seqs@quality@quality
  } else {
    interleaved_quals <- reverse(cur_seqs@quality@quality)
  } 
  names(interleaved_quals) <- paste(as.character(cur_seqs@id), which_pair, sep = '_')
  
  gap_only_cols_cpp_indexing <-
  which(consensusMatrix(aligned_seqs)['-',] == length(aligned_seqs)) - 1
  
  reads_and_qual <- transfer_gaps_cpp(as.character(aligned_seqs),
                                      as.character(interleaved_quals), 
                                      gap_only_cols_cpp_indexing)
  
  qual_mat <- as(FastqQuality(reads_and_qual$quals), 'matrix')
  avg_quals <- trunc(apply(qual_mat, 1, (function(x) {mean(x[x>0])})))
  tweaked_qual_mat <- gapQualityTweaker_non_ol_cpp(reads_and_qual$reads, qual_mat, which_pair, avg_quals)

  aligned_with_qual <-
  ShortReadQ(sread = DNAStringSet(tweaked_qual_mat$reads),
             quality = BStringSet(tweaked_qual_mat$quals),
             id = BStringSet(names(aligned_seqs)))
  return(aligned_with_qual)
}

saveToDisk.alignBinsSP <- function(result, config, seq_dat)
{
  kept <- getKept(result, seq_dat)
  trimmed <- getTrimmed(seq_dat = seq_dat, kept_dat = kept)

  if (length(kept) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_kept_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(kept, tmp_name, compress=F)
  }
  if (length(trimmed) > 0)
  {
    tmp_name <- file.path(result$config$op_dir, 
      paste(config$base_for_names, '_trimmed_', result$config$op_args$name, '.fastq', sep = ''))
    writeFastq(trimmed, tmp_name, compress=F)
  }
  return(result)
}

computeMetrics.alignBinsSP <- function(result, config, seq_dat)
{
  return(result)
}

print.alignBinsSP <- function(result, config)
{
  cat('\n-------------------')
  cat('\nOperation: alignBinsSP')
  cat('\n-------------------')
  cat('\nKept Sequences:\n')
  print(result$summary[,c('parameter', 'k_seqs', 'k_mean_length', 'k_mean_qual')])
  cat('\n-------------------')
  cat('\nTrimmed Sequences:\n')
  print(result$summary[,c('parameter', 't_seqs', 't_mean_length', 't_mean_qual')])
  invisible(result)
}


#  result <- operation_function(all_results, config)
#  result <- saveToDisk(result, config)
#  result <- genReport(result, config)
#  result <- genSummary(result, config)
#  result <- print(result, config)
