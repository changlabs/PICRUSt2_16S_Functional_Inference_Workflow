################################################################################
# Script: decompress_output_files_function.R
# Purpose:
#   - Define reusable functions for recursively finding and decompressing
#     every gzip (.gz) or zip (.zip) file under a given directory, deleting
#     each compressed original once its decompression is confirmed to have
#     succeeded.
#   - Shared by workflow steps that need plain copies of PICRUSt2's compressed
#     tabular outputs (for example, `pred_metagenome_unstrat.tsv.gz`). Step 2
#     uses this helper for its community-level outputs; Step 5 deliberately
#     keeps its much larger contribution tables compressed.
#
# Usage:
#   source(here("R", "functions", "decompress_output_files_function.R"))
#   decompression_summary <- decompress_all_under(some_output_directory)
#
# Exposes three functions:
#   - decompress_gz_file(gz_path, chunk_size)  : decompress one .gz file
#   - decompress_zip_file(zip_path)            : extract one .zip archive
#   - decompress_all_under(search_dir)         : recursively decompress
#                                                 every .gz/.zip file found
#                                                 under search_dir and
#                                                 return a summary data.frame
################################################################################


# ==============================================================================
# Function: decompress_gz_file
# ==============================================================================
# Decompress a Single Gzip File and Report Whether It Succeeded
#
# Streams a .gz file to its decompressed form in fixed-size chunks (rather
# than reading it into memory in one call), so this works regardless of how
# large the decompressed content is. The output path is the input path with
# its .gz extension removed.
#
decompress_gz_file <- function(gz_path, chunk_size = 1e7) {
  target_path <- sub("\\.gz$", "", gz_path, ignore.case = TRUE)

  result <- tryCatch({
    con_in <- gzfile(gz_path, open = "rb")
    on.exit(close(con_in), add = TRUE)

    con_out <- file(target_path, open = "wb")
    on.exit(close(con_out), add = TRUE)

    repeat {
      chunk <- readBin(con_in, what = "raw", n = chunk_size)
      if (length(chunk) == 0) break
      writeBin(chunk, con_out)
    }

    TRUE
  }, error = function(e) {
    warning("Failed to decompress ", gz_path, ": ", conditionMessage(e))
    # Remove a possibly-partial output file so a failed attempt doesn't leave
    # behind a truncated, misleading decompressed file.
    if (file.exists(target_path)) unlink(target_path)
    FALSE
  })

  list(success = isTRUE(result), target_path = target_path)
}


# ==============================================================================
# Function: decompress_zip_file
# ==============================================================================
# Extract a Single Zip Archive and Report Whether It Succeeded
#
# PICRUSt2 itself does not produce .zip archives (only gzip-compressed
# single files), but this is included for completeness in case any tool
# downstream of it in a pipeline using this function does. All contained
# files are extracted into the archive's own directory.
#
decompress_zip_file <- function(zip_path) {
  target_dir <- dirname(zip_path)

  result <- tryCatch({
    utils::unzip(zip_path, exdir = target_dir)
    TRUE
  }, error = function(e) {
    warning("Failed to extract ", zip_path, ": ", conditionMessage(e))
    FALSE
  })

  list(success = isTRUE(result), target_path = target_dir)
}


# ==============================================================================
# Function: decompress_all_under
# ==============================================================================
# Recursively Decompress Every .gz/.zip File Under a Directory
#
# Finds every compressed file under `search_dir`, decompresses each one,
# and deletes the original compressed file once decompression succeeds
# (files that fail to decompress are left in place so no data is lost).
#
decompress_all_under <- function(search_dir) {
  compressed_paths <- list.files(
    search_dir, pattern = "\\.(gz|zip)$",
    recursive = TRUE, full.names = TRUE, ignore.case = TRUE
  )

  if (length(compressed_paths) == 0) {
    return(data.frame(
      `Original File` = character(0), `Decompressed To` = character(0),
      `Original Removed` = logical(0), check.names = FALSE
    ))
  }

  summary_rows <- lapply(compressed_paths, function(path) {
    is_zip <- grepl("\\.zip$", path, ignore.case = TRUE)
    outcome <- if (is_zip) decompress_zip_file(path) else decompress_gz_file(path)

    # Only remove the original once decompression is confirmed successful.
    removed <- FALSE
    if (outcome$success) {
      removed <- isTRUE(unlink(path) == 0)
    }

    data.frame(
      `Original File` = path,
      `Decompressed To` = outcome$target_path,
      `Original Removed` = removed,
      check.names = FALSE
    )
  })

  do.call(rbind, summary_rows)
}
