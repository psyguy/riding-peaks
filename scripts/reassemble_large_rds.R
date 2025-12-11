# Reassemble split RDS files
#
# Large RDS files (>100MB) are split into parts to comply with GitHub's file size limit.
# This script reassembles them into the original files.
#
# Usage: source("scripts/reassemble_large_rds.R")

library(here)

reassemble_rds <- function(base_path, n_parts = 3) {
  # Read all parts
  parts <- lapply(1:n_parts, function(i) {
    part_file <- paste0(base_path, "_part", i, ".rds.split")
    if (!file.exists(part_file)) {
      stop("Missing part file: ", part_file)
    }
    readBin(part_file, "raw", file.info(part_file)$size)
  })

  # Combine parts
  combined <- do.call(c, parts)

  # Write to temp file and read as RDS
  tmp <- tempfile(fileext = ".rds")
  writeBin(combined, tmp)
  obj <- readRDS(tmp)
  unlink(tmp)

  return(obj)
}

# Reassemble brms_pa_n1_strategy1_fits.rds if not present
target_file <- here("fits", "brms_pa_n1_strategy1_fits.rds")
if (!file.exists(target_file)) {
  cat("Reassembling brms_pa_n1_strategy1_fits.rds from split files...\n")

  base_path <- here("fits", "brms_pa_n1_strategy1_fits")
  obj <- reassemble_rds(base_path, n_parts = 3)

  # Optionally save the reassembled file locally
  saveRDS(obj, target_file)
  cat("Saved reassembled file to:", target_file, "\n")
} else {
  cat("File already exists:", target_file, "\n")
}

# Function to load directly without saving (saves disk space)
load_split_rds <- function(base_name) {
  base_path <- here("fits", base_name)
  reassemble_rds(base_path, n_parts = 3)
}

cat("\nTo load the N=1 brms fits directly:\n")
cat("  fits_n1 <- load_split_rds('brms_pa_n1_strategy1_fits')\n")
