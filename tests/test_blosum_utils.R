args <- commandArgs(trailingOnly = TRUE)
script_dir <- if (length(args) > 0) {
  normalizePath(file.path(dirname(args[1]), ".."), mustWork = TRUE)
} else {
  normalizePath(getwd(), mustWork = TRUE)
}
source(file.path(script_dir, "bLactam_MIC_Rscripts", "blosum_utils.R"))

blosum <- read.csv(
  file.path(script_dir, "bLactam_MIC_Rscripts", "BLOSUM62.csv"),
  colClasses = "character",
  check.names = FALSE
)

stopifnot(choose_training_residue("A", c("A", "R"), blosum) == "A")
stopifnot(choose_training_residue("J", c("A", "V"), blosum) == "A")
stopifnot(choose_training_residue("B", c("D", "N"), blosum) == "D")
stopifnot(choose_training_residue("-", c("A", "-"), blosum) == "-")

tie <- data.frame(src = "A", A = "4", R = "4", check.names = FALSE)
stopifnot(choose_training_residue("A", c("A", "R"), tie) == "A")

cat("BLOSUM substitution tests passed\n")
