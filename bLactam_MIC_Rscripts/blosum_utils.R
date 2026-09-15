choose_training_residue <- function(residue, valid_training, blosum) {
  score_source <- substr(residue, 1, 1)
  if (!(score_source %in% blosum$src)) {
    score_source <- "X"
  }

  scores <- as.numeric(blosum[blosum$src == score_source, valid_training, drop = TRUE])
  names(scores) <- valid_training
  if (!length(scores) || anyNA(scores)) {
    stop("BLOSUM62 does not contain scores for the requested training residues")
  }

  names(scores)[which.max(scores)]
}
