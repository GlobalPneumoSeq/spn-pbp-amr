options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install(
    "Biostrings",
    version = "3.20",
    ask = FALSE,
    update = FALSE
)
install.packages(c("randomForest", "iterators", "foreach"), dependencies = NA)
