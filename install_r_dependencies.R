# Keep CRAN dependencies contemporaneous with R 4.4.2 and Bioconductor 3.20.
cran_snapshot <- "https://packagemanager.posit.co/cran/2024-10-31"
options(repos = c(CRAN = cran_snapshot))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

bioc_repos <- BiocManager::repositories(version = "3.20")
bioc_repos["CRAN"] <- cran_snapshot
options(repos = bioc_repos)

BiocManager::install(
    "Biostrings",
    version = "3.20",
    ask = FALSE,
    update = FALSE
)
if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Biostrings installation failed")
}
install.packages(c("iterators", "foreach"), dependencies = NA)
# fit1 was serialized with this version; do not allow the snapshot's newer
# randomForest release to change predictions made by the supplied models.
install.packages(
    "https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-14.tar.gz",
    repos = NULL,
    type = "source"
)
if (as.character(packageVersion("randomForest")) != "4.6.14") {
    stop("randomForest 4.6-14 installation failed")
}
