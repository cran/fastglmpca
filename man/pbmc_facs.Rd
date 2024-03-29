% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbmc_facs.R
\docType{data}
\name{pbmc_facs}
\alias{pbmc_facs}
\title{Mixture of 10 FACS-purified PBMC Single-Cell RNA-seq data}
\format{
\code{pbmc_facs} is a list with the following elements:

\describe{

  \item{counts}{16,791 x 3,774 sparse matrix of UMI counts, with
     rows corresponding to genes and columns corresponding to
     cells (samples). It is an object of class \code{"dgCMatrix"}).}

  \item{counts_test}{UMI counts for an additional test set of 100
    cells.}

  \item{samples}{Data frame containing information about the
    samples, including cell barcode and source FACS population
    (\dQuote{celltype} and \dQuote{facs_subpop}).}

  \item{samples_test}{Sample information for the additional test
     set of 100 cells.}

  \item{genes}{Data frame containing information and the genes,
    including gene symbol and Ensembl identifier.}

  \item{fit}{GLM-PCA model that was fit to the UMI count data in
    the vignette.}}
}
\source{
\url{https://www.10xgenomics.com/resources/datasets}
}
\description{
These data are a selection of the reference
  transcriptome profiles generated via single-cell RNA sequencing
  (RNA-seq) of 10 bead-enriched subpopulations of PBMCs (Donor A),
  described in Zheng \emph{et al} (2017). The data are unique
  molecular identifier (UMI) counts for 16,791 genes in 3,774 cells.
  (Genes with no expression in any of the cells were removed.) Since
  the majority of the UMI counts are zero, they are efficiently
  stored as a 16,791 x 3774 sparse matrix. These data are used in the
  vignette illustrating how \sQuote{fastglmpca} can be used to
  analyze single-cell RNA-seq data. Data for a separate set of 1,000
  cells is provided as a \dQuote{test set} to evaluate out-of-sample
  predictions.
}
\examples{
library(Matrix)
data(pbmc_facs)
cat(sprintf("Number of genes: \%d\n",nrow(pbmc_facs$counts)))
cat(sprintf("Number of cells: \%d\n",ncol(pbmc_facs$counts)))
cat(sprintf("Proportion of counts that are non-zero: \%0.1f\%\%.\n",
            100*mean(pbmc_facs$counts > 0)))

}
\references{
G. X. Y. Zheng \emph{et al} (2017). Massively parallel digital
  transcriptional profiling of single cells. \emph{Nature Communications}
  \bold{8}, 14049. \doi{10.1038/ncomms14049}
}
\keyword{data}
