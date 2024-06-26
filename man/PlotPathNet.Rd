% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BioM2.R
\name{PlotPathNet}
\alias{PlotPathNet}
\title{Network diagram of pathways-level features}
\usage{
PlotPathNet(
  data = NULL,
  BioM2_pathways_obj = NULL,
  FeatureAnno = NULL,
  pathlistDB = NULL,
  PathNames = NULL,
  cutoff = 0.2,
  num = 10
)
}
\arguments{
\item{data}{The input omics data}

\item{BioM2_pathways_obj}{Results produced by BioM2()}

\item{FeatureAnno}{The annotation data stored in a data.frame for probe mapping.
It must have at least two columns named 'ID' and 'entrezID'.
(For details, please refer to data( data("MethylAnno") )}

\item{pathlistDB}{A list of pathways with pathway IDs and their corresponding genes ('entrezID' is used).
For details, please refer to ( data("GO2ALLEGS_BP") )}

\item{PathNames}{A vector.A vector containing the names of pathways}

\item{cutoff}{Threshold for correlation between features within a pathway}

\item{num}{The first few internal features of each pathway that are most relevant to the phenotype}
}
\value{
a ggplot object
}
\description{
Network diagram of pathways-level features
}
