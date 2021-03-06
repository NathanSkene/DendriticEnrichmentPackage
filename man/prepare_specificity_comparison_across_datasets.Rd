% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/prepare_specificity_comparison_across_datasets.r
\name{prepare_specificity_comparison_across_datasets}
\alias{prepare_specificity_comparison_across_datasets}
\title{compare.datasets}
\usage{
prepare_specificity_comparison_across_datasets(ctd, keepCells, species,
  datasetName = NA, useCell, annotLevel = 1,
  sharedName = "sharedCell", datasetGroup)
}
\arguments{
\item{keepCells}{The names of the columns to keep (these should represent celltypes shared across all the datasets you plan to compare)}

\item{species}{The species of the cell type data}

\item{useCell}{The name of the column/celltype that will be compared across datasets}

\item{annotLevel}{The annotation level of the celltype data file to use}

\item{sharedName}{The file name to rename the cell which will be used}

\item{datasetGroup}{The name of the group to which this celltype data file belongs. Assumption is that multiple datasets (i.e. nuclei vs cellbody) will be compared. Should be either 'Cell' or 'Nuclei'}

\item{ctdFilePath}{Path to the cell type data file (formatted as expected for EWCE CTD files)}

\item{name1}{The 'title' of the dataset to be used in figures}
}
\value{
specificity The specificity table of the celltype data, restricted to just the cells in keepCells, and recalculated
}
\description{
\code{compare.datasets} Takes two cell type data files and tests whether the genes contained in 'geneListHGNC' are more specific to one dataset
than the other. This is restricted to the cell types named in 'col1' and 'col2'
}
\examples{
ctdFilePath = "/Users/ns9/Box Sync/Patrick & Julien/celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0.rda"
keepCells = c("astrocytes_ependymal","interneurons","microglia","Oligodendrocytes","pyramidal SS","Oligodendrocyte Precursor")
mouseSS_dist = prepare.specificity.comparison.across.datasets(ctdFilePath,keepCells,species="mouse",datasetName="Cortex KI Mouse",useCell="pyramidal SS",sharedName="Pyramidal Neuron")
}
