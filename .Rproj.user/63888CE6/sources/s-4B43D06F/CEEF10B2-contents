#' compare.datasets
#'
#' \code{compare.datasets} Takes two cell type data files and tests whether the genes contained in 'geneListHGNC' are more specific to one dataset
#' than the other. This is restricted to the cell types named in 'col1' and 'col2'
#'
#' @param ctdFilePath Path to the cell type data file (formatted as expected for EWCE CTD files)
#' @param keepCells The names of the columns to keep (these should represent celltypes shared across all the datasets you plan to compare)
#' @param species The species of the cell type data
#' @param name1 The 'title' of the dataset to be used in figures
#' @param useCell The name of the column/celltype that will be compared across datasets
#' @param annotLevel The annotation level of the celltype data file to use
#' @param sharedName The file name to rename the cell which will be used
#' @param datasetGroup The name of the group to which this celltype data file belongs. Assumption is that multiple datasets (i.e. nuclei vs cellbody) will be compared. Should be either 'Cell' or 'Nuclei'
#' @return specificity The specificity table of the celltype data, restricted to just the cells in keepCells, and recalculated
#' @examples
#' ctdFilePath = "/Users/ns9/Box Sync/Patrick & Julien/celltype_data_allKImouse_wtHypo_MergedStriatal_1to1only_level1_thresh0_trim0.rda"
#' keepCells = c("astrocytes_ependymal","interneurons","microglia","Oligodendrocytes","pyramidal SS","Oligodendrocyte Precursor")
#' mouseSS_dist = prepare.specificity.comparison.across.datasets(ctdFilePath,keepCells,species="mouse",datasetName="Cortex KI Mouse",useCell="pyramidal SS",sharedName="Pyramidal Neuron")
#' @export
#' @import dplyr
#' @import stats
prepare_specificity_comparison_across_datasets <- function(ctd,keepCells,species,datasetName=NA,useCell,annotLevel=1,sharedName="sharedCell",datasetGroup){
    
    if(is.na(datasetName)){stop("datasetName must be provided")}
    if(is.na(datasetGroup)){stop("datasetGroup must be provided")}
    if(!datasetGroup %in% c("Cell","Nuclei")){stop("datasetGroup must be either 'Cell' or 'Nuclei'")}
  
    specificity = ctd[[annotLevel]]$mean_exp
    
    # Check that all keepCells are in specificity
    if(sum(!keepCells %in% colnames(specificity))>0){
        print("keepCells: ")
        print(keepCells)
        print("colnames(specificity)")
        print(colnames(specificity))
        stop("ERROR: not all keepCells are in ctd[[annotLevel]]$mean_exp")
    }
    
    mouseSS_dist = specificity[,keepCells] 
    mouseSS_dist = mouseSS_dist/apply(mouseSS_dist,1,sum)
    mouseSS_dist = mouseSS_dist[!is.na(mouseSS_dist[,useCell]),]
    colnames(mouseSS_dist)[colnames(mouseSS_dist)==useCell] = sharedName
    
    # Prepare outputs
    output = list()
    output$specificity = mouseSS_dist
    output$annotLevel = annotLevel
    #output$ctdFilePath = ctdFilePath
    output$species = species
    output$datasetName = datasetName
    output$useCell = useCell
    output$keepCells = keepCells
    output$datasetGroup = datasetGroup
    return(output)
}