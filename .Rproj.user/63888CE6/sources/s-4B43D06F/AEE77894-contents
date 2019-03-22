#' test_all_comparisons
#'
#' \code{test_all_comparisons} Runs the main dendritic depletion analysis
#'
#' @param allDataSets Created using prepare_specificity_comparison_across_datasets
#' @param orthologs Dataframe containing HGNC.symbol and MGI.symbol as columns
#' @param sharedName The cell type being compared across datasets
#' @return res Results
#' @examples
#' allDataSets = list(mouseSS_dist,tasic_mouse_dist,dronc_mouse_dist,dronc_human_dist,divseq_dist)
#' res = test_all_comparisons(allDataSets,orthologs,sharedName="Pyramidal Neuron")
#' @export
test_all_comparisons <- function(allDataSets,orthologs,sharedName="Pyramidal Neuron"){
    # Load the gene list
    #list_path = sprintf("%s/%s.txt",path,listN)
    #geneListHGNC = load.genelist(list_path,orthologs,speciesWanted="human")
    listName="dendrite_enriched_transcripts_HGNC_1to1only"
    data("dendriticGenesHGNC")
    
    # Find all possible comparisons
    comparisons = get_group_comparisons(allDataSets)
    
    # For each comparison
    comparisons$z = rep(0,length(comparisons$labels))
    comparisons$p = rep(0,length(comparisons$labels))
    for(cc in 1:length(comparisons$labels)){
        dd1 = allDataSets[[comparisons$indexed[1,cc]]]
        dd2 = allDataSets[[comparisons$indexed[2,cc]]]
        a = compare_datasets(dataset1=dd1,dataset2=dd2, geneListHGNC=dendriticGenesHGNC,orthologs=orthologs,sharedName=sharedName)
        comparisons$z[cc] = a$values$z
        comparisons$p[cc] = a$values$p
    }
    
    res = data.frame(labels=comparisons$labels,groupLabels=comparisons$group_labels,z=comparisons$z,p=comparisons$p,list=listName)
    return(res)
}