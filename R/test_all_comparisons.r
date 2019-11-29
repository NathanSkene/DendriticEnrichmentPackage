#' test_all_comparisons
#'
#' \code{test_all_comparisons} Runs the main dendritic depletion analysis
#'
#' @param allDataSets Created using prepare_specificity_comparison_across_datasets
#' @param orthologs Dataframe containing HGNC.symbol and MGI.symbol as columns
#' @param sharedName The cell type being compared across datasets
#' @param listN Name of the text file with the data to be loaded
#' @param reps Number of bootstrap replicates to use
#' @return res Results
#' @examples
#' allDataSets = list(mouseSS_dist,tasic_mouse_dist,dronc_mouse_dist,dronc_human_dist,divseq_dist)
#' res = test_all_comparisons(allDataSets,orthologs,sharedName="Pyramidal Neuron")
#' @export
test_all_comparisons <- function(allDataSets,orthologs,sharedName="Pyramidal Neuron",listN="dendrite_enriched_transcripts_HGNC_1to1only",path=NA,pSides="onesided",reps=10000){
    # Load the gene list
    if(listN=="dendrite_enriched_transcripts_HGNC_1to1only"){
        listName="dendrite_enriched_transcripts_HGNC_1to1only"
        data("dendriticGenesHGNC")
        geneListHGNC=dendriticGenesHGNC
    }else{
        if(is.na(path)){stop("path cannot be NA if genes need to be loaded from a text file")}
        list_path = sprintf("%s/%s.txt",path,listN)
        geneListHGNC = load_genelist(list_path,orthologs,speciesWanted="human")
    }
    
    # Find all possible comparisons
    comparisons = get_group_comparisons(allDataSets)
    
    # For each comparison
    comparisons$z = rep(0,length(comparisons$labels))
    comparisons$p = rep(0,length(comparisons$labels))
    for(cc in 1:length(comparisons$labels)){
        dd1 = allDataSets[[comparisons$indexed[1,cc]]]
        dd2 = allDataSets[[comparisons$indexed[2,cc]]]
        a = compare_datasets(dataset1=dd1,dataset2=dd2, geneListHGNC=geneListHGNC,orthologs=orthologs,sharedName=sharedName,pSides=pSides,reps=reps)
        comparisons$z[cc] = a$values$z
        comparisons$p[cc] = a$values$p
    }
    
    res = data.frame(labels=comparisons$labels,groupLabels=comparisons$group_labels,z=comparisons$z,p=comparisons$p,list=listN,pSides=pSides)
    return(res)
}