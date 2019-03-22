#' compare.datasets
#'
#' \code{compare.datasets} Takes two cell type data files and tests whether the genes contained in 'geneListHGNC' are more specific to one dataset
#' than the other. This is restricted to the cell types named in sbaredName
#'
#' @param dataset1 Cell type data, loaded by prepare.specificity.comparison.across.datasets()
#' @param sharedName A celltype name that is in both ctd specificity files
#' @param dataset1$datasetName The 'title' of the dataset to be used in figures
#' @param dataset2 Cell type data file (formatted as expected for EWCE CTD files)
#' @param col2 The name of the celltype to use from dataset2
#' @param dataset2$datasetName The 'title' of the dataset to be used in figures
#' @param level2annot Array of cell types, with each sequentially corresponding a column in the expression matrix
#' @param reps Number of replacates to sample
#' @param geneListHGNC Gene list (HGNC symbols) to test the specificity of between the two datasets
#' @param orthologs Data table with MGI.symbol and HGNC.symbol columns. Should only contain 1:1 homologs
#' @return p Probility (based on bootstrapping) that there is a difference in the specificity of the gene list between the two datasets
#' @return z Z-score for the difference in the specificity of the gene list between the two datasets
#' @examples
#' a=1+1
#' @export
#' @import dplyr
#' @import stats
#compare.datasets <- function(dataset1,species1,col1,dataset1$datasetName,dataset2,species2,col2,dataset2$datasetName,reps=1000,geneListHGNC=NULL,geneListMGI=NULL,orthologs){
compare_datasets <- function(dataset1,dataset2,reps=1000,geneListHGNC=NULL,geneListMGI=NULL,orthologs,sharedName){
    # Check the orthologs table
    check_1to1_orthologs(orthologs)
    
    # Check that the sharedName column is in both datasets
    if(!sharedName %in% colnames(dataset1$specificity)){stop("ERROR: sharedName is not found in columns of dataset1$specificity")}
    if(!sharedName %in% colnames(dataset2$specificity)){stop("ERROR: sharedName is not found in columns of dataset2$specificity")}
    
    # Check what species the gene list is
    if(is.null(geneListHGNC) & is.null(geneListMGI)){
        stop("ERROR: gene list must be passed to either geneListHGNC or geneListMGI")
    }
    if(length(geneListHGNC)){
        list_species = "HUMAN"
    }
    if(length(geneListMGI)){
        list_species = "MOUSE"
    }
    #if(sum(geneListHGNC %in% orthologs$HGNC.symbol) > sum(geneListHGNC %in% orthologs$MGI.symbol)){
    #    list_species="HUMAN"}else{list_species="MOUSE"
    #}
    
    all_orth = orthologs
    if(dataset1$species=="mouse"){all_orth=dplyr::filter(all_orth,MGI.symbol %in% rownames(dataset1$specificity))}
    if(dataset2$species=="mouse"){all_orth=dplyr::filter(all_orth,MGI.symbol %in% rownames(dataset2$specificity))}
    if(dataset1$species=="human"){all_orth=dplyr::filter(all_orth,HGNC.symbol %in% rownames(dataset1$specificity))}
    if(dataset2$species=="human"){all_orth=dplyr::filter(all_orth,HGNC.symbol %in% rownames(dataset2$specificity))}
    if(list_species=="HUMAN"){
        dendritic_orth = all_orth[all_orth$HGNC.symbol %in% geneListHGNC,] 
        not_dendritic_orth = all_orth[!all_orth$HGNC.symbol %in% geneListHGNC,] 
    }else{
        dendritic_orth = all_orth[all_orth$MGI.symbol %in% geneListMGI,] 
        not_dendritic_orth = all_orth[!all_orth$MGI.symbol %in% geneListMGI,] 
    }
    if(dataset1$species=="mouse"){
        data1=dataset1$specificity[dendritic_orth$MGI.symbol,sharedName]
        alldata1=dataset1$specificity[not_dendritic_orth$MGI.symbol,sharedName]
    }
    if(dataset2$species=="mouse"){
        data2=dataset2$specificity[dendritic_orth$MGI.symbol,sharedName]
        alldata2=dataset2$specificity[not_dendritic_orth$MGI.symbol,sharedName]
    }
    if(dataset1$species=="human"){
        data1=dataset1$specificity[dendritic_orth$HGNC.symbol,sharedName]
        alldata1=dataset1$specificity[not_dendritic_orth$HGNC.symbol,sharedName]
    }
    if(dataset2$species=="human"){
        data2=dataset2$specificity[dendritic_orth$HGNC.symbol,sharedName]
        alldata2=dataset2$specificity[not_dendritic_orth$HGNC.symbol,sharedName]
    }
    diff = mean(data1) - mean(data2)
    diffs=rep(0,reps)
    for(i in 1:reps){
        randomgenes = sample(1:dim(all_orth)[1],dim(dendritic_orth))
        if(dataset1$species=="mouse"){randdata1=dataset1$specificity[all_orth[randomgenes,]$MGI.symbol,sharedName]}
        if(dataset2$species=="mouse"){randdata2=dataset2$specificity[all_orth[randomgenes,]$MGI.symbol,sharedName]}
        if(dataset1$species=="human"){randdata1=dataset1$specificity[all_orth[randomgenes,]$HGNC.symbol,sharedName]}
        if(dataset2$species=="human"){randdata2=dataset2$specificity[all_orth[randomgenes,]$HGNC.symbol,sharedName]}        
        diffs[i] = mean(randdata1) - mean(randdata2)
    }
    
    # Plot the histogram data
    hist_data = data.frame(specificity=c(data1,data2),Dataset=c(rep(dataset1$datasetName,length(data1)),rep(dataset2$datasetName,length(data2))))
    graph_theme = theme_bw(base_size = 12, base_family = "Helvetica") +
        theme(legend.position = c(0.7, 0.85), panel.grid.major = element_line(size = .5, color = "grey"),
              axis.line = element_line(size=.7, color = "black"), text = element_text(size=14),
              axis.title.x = element_text(vjust = -0.35), axis.title.y = element_text(vjust = 0.6))+
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"))
    hist_plot = ggplot(hist_data)+geom_density(aes(x=specificity,fill=Dataset),alpha=0.3)+graph_theme+
        ylab("Frequency Density")
    
    # Plot the histogram data
    not_hist_data = data.frame(specificity=c(alldata1,alldata2),Dataset=c(rep(dataset1$datasetName,length(alldata1)),rep(dataset2$datasetName,length(alldata2))))
    comb_hist_data = rbind(cbind(hist_data,type="Dendritic"),cbind(not_hist_data,type="Not dendritic"))
    comb_hist_data$joint = sprintf("%s - %s",comb_hist_data$Dataset,comb_hist_data$type)
    alpha_vals = c(0.2,0.6)
    names(alpha_vals)=c("Dendritic","Not dendritic")
    comb_hist_plot = ggplot(comb_hist_data)+geom_density(aes(x=specificity,fill=Dataset,alpha=type))+graph_theme+
        ylab("Frequency Density")+scale_alpha_manual(values=alpha_vals)
    
    
    fill_vals = c("plum1","steelblue1","red","navy")
    names(fill_vals)=unique(comb_hist_data$joint)
    comb_hist_plot = ggplot(comb_hist_data)+geom_density(aes(x=specificity,fill=joint),alpha=0.5)+graph_theme+
        ylab("Frequency Density")+scale_fill_manual(values=fill_vals)
    
    # Plot the bootstrap analysis
    boot_data = data.frame(specificity=c(diffs),type=c(rep("Random",reps)))
    boot_plot = boot_plot = ggplot(boot_data)+geom_density(aes(x=specificity,fill=type))+geom_vline(xintercept=diff,col="red") + 
        guides(fill=FALSE) + graph_theme + ylab("Frequency Density") + xlab("Reduction in Specificity")
    
    # Plot together
    #plot_grid(hist_plot, comb_hist_plot, boot_plot, labels = c("A", "B", "C"), align = "h")
    
    # Get statistics
    boot_mean = mean(boot_data$specificity)
    boot_sd = sd(boot_data$specificity)
    boot_z  = (diff-boot_mean)/boot_sd
    p=sum(diff<diffs)/length(diffs)
    
    # Comparison name
    if(nchar(as.character(dataset1$datasetName))>6 | nchar(as.character(dataset2$datasetName))>6){
        comp_name = sprintf("%s vs\n %s",dataset1$datasetName,dataset2$datasetName)
    }else{comp_name = sprintf("%s vs %s",dataset1$datasetName,dataset2$datasetName)}
    
    # Prepare results
    results = list()
    results$values = data.frame(p=p,z=boot_z,boot_sd=boot_sd,boot_mean=boot_mean,name1=dataset1$datasetName,name2=dataset2$datasetName,comparison_name=comp_name)
    results$hist_plot  = hist_plot
    results$comb_hist_plot  = comb_hist_plot
    results$boot_plot = boot_plot
    
    return(results)
}