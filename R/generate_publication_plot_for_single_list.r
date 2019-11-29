#' generate_publication_plot_for_single_list
#'
#' \code{generate_publication_plot_for_single_list} Generates the main plot
#'
#' @param allRes Output of test_all_comparisons
#' @param listName Used to select a particular list from listName (not relevant to current version of package)
#' @param save_folder Folder where output PDFs should be stored
#' @param angle for xaxis labels, either 90 or 45 is probably best
#' @return thePlot (as a ggplot)
#' @examples
#' allDataSets = list(mouseSS_dist,tasic_mouse_dist,dronc_mouse_dist,dronc_human_dist,divseq_dist)
#' res = test_all_comparisons(allDataSets,orthologs,sharedName="Pyramidal Neuron")
#' thePlot = generate_publication_plot_for_single_list(res)
#' @export
#' @import cowplot
#' @import ggplot2
generate_publication_plot_for_single_list <- function(allRes,listName="dendrite_enriched_transcripts_HGNC_1to1only",save_folder="",xaxis_angle=90){
    #load("~/Box Sync/DendriticEnrichment/allRes_LotsOfListsWtP.Rda")
    if(!is.na(listName)){
        subRes = allRes[allRes$list == listName,]
        print("using restrictor")
    }else{
        subRes=allRes
        print("using allRes")
    }
    subRes$q = p.adjust(subRes$p,method="BH")
    subRes$signif = "Significant"
    subRes$signif[subRes$q>0.05] = "Nonsignificant"
    group.colors <- c(Significant = "#b30000", Nonsignificant = "#000000")
    library("cowplot")
    thePlot = ggplot(subRes)+geom_bar(aes(x=labels,y=z,fill=signif),stat="identity",position="dodge", width=0.9)+
        #facet_wrap(~groupLabels, strip.position = "bottom", scales = "free_x")+
        #facet_grid(.~groupLabels, scales = "free", space = "free", switch="x")+
        theme_cowplot()+
        facet_grid(list~groupLabels, scales = "free", space = "free", switch="x")+
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        ylab("Z-Score")+xlab("")+
        scale_fill_manual(values=group.colors)+
        scale_y_continuous(breaks=seq(-20,20,1)) + guides(fill=FALSE)
    #if(length(unique(allRes$list))>1){
    #    thePlot = thePlot + facet_wrap(~list,ncol=1)
    #}
    for(widths in seq(from=3,to=12,by=3)){
        for(height in seq(from=3,to=9,by=3)){
            pdf(sprintf("%s/Fig_DendriticEnrichment_%s_%s_w%s_h%s.pdf",save_folder,listName,Sys.Date(),widths,height),width=widths,height=height)
            print(thePlot)
            dev.off()
        }
    }
    widths=10
    height=5
    pdf(sprintf("%s/Fig_DendriticEnrichment_%s_%s_w%s_h%s.pdf",save_folder,listName,Sys.Date(),widths,height),width=widths,height=height)
    print(thePlot)
    dev.off()
    return(thePlot)
}