#' load_genelist
#'
#' \code{load_genelist} Loads a gene list from a file
#' 
#' Expects the first line to be descriptive
#'
#' @param path Folder where it is stored
#' @param lookuptable Dataframe containing HGNC.symbol and MGI.symbol as columns
#' @param speciesWanted Either 'human' or 'mouse'
#' @return res Results
#' @examples
#' # list_path = sprintf("%s/%s.txt",path,listN)
#' # geneListHGNC = load.genelist(list_path,orthologs,speciesWanted="human")
#' @export
load_genelist <- function(path,lookuptable,speciesWanted="human"){
    # Check speciesWanted is either human or mouse
    if(sum(speciesWanted %in% c("human","mouse"))!=1){
        stop("ERROR: speciesWanted is not equal to either 'human' or 'mouse'")
    }
    
    # Check lookup table has HGNC.symbol and MGI.symbol
    if(sum(c("MGI.symbol","HGNC.symbol") %in% colnames(lookuptable))!=2){
        stop("ERROR: lookup table must have MGI.symbol and HGNC.symbol as column headers")
    }
    
    geneList = read.csv(list_path,stringsAsFactors = F)[,1][-1]
    
    # Find the species
    species="human"
    if(sum(geneList %in% orthologs$MGI.symbol) > sum(geneList %in% orthologs$HGNC.symbol)){species = "mouse"}
    
    if(speciesWanted=="human" & species=="human"){return(geneList)}
    if(speciesWanted=="mouse" & species=="mouse"){return(geneList)}
    if(speciesWanted=="human" & species=="mouse"){
        geneListHGNC = orthologs[orthologs$MGI.symbol %in% geneList,]$HGNC.symbol
        return(geneListHGNC)
    }
    if(speciesWanted=="mouse" & species=="human"){
        geneListMGI = orthologs[orthologs$HGNC.symbol %in% geneList,]$MGI.symbol
        return(geneListMGI)
    }
    stop("ERROR: I don't understand how the function got this far without returning")
}