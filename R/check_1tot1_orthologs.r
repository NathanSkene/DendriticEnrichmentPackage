#' check.1to1.orthologs
#'
#' \code{check.1to1.orthologs} Confirms that the datatable passed as input has columns for MGI.symbol, HGNC.symbol and only contains 1:1 orthologs
#'
#' @param orthologs Data table with MGI.symbol and HGNC.symbol columns. Should only contain 1:1 homologs
#' @examples
#' check.1to1.orthologs(orth)
#' @export
check_1to1_orthologs <- function(orth){
    if(sum(c("HGNC.symbol","MGI.symbol") %in% colnames(orth))!=2){
        stop("ERROR: orthologs data table does not have the correct column headers. Must contain 'HGNC.symbol' and 'MGI.symbol'")
    }
    orth2 = unique(orth[,c("HGNC.symbol","MGI.symbol")])
    if(sum(duplicated(orth2$MGI.symbol))>0){
        stop("ERROR: orthologs data table has duplicated MGI.symbols. Should only contain 1:1 orthologs")
    }
    if(sum(duplicated(orth2$HGNC.symbol))>0){
        stop("ERROR: orthologs data table has duplicated HGNC.symbols. Should only contain 1:1 orthologs")
    }    
}