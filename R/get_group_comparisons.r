get_group_comparisons <- function(datasets){
    # Find which group each dataset was assigned to
    groupTypes = rep("",length(datasets))
    datasetNames = rep("",length(datasets))
    for(i in 1:length(datasets)){
        groupTypes[i]    = datasets[[i]]$datasetGroup
        datasetNames[i]  = datasets[[i]]$datasetName
    }
    uniqueGroupTypes = unique(groupTypes)
    
    # Confirm there are only two groups
    if(length(uniqueGroupTypes)!=2){
        stop("ERROR: There should only be two datasetGroup's")
    }
    
    # Find all possible combinations
    allCombs = combn(1:length(datasets),m=2)
    # Map onto them the groupType
    combGroups = matrix(groupTypes[allCombs],nrow=2,ncol=dim(allCombs)[2])
    
    # Ensure all mixed group comparisons are in the same direction
    for(j in 1:dim(combGroups)[2]){
        if(length(unique(combGroups[,j]))==2){
            if(combGroups[1,j]==uniqueGroupTypes[2]){
                tmp1 = allCombs[2,j]
                tmp2 = allCombs[1,j]
                allCombs[1,j] = tmp1
                allCombs[2,j] = tmp2
            }
        }
    }
    
    # Re-map onto them the groupType
    combGroups = matrix(groupTypes[allCombs],nrow=2,ncol=dim(allCombs)[2])
    combNames  = matrix(datasetNames[allCombs],nrow=2,ncol=dim(allCombs)[2])
    getLabel <- function(x){
        if(length(x)!=2){stop("ERROR: Function only works with inputs of length 2")}
        y = sprintf("%s vs %s",x[1],x[2])
        return(y)
    }
    
    comparisons = list()
    comparisons$indexed = allCombs
    comparisons$groups  = combGroups
    comparisons$group_labels = apply(combGroups,2,getLabel)
    comparisons$labels = apply(combNames,2,getLabel)
    return(comparisons)
}