###############################
########## FUNCTIONS ##########
###############################

## quantile normalize the gene expression values. takes a dataframe
## with genes as rows and samples as columns
quant_norm_mean = function(df){

    means = sort(apply(df,1,mean))
    df.quant_norm = data.frame(lapply(df, function(x) means[rank(x,ties.method="min")]))

    return(df.quant_norm)
}

# quantile normalize to the normal distribution by row.    
quant_norm_norm = function(df){

    ncols = ncol(df)
    norm_vec = qnorm(ppoints(ncols))

    df.rank = t(apply(df,1,rank,ties.method="min"))
    df.norm = data.frame(matrix(norm_vec[df.rank],ncol=ncol(df),nrow=nrow(df)))
    
    return(df.norm)
}

## run PEER
run_peer = function(mat,num_factors){

    set.seed(4)

    ## create model object
    model = PEER()

    ## set observed data 
    PEER_setPhenoMean(model, mat)

    ## set number of hidden factors to infer
    PEER_setNk(model,num_factors)

    ## perform inference
    PEER_update(model)

    ## get peer factors
    factors = PEER_getX(model)

    return(factors)
}


## read a tabix file and return a dataframe
tabix2df = function(file, range, colnames){

    data.frame(lines=tabix.read(file, range)) %>%
        separate(lines, into=colnames, sep="\t", convert=T)
}
