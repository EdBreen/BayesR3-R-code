makeFixedDesign = function(fixedEffects, toFile=NULL)
{	
    # Create a fixed effects file suitable for bayesR3
    # this function assumes that 'fixedEffects' is a data frame containing 
    # 1 or more factor columns of fixed effects. 
    # Each factor column should contain 2 or more levels.
    # It is also assumed that column 1 of 'fixedEffects' contains the record IDs.
    # The first column of the design matrix is for the mean
    # Returns a data frame containing the design matrix
    
    fix = data.frame(fixedEffects[, -1])  ## strip off the ID column
    id = which(sapply(fix, is.factor)) # find factor columns
    df2 = NULL
    if (length(id)) {
        for (i in id) {
            t = table(fix[, i])
            fix[, i] = factor(fix[, i], levels = names(t)[order(t, decreasing = T)])
        }
        dm = model.matrix( ~ ., fix)
        df2 = cbind(ID = fixedEffects[, 1], as.data.frame(dm))
        if (!is.null(toFile))
            write.table(df2, file = toFile,row.names = F,col.names = F,sep = ' ',quote = F)
    }
    invisible(df2)
}
