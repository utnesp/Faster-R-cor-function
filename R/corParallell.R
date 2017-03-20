## Created 2017-03-17
## Authour: Peter Utnes
## All rights reserved

#' @title Faster R cor() function using parallellization
#' @author Peter Utnes \email{utnesp@gmail.com}
#' @description #' This function iterates one row at a time and then writes out each correlation to a file. It is done in a manner that requires low memory due to the fact that only one row is is read and calculated at a time, and each single correlation is then appended to a file (instead of heaping up in memory).
#' @param df a numeric data frame or matrix with rows and columns corresponding to variables and samples, respectively.
#' @param var variable to do the correlation with. Must be part of df.
#' @param var.list which variables to correlate to. Defaults to row.names.
#' @param file txt file for storing results
#' @param correlation_type correlation methods may be one of "pearson" (default), "kendall", "spearman". 
#' @param no_cores number of cores used may be specified manually or it will be designated using all available cores - 1 (default)
#' @param read.file defaults to FALSE. If TRUE, then assigns object to global environment using the name specified in the file
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". See detailed description in ?cor
#' @export
#' @import foreach
#' @import doParallel
#' @import easybiomart
#' @examples cor.parallell(counts, "ENSG00000134323", file = "/path/to/file/MYCN.cor.txt"). \cr # Default mart \cr nif ( exists("mart") == "FALSE") { \cr    mart = useMart("ENSEMBL_MART_ENSEMBL", dataset='hsapiens_gene_ensembl') \cr }
cor.parallell <- function(df, var, var.list = NULL, file = "test.txt", correlation_type = "pearson", read.file = F, annotate = F, no_cores = "", use = "na.or.complete") {
    start = Sys.time()
    
    if (class(df) != "matrix") {
        cat(paste("\nConverting ", class(df), "using as.matrix:", sep = ""))
        cat(paste("\n> head(", class(df), ")\n", sep = ""))
        tryCatch(print(df[1:5, 1:5]))
        df <- as.matrix(df)
        cat(paste("\n> head(", class(df), ")\n", sep = ""))
        tryCatch(print(df[1:5, 1:5]))
    }
    
    if (is.null(var.list)) {
        var.list <- rownames(df)
    }
    
    var_one <- df[rownames(df) == var, ]
    
    if (no_cores == "") {
        no_cores <- detectCores() - 1
        registerDoParallel(no_cores)
    }
    
    # Creates and overwrites file if it already exists
    dir.create(dirname(file), showWarnings = F)
    # system(paste("mkdir -p", dirname(file)), ignore.stdout = T, ignore.stderr = T)
    # system(paste("trash", file), ignore.stdout = T, ignore.stderr = T) # on my own system I use this as it is more safe. 
    # The trash command is similar to rm, but instead of deleting the file permanently it moves the file to the trash can.
    # It is not in use here because not everybody is able to use this command. To implement this type of command create a
    # sh file like this:
        # #!/bin/bash
        # mv "$@" /Users/put001/Desktop/Trashcan
    # then create a symlink: ln -s trash /path/to/sh/file 
    # make sure the symlink is in your $PATH
    file.create(file, overwrite = T)
    
    nr_times = floor(nrow(df) / no_cores) - 1
    remainder = nrow(df) %% no_cores
    
    ## Clean-up any existing temp files
     for (i in 1:(no_cores+1)) {
        tryCatch(unlink(paste(file, i, "temp.txt", sep = "_"), force = T))
    }
    
    for (i in 0:nr_times) {
        foreach (z = 1:no_cores,
                 .combine = rbind)  %dopar%
            cat(
                round(cor(var_one, # correlate gene of interest to ...
                          df[rownames(df) == var.list[((i * no_cores) + z)], ], method = correlation_type, use = use), 2), #  gene x in the gene list, and round
                "\t", var.list[((i * no_cores) + z)], "\n", 
                sep = "", file = paste(file, z, "temp.txt", sep = "_"), append = T) # add tab and newline, then append to file
    }
    
    
        if (remainder != 0) {
                z <- length(grep(paste(file, ".", "temp.txt", sep = "_"), list.files(dirname(file), full.names = T))) + 1
            
                file.rem = paste(file, z, "temp.txt", sep = "_")
            
                # add in the rest
                for (i in 1:remainder) {
                    cat(
                        round(cor(var_one, # Pearson correlation with var one ...
                                  df[rownames(df) == var.list[((nr_times + 1) * no_cores) + i], ], method = correlation_type), 2), # compared against gene x in the gene list, and round
                        "\t", var.list[((nr_times + 1) * no_cores) + i], "\n",
                        sep = "", file = file.rem, append = T)
                }
        } 
    
    stopImplicitCluster()
    
    # file.cat <- gsub(paste(z, "temp.txt", sep = "_"), "*_temp.txt", file.rem)
    file.cat <- paste(file, '*', "temp.txt", sep = "_")
    
    cat("\nNumber of correlations in temp files:\n")
    tryCatch(system(paste("wc -l", file.cat)))
    
    cat(paste("\nNumber of variables provided to this function:", length(var.list)))
    
    ## combine all files 
    system(paste("cat", file.cat, ">", file))
    ## clean up temp files
    unlink(file.cat, force = T)
    
    if (read.file == T | annotate == T) {
        assign.name <- gsub(paste(dirname(file), "/", sep = ""), "", file)
        assign.name <- gsub(".txt", "", assign.name)
        assign(assign.name, read.delim(file, header = F))
        
        if (annotate == T) {
            t <- get(assign.name)
            colnames(t) <- c("correlation_type", "gene_name")
            
            t <- ensg2ext_name_biotype(t$gene_name, combine = T)
            
            t <- t[order(-abs(t$correlation_type)), ]
            
            colnames(t) <- c("ensembl_gene_id", "external_gene_name", "gene_biotype", correlation_type)
            assign(assign.name, t, envir = .GlobalEnv)
            
            rm(t)    
        } else {
            t <- get(assign.name)
            colnames(t) <- c("correlation_type", "gene_name")
            t <- t[order(-abs(t$correlation_type)), ]
            colnames(t) <- c(correlation_type, "gene_name")
            assign(assign.name, t, envir = .GlobalEnv)
            rm(t)
        }
        cat("\n", assign.name, "assigned as object to .GlobalEnv\n") 
    }
    
    cat(paste("\nNumber of correlations in final file:\n", system(paste("wc -l", file), intern = T)))
    end = Sys.time()
    cat(paste("\nFinished in", round((end - start) / 60, 1), "minutes"))
}