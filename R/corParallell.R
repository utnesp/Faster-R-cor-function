#' @title Faster R cor() function using parallellization
#' 
#' @description #' This function iterates one row at a time and then writes out each correlation to a file. It is done in a manner that requires low memory due to the fact that only one row is is read and calculated at a time, and each single correlation is then appended to a file (instead of heaping up in memory)
#' 
#' @param df Must be a data frame or matrix with row.names = gene names, colnames = samples
#' @param gene_one_id The variable to do the correlation with
#' @param gene_list_ids Which variables to correlate to. Defaults to row.names
#' @param file Where to store the results
#' @param correlation_type Correlation methods may be one of "pearson" (default), "kendall", "spearman". 
#' @param no_cores Number of cores used may be specified manually or it will be designated using all available cores - 1 (default)
#' @param read.file Defaults to FALSE. If TRUE, then assigns object to global environment using the name specified in file
#' @param use an optional character string giving a method for computing covariances in the presence of missing values. This must be (an abbreviation of) one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs". See detailed description in ?cor
#' @export
#' @import foreach
#' @import doParallel
#' @import easybiomart
#' @examples
#' cor.parallell()

cor.parallell <- function(df, gene_one_id, gene_list_ids = "", file = "test.txt", correlation_type = "pearson", read.file = F, annotate = F, no_cores = "", use = "na.or.complete") {
    start = Sys.time()
    
    if (class(df) != "matrix") {
        print(paste("Converting ", class(df), "using as.matrix", sep = ""))
        print(paste("Head ", class(df), ":", sep = ""))
        tryCatch(print(df[1:5, 1:5]))
        tryCatch(print(str(df[1:5, 1:5])))
        df <- as.matrix(df)
        print(paste("Head ", class(df), ":", sep = ""))
        tryCatch(print(df[1:5, 1:5]))
        tryCatch(print(str(df[1:5, 1:5])))
    }
    
    if (gene_list_ids == "") gene_list_ids = rownames(df)
    
    gene_one <- df[rownames(df) == gene_one_id, ]
    
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
                round(cor(gene_one, # correlate gene of interest to ...
                          df[rownames(df) == gene_list_ids[((i * no_cores) + z)], ], method = correlation_type, use = use), 2), #  gene x in the gene list, and round
                "\t", gene_list_ids[((i * no_cores) + z)], "\n", 
                sep = "", file = paste(file, z, "temp.txt", sep = "_"), append = T) # add tab and newline, then append to file
    }
    
    
        if (remainder != 0) {
                z <- length(grep(paste(file, ".", "temp.txt", sep = "_"), list.files(dirname(file), full.names = T))) + 1
            
                file.rem = paste(file, z, "temp.txt", sep = "_")
            
            
                # add in the rest
                for (i in 1:remainder) {
                    cat(
                        round(cor(gene_one, # Pearson correlation with gene one ...
                                  df[rownames(df) == gene_list_ids[((nr_times + 1) * no_cores) + i], ], method = correlation_type), 2), # compared against gene x in the gene list, and round
                        "\t", gene_list_ids[((nr_times + 1) * no_cores) + i], "\n",
                        sep = "", file = file.rem, append = T)
                }
        } 
    
    stopImplicitCluster()
    
    # file.cat <- gsub(paste(z, "temp.txt", sep = "_"), "*_temp.txt", file.rem)
    file.cat <- paste(file, '*', "temp.txt", sep = "_")
    print("Number of correlations in temp files:")
    system(paste("wc -l", file.cat))
    print(paste("Number of genes provided to this function:", length(gene_list_ids)))
    
    ## combine all files 
    system(paste("cat", file.cat, ">", file))
    ## clean up temp files
    unlink(file.cat, force = T)
    
    if (read.file == T || annotate == T) {
        assign.name <- gsub(paste(dirname(file), "/", sep = ""), "", file)
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
    }
    
    print(paste("Number of correlations in final file:", system(paste("wc -l", file), intern = T), " vs ", nrow(df), " correlations to make in your input", sep = ""))
    end = Sys.time()
    print(paste("Finished in", (end - start) / 60, "minutes"))
}
