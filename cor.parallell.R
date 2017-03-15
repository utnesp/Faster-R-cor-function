cor.parallell <- function(df, gene_one_id, gene_list_ids = "", file = "test.txt", correlation_type = "pearson", read.file = F, annotate = F, no_cores = "", use = "na.or.complete") {
    start = Sys.time()
    ## This function iterates one row at a time and then writes out each correlation to a file. Its input requires a data frame (can not be a matrix) 
    ## It is done in a manner that requires low memory due to the fact that only one row is is read and calculated at a time, and each single correlation is then appended 
    ## to a file (instead of heaping up in memory)
    ## All your cores, except for one will be used -- unless specified manually.
    
    library(foreach)
    library(doParallel)
    
    print("Correlation methods may be one of pearson, kendall, spearman")
    print("Input (df) must be a data frame or matrix with row.names = gene names, colnames = samples")
    print("No cores used may be specified manually or it will be automatically designated using all available cores - 1")
    
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
    system(paste("mkdir -p", dirname(file)), ignore.stdout = T, ignore.stderr = T)
    system(paste("trash", file), ignore.stdout = T, ignore.stderr = T)
    file.create(file)
    
    nr_times = floor(nrow(df) / no_cores) - 1
    remainder = nrow(df) %% no_cores
    
    ## Clean-up any existing temp files
     for (i in 1:(no_cores+1)) {
        tryCatch(system(paste("trash", paste(file, i, "temp.txt", sep = "_")), ignore.stdout = T, ignore.stderr = T))
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
    system(paste("trash", file.cat), ignore.stdout = T, ignore.stderr = T)
    
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
