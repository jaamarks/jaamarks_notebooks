probabel_phenotype <- function(cov_file, pheno_file, study, pop, pcs, other_covs, pheno){
    #Convert your phenotype PED files (as required for rvtests) to the
    #ProbABEL format.
    
    pc_paste <- paste(pcs, collapse="+")
    othercov_paste <- paste(other_covs, collapse="_")
    out_name <- paste(study, pop, pheno, othercov_paste, pc_paste, sep="_")
    out_name <- paste0(out_name, ".txt")

    all_covs <- c(other_covs, pcs)

    # read in data files
    my_cov <- read.table(cov_file, header=T)
    my_pheno <- read.table(pheno_file, header=T)

    num.subjects <- length(my_cov[,1])

    # create new phenotype file (dataframe)
    num_cols = 2 + length(other_covs) # iid + pheno + covs
    covar.data <- data.frame(matrix(ncol = num_cols, nrow = num.subjects))
    names(covar.data) <- c("iid", pheno, other_covs)
    covar.data[1] <- my_pheno[,"iid"]
    covar.data[2] <- my_pheno[, pheno]

    for (cov in all_covs){
        covar.data[, cov] <- my_cov[,cov]
    }
    #head(covar.data)

    write.table(x = covar.data, file = out_name, quote = F, row.names = F)
    write.table(x = covar.data[,1], file = paste0("phenotype_ids_", pop), 
                quote = F, row.names = F, col.names=F)
}


# Example use below, in a Jupyter Notebook
################################################################################
################################################################################
#setwd("/shared/jmarks/hiv/uhs1234/acquisition_gwas/phenotype/final")
#
### AA
#cov_file <- "../processing/UHS1234_NGCW1_AFR_cov.ped"
#pheno_file <- "../processing/UHS1234_NGCW1_AFR_phen.ped"
#study <- "uhs1234"
#pop <- "aa"
#pcs = c("PC10", "PC9", "PC2", "PC6")
#other_covs <- c("age", "gender")
#pheno <- "hiv"
#probabel_phenotype(cov_file, pheno_file, study, pop, pcs, other_covs, pheno)
#
#
### EA
#cov_file <- "../processing/UHS1234_NGCW1_EUR_cov.ped"
#pheno_file <- "../processing/UHS1234_NGCW1_EUR_phen.ped"
#study <- "uhs1234"
#pop <- "ea"
#pcs = c("PC1", "PC9", "PC10")
#other_covs <- c("age", "gender")
#pheno <- "hiv"
#probabel_phenotype(cov_file, pheno_file, study, pop, pcs, other_covs, pheno)
