
library(purrr)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(broom)

target_1 <- "adhd2019"
target_2 <- "adhd2016"

brain <- c(
    "Brain_Amygdala",
    "Brain_Anterior_cingulate_cortex_BA24",
    "Brain_Caudate_basal_ganglia",
    "Brain_Cerebellar_Hemisphere",
    "Brain_Cerebellum",
    "Brain_Cortex",
    "Brain_Frontal_Cortex_BA9",
    "Brain_Hippocampus",
    "Brain_Hypothalamus",
    "Brain_Nucleus_accumbens_basal_ganglia",
    "Brain_Putamen_basal_ganglia",
    "Brain_Spinal_cord_cervical_c-1",
    "Brain_Substantia_nigra"
)

smultxican_function  <- function(x, y) {
    x <- fread(x)
    scz.gene.david <- fread(file = "/home/nick/magma-library-uniq-ensembl.csv")
    x <- x[complete.cases(x$pvalue), ]
    x$gene <- str_extract(x$gene, pattern = "[A-Z0-9]+")
    colnames(x) <- gsub("gene", "ENSEMBL_GENE_ID", colnames(x))
    x <- left_join(x, scz.gene.david, by = "ENSEMBL_GENE_ID")
    x <- cbind(x[, 1:18], apply(x[, 19:39], 2, function(x) ifelse(is.na(x), 0, x)))

    y <- get(load(y))
    qvalue <- p.adjust(x$pval, method = "fdr", n = length(x$pval))
    target_variable <- transform(x, qvalue = qvalue, FDRtheoretical = y$FDR)
    return(target_variable)
}

gene_function <- function(x, y) {
    x <- fread(x)
    scz.gene.david <- fread("/home/nick/magma-library-uniq-entrez.csv")
    colnames(x) <- gsub("GENE", "ID", colnames(x))
    x <- left_join(x, scz.gene.david, by = "ID")
    x <- cbind(x[, 1:10], apply(x[, 11:31], 2, function(x) ifelse(is.na(x), 0, x)))
    
    y <- get(load(y))
    qvalue <- p.adjust(x$P, method = "fdr", n = length(x$P))
    target_variable <- transform(x, qvalue = qvalue, FDRtheoretical = y$FDR)
    return(target_variable)
}

predixcan_function <- function(x, y) {
    x <- fread(x)
    x <- x[complete.cases(x$pvalue), ]
    scz.gene.david <- fread(file = "/home/nick/magma-library-uniq-ensembl.csv")
    x$gene <- str_extract(x$gene, pattern = "[A-Z0-9]+")
    colnames(x) <- gsub("gene", "ENSEMBL_GENE_ID", colnames(x[1]))
    x <- left_join(x, scz.gene.david, by = "ENSEMBL_GENE_ID")
    x <- cbind(x[, 1:12], apply(x[, 13:33], 2, function(x) ifelse(is.na(x), 0, x)))

    y <- get(load(y))
    qvalue <- p.adjust(x$pvalue, method = "fdr", n = length(x$pvalue))
    target_variable <- transform(x, qvalue = qvalue, FDRtheoretical = y$FDR)

    return(target_variable)
}

addition <- as_mapper(~ {
    param_forjoin <- ..2
    param_foroutput <- ..3
    a <- lapply(1:9, function(y) {
        switch(y, subset(..1[[1]], qvalue < 0.05),
            subset(..1[[2]], FDRtheoretical < 0.05),
            subset(..1[[2]], qvalue < 0.05),
            subset(..1[[3]], qvalue < 0.05),
            subset(..1[[4]], FDRtheoretical < 0.05),
            subset(..1[[4]], qvalue < 0.05),
            subset(..1[[5]], qvalue < 0.05),
            subset(..1[[6]], FDRtheoretical < 0.05),
            subset(..1[[6]], qvalue < 0.05)
        )
    })
    c <- map2_dbl(
        list(a[[1]], a[[1]], a[[4]], a[[4]], a[[7]], a[[7]]),
        list(a[[3]], a[[2]], a[[6]], a[[5]], a[[9]], a[[8]]),
        ~ inner_join(.x, .y, by = param_forjoin) %>% nrow()
    )
    if (nrow(a[[1]]) == 0) {
        result1 <- matrix(c(c[1], c[2], (c[2] / c[1]), nrow(a[[1]]), NA), ncol = 5, byrow = F)
        result2 <- matrix(c(c[3], c[4], (c[4] / c[3]), nrow(a[[1]]), NA), ncol = 5, byrow = F)
        result3 <- matrix(c(c[5], c[6], (c[6] / c[5]), nrow(a[[1]]), NA), ncol = 5, byrow = F)
    } else {
        result1 <- matrix(c(c[1], c[2], (c[2] / c[1]), nrow(a[[1]]), (tidy(prop.test(c(c[2], c[1]), c(nrow(a[[1]]), nrow(a[[1]])), alternative = "greater"))$p.value)), ncol = 5, byrow = F)
        result2 <- matrix(c(c[3], c[4], (c[4] / c[3]), nrow(a[[1]]), (tidy(prop.test(c(c[4], c[3]), c(nrow(a[[4]]), nrow(a[[4]])), alternative = "greater"))$p.value)), ncol = 5, byrow = F)
        result3 <- matrix(c(c[5], c[6], (c[6] / c[5]), nrow(a[[1]]), (tidy(prop.test(c(c[6], c[5]), c(nrow(a[[7]]), nrow(a[[7]])), alternative = "greater"))$p.value)), ncol = 5, byrow = F)
    }
    return(list(result1, result2, result3))
})

for (j in c("gene", "smultixcan","predixcan")) {
    if (j == "gene") {
        x <- vapply(c(target_1, target_2), function(x) {
            paste0("/home/nick/yuping/", x, "/gene/", x, ".genes.out")
        }, character(1), USE.NAMES = FALSE)

        x <- flatten_chr(replicate(3, x, simplify = FALSE))
        y <- flatten_chr(sapply(c("Nolasso-theore.Rdata", "BioNolasso-theore.Rdata", "BioLasso-theore.Rdata"), function(x) {
            two_name <- c(paste0("/home/nick/yuping/", target_1, "/gene/FDRreg/", x), paste0("/home/nick/yuping/", target_2, "/gene/FDRreg/", x))
            return(list(two_name))
        }))

        b <- map2(x, y, gene_function)
        result <- reduce(addition(b, "ID", "gene"), rbind)
        write.csv(result, file = paste0("/home/nick/yuping/validation/adhd_gene_valiation.csv"), quote = F)
        rm(result)
    } else if (j == "smultixcan") {
        x <- vapply(c(target_1, target_2), function(x) {
            paste0("/home/nick/yuping/", x, "/smultixcan/", x, ".allbrain.txt")
        }, character(1), USE.NAMES = FALSE)

        x <- flatten_chr(replicate(3, x, simplify = FALSE))
        y <- flatten_chr(sapply(c("Nolasso-theore.Rdata", "BioNolasso-theore.Rdata", "BioLasso-theore.Rdata"), function(x) {
            two_name <- c(paste0("/home/nick/yuping/", target_1, "/smultixcan/FDRreg/", x), paste0("/home/nick/yuping/", target_2, "/smultixcan/FDRreg/", x))
            return(list(two_name))
        }))

        b <- map2(x, y, smultxican_function)
        result <- reduce(addition(b, "ENSEMBL_GENE_ID", "smultxican"), rbind)
        write.csv(result, file = paste0("/home/nick/yuping/validation/adhd_smultixcan_valiation.csv"), quote = F)
        rm(result)
    } else if(j == "predixcan") {
      for (i in brain) {
          x <- vapply(c(target_1, target_2), function(x) {
              paste0("/home/nick/yuping/", x, "/predixcan/gtex_v7_", str_to_lower(x), "_for_", i, ".csv")
          }, character(1), USE.NAMES = FALSE)
      
          x <- flatten_chr(replicate(3, x, simplify = FALSE))
          y <- flatten_chr(sapply(c("Nolasso-theore.Rdata", "BioNolasso-theore.Rdata", "BioLasso-theore.Rdata"), function(x) {
              two_name <- c(paste0("/home/nick/yuping/", target_1, "/predixcan/FDRreg/", i, x), paste0("/home/nick/yuping/", target_2, "/predixcan/FDRreg/", i, x))
              return(list(two_name))
          }))
          
          b <- map2(x, y, predixcan_function)
          result <- reduce(addition(b, "ENSEMBL_GENE_ID", "predixcan"), rbind)
          if (exists("final")) {
              final <- rbind(final, result)
          } else {
              final <- result
          }
        }
        write.csv(final, file = "/home/nick/yuping/validation/adhd_predixcanvalidation.csv", quote = F)
    }
}