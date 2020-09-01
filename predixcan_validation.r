library(purrr)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(broom)
library(fmsb)


target_1 <- "adhd2019"
target_2 <- "adhd2017"
scz.gene.david <- fread(file = "/home/nick/magma-library-uniq-ensembl.csv")


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


for (i in brain) {
    a <- function(x, y) {
        x <- fread(x)
        x <- x[complete.cases(x$pvalue), ]
        x$gene <- str_extract(x$gene, pattern = "[A-Z0-9]+")
        colnames(x) <- gsub("gene", "ENSEMBL_GENE_ID", colnames(x[1]))
        x <- left_join(x, scz.gene.david, by = "ENSEMBL_GENE_ID")
        x <- cbind(x[, 1:12], apply(x[, 13:33], 2, function(x) ifelse(is.na(x), 0, x)))

        y <- get(load(y))
        qvalue <- p.adjust(x$pvalue, method = "fdr", n = length(x$pvalue))
        target_variable <- transform(x, qvalue = qvalue, FDRtheoretical = y$FDR)

        return(target_variable)
    }



    x <- unname(vapply(
        c(target_1, target_2, target_1, target_2, target_1, target_2), function(x) {
            paste0("/home/nick/yuping/", x, "/predixcan/gtex_v7_", str_to_lower(x), "_for_", i, ".csv")
        },
        character(1)
    ))


    y <- list(
        paste0("/home/nick/yuping/", target_1, "/predixcan/FDRreg/", i, "Nolasso-theore.Rdata"),
        paste0("/home/nick/yuping/", target_2, "/predixcan/FDRreg/", i, "Nolasso-theore.Rdata"),
        paste0("/home/nick/yuping/", target_1, "/predixcan/FDRreg/", i, "BioNolasso-theore.Rdata"),
        paste0("/home/nick/yuping/", target_2, "/predixcan/FDRreg/", i, "BioNolasso-theore.Rdata"),
        paste0("/home/nick/yuping/", target_1, "/predixcan/FDRreg/", i, "BioLasso-theore.Rdata"),
        paste0("/home/nick/yuping/", target_2, "/predixcan/FDRreg/", i, "BioLasso-theore.Rdata")
    )



    b <- map2(x, y, a)

    a <- lapply(1:6, function(y) {
        switch(y, subset(b[[1]], qvalue < 0.05),
            subset(b[[2]], FDRtheoretical < 0.05),
            subset(b[[2]], qvalue < 0.05),
            subset(b[[3]], qvalue < 0.05),
            subset(b[[4]], FDRtheoretical < 0.05),
            subset(b[[4]], qvalue < 0.05)
        )
    })


    c <- map2_dbl(
        list(a[[1]], a[[1]], a[[4]], a[[4]], a[[7]], a[[7]]),
        list(a[[3]], a[[2]], a[[6]], a[[5]], a[[9]], a[[8]]),
        ~ inner_join(.x, .y, by = "ENSEMBL_GENE_ID") %>% nrow()
    )


    if (nrow(a[[1]]) == 0) {
        result1 <- matrix(c(c[1], c[2], (c[2] / c[1]), nrow(a[[1]]), NA), ncol = 5, byrow = F)
        result2 <- matrix(c(c[3], c[4], (c[4] / c[3]), nrow(a[[1]]), NA), ncol = 5, byrow = F)
        result3 <- matrix()
    } else {
        result1 <- matrix(c(c[1], c[2], (c[2] / c[1]), nrow(a[[1]]), (tidy(prop.test(c(c[2], c[1]), c(nrow(a[[1]]), nrow(a[[1]]))))$p.value)), ncol = 5, byrow = F)
        result2 <- matrix(c(c[3], c[4], (c[4] / c[3]), nrow(a[[1]]), (tidy(prop.test(c(c[4], c[3]), c(nrow(a[[4]]), nrow(a[[4]]))))$p.value)), ncol = 5, byrow = F)
    }


    result <- rbind(result1, result2)
    if (exists("final")) {
        final <- rbind(final, result)
    } else {
        final <- result
    }
}

write.csv(final, file = "/home/nick/yuping/validation/asd_predixcanvalidation.csv", quote = F)