
library(purrr)
library(data.table)
library(dplyr)
library(tidyverse)
library(stringr)
library(broom)

smultxican_function  <- function(x, y) {
    x <- fread(x)
    x <- x[complete.cases(x$pvalue), ]
    x$gene <- str_extract(x$gene, pattern = "[A-Z0-9]+")
    colnames(x) <- gsub("gene", "ENSEMBL_GENE_ID", colnames(x))
    scz.gene.david <- fread(file = "/home/nick/magma-library-uniq-ensembl.csv")
    x <- left_join(x, scz.gene.david, by = "ENSEMBL_GENE_ID")
    x <- cbind(x[, 1:18], apply(x[, 19:39], 2, function(x) ifelse(is.na(x), 0, x)))

    y <- get(load(y))
    qvalue <- p.adjust(x$pval, method = "fdr", n = length(x$pval))
    target_variable <- transform(x, qvalue = qvalue, FDRtheoretical = y$FDR)
    return(target_variable)
}


a <- function(x, y) {
    x <- fread(x)
    scz.gene.david <- fread("/home/nick/magma-library-uniq-entrez.csv")
    colnames(x) </- gsub("GENE", "ID", colnames(x))
    x <- left_join(x, scz.gene.david, by = "ID")
    x <- cbind(x[, 1:10], apply(x[, 11:31], 2, function(x) ifelse(is.na(x), 0, x)))

    y <- get(load(y))
    qvalue <- p.adjust(x$P, method = "fdr", n = length(x$P))
    target_variable <- transform(x,
        qvalue = qvalue,
        FDRtheoretical = y$FDR
    )
    return(target_variable)
}






x <- list("/home/nick/yuping/asd2019/smultixcan/asd2019.allbrain.txt",
          "/home/nick/yuping/asd2015.pgc/smultixcan/asd2015.pgc.allbrain.txt",
          "/home/nick/yuping/asd2019/smultixcan/asd2019.allbrain.txt",
          "/home/nick/yuping/asd2015.pgc/smultixcan/asd2015.pgc.allbrain.txt",
          "/home/nick/yuping/asd2019/smultixcan/asd2019.allbrain.txt",
          "/home/nick/yuping/asd2015.pgc/smultixcan/asd2015.pgc.allbrain.txt")



y <- list("/home/nick/yuping/asd2019/smultixcan/FDRreg/Nolasso-theore.Rdata",
          "/home/nick/yuping/asd2015.pgc/smultixcan/FDRreg/Nolasso-theore.Rdata",
          "/home/nick/yuping/asd2019/smultixcan/FDRreg/BioNolasso-theore.Rdata",
          "/home/nick/yuping/asd2015.pgc/smultixcan/FDRreg/BioNolasso-theore.Rdata",
          "/home/nick/yuping/asd2019/smultixcan/FDRreg/BioLasso-theore.Rdata",
          "/home/nick/yuping/asd2015.pgc/smultixcan/FDRreg/BioLasso-theore.Rdata")








b <- map2(x, y, a)

a <- lapply(1:9,function(y){
       switch(y,subset(b[[1]], qvalue < 0.05),
                subset(b[[2]], FDRtheoretical < 0.05),
                subset(b[[2]], qvalue < 0.05),
                subset(b[[3]], qvalue < 0.05),
                subset(b[[4]], FDRtheoretical < 0.05),
                subset(b[[4]], qvalue < 0.05),
                subset(b[[5]], qvalue < 0.05),
                subset(b[[6]], FDRtheoretical < 0.05),
                subset(b[[6]], qvalue < 0.05) )})


c <- map2_dbl(
  list(a[[1]], a[[1]],a[[4]],a[[4]],a[[7]],a[[7]]),
  list(a[[3]], a[[2]],a[[6]],a[[5]],a[[9]],a[[8]]),
  ~ inner_join(.x, .y, by = "ENSEMBL_GENE_ID") %>% nrow()
)

result1 <- matrix(c(c[1],c[2],(c[2]/c[1]),nrow(a[[1]]),(tidy(prop.test(c(c[2],c[1]),c(nrow(a[[1]]),nrow(a[[1]]))))$p.value)),ncol=5,byrow=F)
result2 <- matrix(c(c[3],c[4],(c[4]/c[3]),nrow(a[[1]]),(tidy(prop.test(c(c[4],c[3]),c(nrow(a[[4]]),nrow(a[[4]]))))$p.value)),ncol=5,byrow=F)
result3 <- matrix(c(c[5],c[6],(c[6]/c[5]),nrow(a[[1]]),(tidy(prop.test(c(c[6],c[5]),c(nrow(a[[7]]),nrow(a[[7]]))))$p.value)),ncol=5,byrow=F)





result <- data.frame(rbind(result1,result2,result3))
write.csv(result,file="/home/nick/yuping/validation/asd_smultixcanvalidation.csv", quote = F)