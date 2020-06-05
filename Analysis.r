#!/usr/bin/env Rscript
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)


check = FALSE
argument_s = commandArgs(trailingOnly=TRUE)
for (ar in argument_s){
    if (check==TRUE){
        data = read.table(ar, header=TRUE)
         = data[,c(1,7)]
        check = FALSE
    }
    else{
        read = read.table(ar, header=TRUE)
        read = read[,c(1,7)]
        data = inner_join(data, read, by='Geneid')
    }
}
colnames(data) = gsub("\\.bam$", "", colnames(data))
colnames(data) = gsub("_NameSorted", "", colnames(data))
colnames(data) = gsub("out\\.BAMs\\.MAPPED_", "", colnames(data))
rownames(data) = data$Geneid
data = data[,2:ncol(data)]
key = c("HBR", "UHRR")
types = c("Collibri", "KAPA")

for (t in types){
    m = select(data,contains(t))
    v = c()
    for (name in colnames(data)){
        if (grepl(key[1], name, fixed=T)){
            v = append(v, key[1])
        }
        else if(grepl(key[2], name, fixed=T)){
            v = append(v, key[2])
        }
        else{
            print("key NOT FOUND")
        }
    }
    K = factor(v)
    coldata = data.frame(row.names=colnames(m), K)
    dds = DESeqDataSetFromMatrix(data=m, colData=coldata, design=~K)
    dds = DESeq(dds)
    res = results(dds)
    res = res[order(res$padj), ]

    plot = EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'pvalue', xlim = c(-5, 8))
    ggsave(paste(t,"_vulcano.png", sep=""), plot=plot, device="png")

}
