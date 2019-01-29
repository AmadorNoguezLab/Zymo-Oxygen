##edgeR

setwd("C:/Users/Amador-Noguez Lab/Desktop/Oxygen Paper/Submission/RNAseq_analysis")

library(edgeR)

data_HTSeq = read.delim("HTSeq-raw-counts.txt")
data_RSEM = read.delim("RSEM-raw-counts.txt")
data_rockhop = read.delim("O2_timecourse/rockhopper-raw-counts.txt")

row.names(data_HTSeq) = data_HTSeq[,1]
row.names (data_rockhop) = paste(row.names(data_rockhop),data_rockhop[,1],":",data_rockhop[,2])
row.names(data_RSEM) = data_RSEM[,1]

data_HTSeq = data_HTSeq[,2:length(data_HTSeq)]
data_RSEM = data_RSEM [,2:length(data_RSEM)]
data_rockhop = data_rockhop [,3:length(data_rockhop)]


calculateDEG <- function(data, name){

samples = c(1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,4,4,5,5,5,5,6,6,6,6,7,7,7,7,8,8,8,9,9,9,10,10,10)
reps = c (1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,4,5,1,2,3,4,1,2,3,4,1,2,3,4,1,2,3,1,2,3,1,2,3)

DGE = DGEList(counts = data, group = samples)

DGE = calcNormFactors(DGE)

design = model.matrix(~0+group, data = DGE$samples)

yglm = estimateDisp(DGE, design)
yglm = estimateGLMTrendedDisp(yglm, design)
yglm = estimateGLMTagwiseDisp(yglm, design)


fit = glmFit(yglm, design)
O2_5 = glmLRT(fit, contrast = c(-1,1,0,0,0,0,0,0,0,0))
O2_15 = glmLRT(fit, contrast = c(-1,0,1,0,0,0,0,0,0,0))
O2_30 = glmLRT(fit, contrast = c(-1,0,0,1,0,0,0,0,0,0))
O2_45 = glmLRT(fit, contrast = c(-1,0,0,0,1,0,0,0,0,0))
O2_60 = glmLRT(fit, contrast = c(-1,0,0,0,0,1,0,0,0,0))
O2_120 = glmLRT(fit, contrast = c(-1,0,0,0,0,0,1,0,0,0))
An_60 = glmLRT(fit, contrast = c(0,0,0,0,0,0,0,1,-1,0))
An_120 = glmLRT(fit, contrast = c(0,0,0,0,0,0,0,1,0,-1))


table_O2_5 = topTags(O2_5, n = length(row.names(O2_5$table)))
table_O2_15 = topTags(O2_15, n = length(row.names(O2_15$table)))
table_O2_30 = topTags(O2_30, n = length(row.names(O2_30$table)))
table_O2_45 = topTags(O2_45, n = length(row.names(O2_45$table)))
table_O2_60 = topTags(O2_60, n = length(row.names(O2_60$table)))
table_O2_120 = topTags(O2_120, n = length(row.names(O2_120$table)))
table_An_60 = topTags(An_60, n = length(row.names(An_60$table)))
table_An_120 = topTags(An_120, n = length(row.names(An_120$table)))

write.table(table_O2_5, paste(name,"/O2_5.txt",sep = ""), sep = "\t")
write.table(table_O2_15, paste(name,"/O2_15.txt",sep = ""), sep = "\t")
write.table(table_O2_30, paste(name,"/O2_30.txt",sep = ""), sep = "\t")
write.table(table_O2_45, paste(name,"/O2_45.txt",sep = ""), sep = "\t")
write.table(table_O2_60, paste(name,"/O2_60.txt",sep = ""), sep = "\t")
write.table(table_O2_120, paste(name,"/O2_120.txt",sep = ""), sep = "\t")
write.table(table_An_60, paste(name,"/An_60.txt",sep = ""), sep = "\t")
write.table(table_An_120, paste(name,"/An_120.txt",sep = ""), sep = "\t")


names = rownames(table_O2_5$table)
O2_5_logFC = table_O2_5$table$logFC
names(O2_5_logFC) = names
names = rownames(table_O2_15$table)
O2_15_logFC = table_O2_15$table$logFC
names(O2_15_logFC) = names
names = rownames(table_O2_30$table)
O2_30_logFC = table_O2_30$table$logFC
names(O2_30_logFC) = names
names = rownames(table_O2_45$table)
O2_45_logFC = table_O2_45$table$logFC
names(O2_45_logFC) = names
names = rownames(table_O2_60$table)
O2_60_logFC = table_O2_60$table$logFC
names(O2_60_logFC) = names
names = rownames(table_O2_120$table)
O2_120_logFC = table_O2_120$table$logFC
names(O2_120_logFC) = names
names = rownames(table_An_60$table)
An_60_logFC = table_An_60$table$logFC
names(An_60_logFC) = names
names = rownames(table_An_120$table)
An_120_logFC = table_An_120$table$logFC
names(An_120_logFC) = names

O2_5_logFC = O2_5_logFC[order(names(O2_5_logFC))]

O2_15_logFC = O2_15_logFC[order(names(O2_15_logFC))]

O2_30_logFC = O2_30_logFC[order(names(O2_30_logFC))]

O2_45_logFC = O2_45_logFC[order(names(O2_45_logFC))]

O2_60_logFC = O2_60_logFC[order(names(O2_60_logFC))]

O2_120_logFC = O2_120_logFC[order(names(O2_120_logFC))]

An_60_logFC = An_60_logFC[order(names(An_60_logFC))]

An_120_logFC = An_120_logFC[order(names(An_120_logFC))]

output= cbind(O2_5_logFC,O2_15_logFC,O2_30_logFC,O2_45_logFC,O2_60_logFC,O2_120_logFC,An_60_logFC,An_120_logFC)
write.table(output, "RSEM-L2FC-all.txt", sep = "\t", row.names = T, col.names = F)

}
calculateDEG(data_HTSeq, "HTSeq")
calculateDEG(data_RSEM, "RSEM")
calculateDEG(data_rockhop, "rockhopper")

