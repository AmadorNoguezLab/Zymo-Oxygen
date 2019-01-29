#######################################################################################################
# Proteomics data processing script.
# Takes as input raw data downloaded from the Coon Lab data-sharing website.
# Two files: aerobic and anaerobic with Label-Free Quantitation (log-transformed signals)
# for each replicate. Protein ID is a unique number/letter code, protein name gives annotation info,
# gene name is either annotated gene name (e.g. thyA) of Z. mobilis locus tag (e.g. ZMO1353)
#
# Out-put is a table with all values - rows are samples, columns are protein signals
# factor columns are added to indicate condition, time, and replicate
#######################################################################################################

setwd("C:/Users/Amador-Noguez Lab/Desktop/Oxygen Paper/Proteomics")

#import files: rows are proteins, columns are info or samples.
aerobic = read.delim("aerobic.txt")
anaerobic = read.delim("anaerobic.txt")

#import conversion from protein ID to locus tag
locusTags = read.delim("protein_locus_tags.txt")


#combine aerobic and anaerobic files
data = merge(anaerobic,aerobic, all.x = F, all.y = F)
#data2 = merge(anaerobic,aerobic, all.x = T, all.y = T)

data = merge(locusTags, data)

#store protein and sample information
protein.IDs = data[,1]
locus.tags = data[,2]
protein.names = data[,3]
gene.names = data[,4]
sampleNames = colnames(data)[6:ncol(data)]
n.sample = length(sampleNames)

#transpose: rows now samples, columns are proteins. Use protein ID if there is no locus tag.
transposed = as.data.frame(t(data[,6:ncol(data)]))
colnames(transposed) = locus.tags
colnames(transposed)[which(locus.tags == 0)] = protein.IDs[which(locus.tags == 0)]

#create factors based on sample names to specify condition, time and replicate.

####condition########
condition = rep("", n.sample)
An = which(grepl("ANA", sampleNames))
O2 = which(grepl("OX", sampleNames))

condition[An] = "An"
condition[O2] = "O2"

condition = factor(condition, labels = c("Control", "Treatment"))
####################


#######time###########
time = rep(0, n.sample)
zero = which(grepl("_0.", sampleNames))
five = which(grepl("_5.", sampleNames))
fifteen = which(grepl("_15.", sampleNames))
thirty = which(grepl("_30.", sampleNames))
sixty = which(grepl("_60.", sampleNames))
one.twenty = which(grepl("_120.", sampleNames))

time[zero] = 0
time[five] = 5
time[fifteen] = 15
time[thirty] = 30
time[sixty] = 60
time[one.twenty] = 120

time = factor(time, labels = c("0","5","15","30","60","120"))
#######################


########replicate######
replicate = rep("", n.sample)
A = which(grepl("_A", sampleNames))
B = which(grepl("_B", sampleNames))
C = which(grepl("_C", sampleNames))

replicate[A] = "A"
replicate[B] = "B"
replicate[C] = "C"

replicate = factor(replicate, labels = c("A","B","C"))
########################

#add factors to the data frame
processed = cbind(condition, time, replicate, transposed)

#shorten names
shortNames = paste(condition, time, replicate)
rownames(processed) = shortNames

#export as tab-delimited text file
write.table(processed, "Proteomics_re-shaped.txt", sep = "\t")
