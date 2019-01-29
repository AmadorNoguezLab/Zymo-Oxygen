########## RNAseq statistical analysis #############

setwd("C:/Users/Amador-Noguez Lab/Desktop/Oxygen Paper/Submission/RNAseq_analysis")

data = read.delim("HTSeq-raw-counts.txt")
data = read.delim("RSEM-TPM.txt")
data = read.csv("HTSeq-counts-noStranded.csv")

IDs = as.vector(data[,1])
samples = colnames(data)[-1]

transposed = as.data.frame(t(data[,-1]))
colnames(transposed) = IDs
num = as.data.frame(apply(transposed, c(1,2), as.numeric))
log2 = log(transposed,2)

#create factors based on sample names to specify condition, time and replicate.

n.sample = length(samples)

####condition########
condition = rep("", n.sample)
An = which(grepl("An", samples))
O2 = which(grepl("O2", samples))

condition[An] = "An"
condition[O2] = "O2"

condition = factor(condition, labels = c("Control", "Treatment"))
####################


#######time###########
time = rep(0, n.sample)
zero = which(grepl("_0_", samples))
five = which(grepl("_5_", samples))
fifteen = which(grepl("_15_", samples))
thirty = which(grepl("_30_", samples))
forty.five = which(grepl("_45_",samples))
sixty = which(grepl("_60_", samples))
one.twenty = which(grepl("_120_", samples))

time[zero] = 0
time[five] = 5
time[fifteen] = 15
time[thirty] = 30
time[forty.five] = 45
time[sixty] = 60
time[one.twenty] = 120

time = factor(time, labels = c("0","5","15","30","45","60","120"))
#######################


########replicate######
replicate = rep("", n.sample)
A = which(grepl("_A", samples))
B = which(grepl("_B", samples))
C = which(grepl("_C", samples))
D = which(grepl("_D", samples))
E = which(grepl("_E", samples))

replicate[A] = "A"
replicate[B] = "B"
replicate[C] = "C"
replicate[D] = "D"
replicate[E] = "E"

replicate = factor(replicate, labels = c("A","B","C","D","E"))
########################

#add factors to the data frame
processed = cbind(condition, time, replicate, log2)

data = processed[processed$condition=="Treatment",] 
data = data[order(data$time),]

require(lme4); require(lmerTest); require(qqman)

#template for the formula to be used to generate a generalized linera model
model.formula = paste(IDs," ~ time + (1|replicate)",sep="")

# Number of protein.ids
n = length(IDs)

# Generate a vector to save p-values
pvals = rep(NA,n)

for (i in 1:n){
  now.model <- lmer(as.formula(model.formula[i]), data=data)
  now.summary <- anova(now.model)
  
  if (!is.null(now.summary$`Pr(>F)`)) {pvals[i] <- now.summary$`Pr(>F)`}
  else {pvals [i] = 0.99; print(paste(IDs[i],"  null!!!"))}
  print(i)
#  print(pvals[i])
}

names(pvals) = IDs
test["ZMO0428"]
# Bejamini-Hotchberg Adjustment
adj.pvals <- p.adjust(pvals,method="BH")


pvals = pvals[order(names(pvals))]
adj.pvals = adj.pvals[order(names(adj.pvals))]

fold.change = read.delim("RSEM-L2FC-all.txt")
#fold.change = fold.change[order(fold.change[,1]),]

#fold.change[,1] == names(pvals)[F]

results = cbind(fold.change,pvals,adj.pvals)
results["ZMO0428",]

write.table(results, "RNA-TPM-L2FC-FDR.txt", sep = "\t", col.names = T, row.names = F)
