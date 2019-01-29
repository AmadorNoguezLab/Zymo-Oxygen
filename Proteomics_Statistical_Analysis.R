######## Proteomics statistical analysis ########

require(lme4); require(lmerTest); require(qqman)

setwd("C:/Users/Amador-Noguez Lab/Desktop/Oxygen Paper/Submission/Proteomics")

#import data
data = read.table("Proteomics_re-shaped.txt")

# Extract id names
id = colnames(data)[-c(1,2,3)]

# Extract treatment data only, order by time and convert to factor.
data = data[data$condition=="Treatment",-1] 
data = data[order(data$time),]
data$time = as.factor(data$time)

#template for the formula to be used to generate a generalized linera model
model.formula = paste(id," ~ time + (1|replicate)",sep="")

# Number of protein.ids
n = length(id)

# Generate a vector to save p-values
pvals = rep(NA,n)

for (i in 1:n){
  now.model <- lmer(as.formula(model.formula[i]), data=data)
  now.summary <- anova(now.model)
  
  if (!is.null(now.summary$`Pr(>F)`)) {pvals[i] <- now.summary$`Pr(>F)`}
  else {pvasl [i] = 0.99; print(paste(id[i],"  null!!!"))}
  print(i)
  print(pvals[i])
}

names(pvals) = id

# Bejamini-Hotchberg Adjustment
adj.pvals <- p.adjust(pvals,method="BH")


#import fold change data
fold.change = read.csv("proteomics.csv")

pvals = pvals[order(names(pvals))]
#pvals[1:20]
adj.pvals = adj.pvals[order(names(adj.pvals))]
#adj.pvals[1:20]
fold.change = fold.change[order(fold.change$name),]

proteomics = cbind(fold.change,pvals,adj.pvals)
write.table(proteomics,"protein-L2FC-FDR.txt",sep = "\t",col.names = T, row.names = F)
