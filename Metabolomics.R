##########################################################################################################
# Metabolomics data analysis script.
# Input: one or more tab-delimited text files containing LC-MS data exported from Maven.
# Samples are columns, first column is metabolite name / labeled form.
# ONE SINGLE BLANK LINE between different metabolites.
# specific to the oxygen exposure timecourse experimental design.
#
# Output: tab-delimited text file with average log2 fold change relative to time zero for each timepoint

#clear variable list
#rm(list = ls())

#set working directory
#getwd()
setwd("C:/Users/Amador-Noguez Lab/Desktop/Oxygen Paper/Submission/metabolomics")

#######################################################################################
## import data from tab-delimited text file as a matrix
## One row represents one labeled form of one metabolite
## Each column is one repicalte of a single sample condition/timepoint (e.g. An 0 A)
## First column is metabolite or labeled form name
#######################################################################################
## output data normalized by both internal labeled standard and OD,
## and calculated as a fold-change relative to time zero
## adding condition, time, and replicate information for each sample
## transposed such that columns are Log2-fold-change metabolite abundances 
## and each row is one sample
import.metabolomics.data <- function(dataFile, ODnormalization){
  
  #dataFile = "O2_timecourse/intracellular_metabolites.txt" 
  #ODnormalization= "O2_timecourse/OD-intracellular.txt"
  
  #dataFile = "O2_timecourse/RNAseq_metabolites.txt"
  #ODnormalization= "O2_timecourse/OD-RNA.txt"
  
  #dataFile = "O2_timecourse/IPP-HMBPP_metabolites.txt"
  #ODnormalization = "O2_timecourse/OD-IPP.txt"
  dataTable = as.matrix(read.delim(dataFile, header = FALSE))
  
  ##Add a final blank line to the matrix (for finding the last fully-labeled form)
  blankRow = rep("", ncol(dataTable))
  dataTable = rbind(dataTable,blankRow)
  
  #store sample names - exclude first column with metabolite names
  row1 = dataTable[1,]
  sampleNames = row1[2:length(row1)]
  
  #store the number of samples in the data set
  n.sample = length(sampleNames)
  
  #store vector with all metabolite names (and labeled forms)
  col1 = dataTable[,1]
  
  #Extract only metabolite names (and exclude blank lines and labeled forms)
  blanks = which(col1 == "")
  parent = which(grepl("PARENT", col1))
  labeled = which(grepl("label", col1))
  exclude = c(blanks, parent, labeled)
  metaboliteNames = col1[-exclude]
  
  #store the number of metabolites in the data set
  n.met = length(metaboliteNames)

  #make one matrix for the unlabeled form (from the samples)
  #need to convert from character to numeric and then reassemble the matrix.
  unlabeled = mapply(dataTable[parent,2:ncol(dataTable)], FUN = as.numeric)
  unlabeled = matrix(data = unlabeled, nrow = n.met, ncol = n.sample)
  
  #add names
  rownames(unlabeled) = metaboliteNames
  colnames(unlabeled) = sampleNames
  
  #make a separate matrix, with the same dimensions, for the fully-labeled forms (internal standard)
  #again, need to convert from character to numeric and then reassemble the matrix.
  full = blanks-1
  fullyLabeled = mapply(dataTable[full,2:ncol(dataTable)], FUN = as.numeric)
  fullyLabeled = matrix(data = fullyLabeled, nrow = n.met, ncol = n.sample)
  
  #add names
  rownames(fullyLabeled) = metaboliteNames
  colnames(fullyLabeled) = sampleNames
  
  ####################################################################################################
  #remove metabolites with low signals
  ####################################################################################################
  threshold = 4000
  remove = NULL
  
  for(r in 1:n.met){
    u.bad = length(which(unlabeled[r,] < 100)) >= 1 || length(which(unlabeled[r,] < threshold)) > 3
    f.bad = length(which(fullyLabeled[r,] < 100)) >= 1 || length(which(fullyLabeled[r,] < threshold)) > 3

    if (u.bad || f.bad){
       remove = c(remove,r)
    }
  }
  
  if(length(remove) > 0){
    cat("Removed", length(metaboliteNames[remove])  ,"metabolites:", metaboliteNames[remove], "- Signal below threshold of", threshold)
  
    #unlabeled[remove,] = NA
    #fullyLabeled[remove,] = NA
    unlabeled = unlabeled[-remove,]
    fullyLabeled = fullyLabeled[-remove,] 
    metaboliteNames = metaboliteNames[-remove]
    n.met = length(metaboliteNames)
  }
  ####################################################################################################
  
  
  #obtain normalized singnals for each metabolite by dividing unlabeled signal by internal standard
  normalized = unlabeled/fullyLabeled
  
  ######################################################################
  #normalize by OD
  ######################################################################
  OD = as.matrix(read.delim(ODnormalization, header = T, row.names = 1))
  
  colnames(OD)
  colnames(normalized)
  
  normalized = rbind(OD, normalized)
  
  ODnormalized = apply(normalized, 1, function(x) x/normalized[1,])
  ODnormalized = t(ODnormalized[,2:ncol(ODnormalized)])

  ######################################################################
  
  ####################################################################################################
  #create vectors based on sample name to specify time, condition, and replicate
  ####################################################################################################
  condition = rep("", n.sample)
  An = which(grepl("An", sampleNames))
  O2 = which(grepl("O2", sampleNames))
  
  condition[An] = "An"
  condition[O2] = "O2"
  
  conds = factor(condition)
  
  time = rep(0, n.sample)
  zero = which(grepl("_0", sampleNames))
  one = which(grepl("_1", sampleNames))
  five = which(grepl("_5", sampleNames))
  ten = which(grepl("_10", sampleNames))
  fifteen = which(grepl("_15", sampleNames))
  thirty = which(grepl("_30", sampleNames))
  fourty.five = which(grepl("_45", sampleNames))
  sixty = which(grepl("_60", sampleNames))
  one.twenty = which(grepl("_120", sampleNames))
  
  time[zero] = 0
  time[one] = 1
  time[five] = 5
  time[ten] = 10
  time[fifteen] = 15
  time[thirty] = 30
  time[fourty.five] = 45
  time[sixty] = 60
  time[one.twenty] = 120
  
  replicate = rep("", n.sample)
  A = which(grepl("_A", sampleNames))
  B = which(grepl("_B", sampleNames))
  C = which(grepl("_C", sampleNames))
  D = which(grepl("_D", sampleNames))
  
  replicate[A] = "A"
  replicate[B] = "B"
  replicate[C] = "C"
  replicate[D] = "D"
  
  reps = factor(replicate)

  ####################################################################################################
  
  #combine sample info
  info = rbind(condition, time, replicate)
  colnames(info) = sampleNames
  
  factors = rbind(condition = factor(condition, labels = c("An","O2")), time, replicate = reps)
  colnames(factors) = sampleNames
  
  #For each replicate (and each condition),
  #calcuate fold change relative to time zero
  #add it to the foldChange data frame
  foldChange = data.frame(0)

    
  for(con in levels(conds)){
    for(rep in levels(reps)){
      
      subset = ODnormalized[ , info["replicate",] == rep & info["condition", ] == con ]

      if(length(subset) > 0){
        fc = apply(subset, 2, function(x) x/subset[,1])
        foldChange = cbind(foldChange, fc)
      }
      
    }
  }
  
  #remove inital zero column
  foldChange = foldChange[,2:ncol(foldChange)]
  
  #take log base 2 of fold change
  log2FC = log(foldChange, base = 2)
  
  #add sample information
  data = rbind.data.frame(log2FC, factors)
  
  #transpose - rows are samples and columns are metabolites.
  data = t(data)
  
  #label An and O2 conditions
  data = as.data.frame(data)
  data$condition = factor(data$condition, labels = c("An","O2"))
  
  return(data)

}

#import data
RNA = import.metabolomics.data(dataFile = "RNAseq_metabolites.txt", ODnormalization= "OD-RNA.txt")
intra = import.metabolomics.data(dataFile = "intracellular_metabolites.txt", ODnormalization= "OD-intracellular.txt")
prot = import.metabolomics.data(dataFile = "proteomics_metabolites.txt", ODnormalization= "OD-proteomics.txt")
IPP = import.metabolomics.data(dataFile = "IPP_full_metabolites.txt", ODnormalization= "OD-IPP.txt")

#obtain sorted list of unique timepoints across all experiments
time = c(RNA$time, intra$time, prot$time, IPP$time)
time = unique(time)
time = sort(time)

#obtain a list of all metabolites that appear at least once
metabolites = c(colnames(RNA)[1:(ncol(RNA)-3)],colnames(intra)[1:(ncol(intra)-3)], colnames(prot)[1:(ncol(prot)-3)], colnames(IPP)[1:(ncol(IPP)-3)])
metabolites = unique(metabolites)

#all experiments have only An and O2 conditions.
conditions = levels(RNA$condition)

#construct a vector "samples" with all possible time/condition combinations
samples = NULL
n.sample = 0
valid.time = NULL
valid.cond = NULL

for (c in conditions){
  for(t in time){
  
    rna = c %in% RNA$condition[RNA$time == t]
    int = c %in% intra$condition[intra$time == t]
    pro = c %in% prot$condition[prot$time == t]
    ipp = c %in% IPP$condiion[IPP$time == t]
    
    if(rna || int || pro || ipp){
      n.sample = n.sample + 1
      samples[n.sample] = paste(c, t)
      valid.time[n.sample] = t
      valid.cond[n.sample] = c
    }
    
  }
}


#initialize 2 data frames to fill in with means and standard error
means = data.frame(matrix(NA, nrow = n.sample, ncol = length(metabolites)))
colnames(means) = metabolites
rownames(means) = samples

stderr = means
n.s = means

sets = rep( list(list()), length(metabolites) )
names(sets) = metabolites
thirty = ""



#calculate the average Log2 fold change per metabolite
for(m in metabolites){
  for(i in 1:n.sample){
    t = valid.time[i]
    c = valid.cond[i]
    
    
    rna = RNA[RNA$condition == c & RNA$time ==t, m]
    int = intra[intra$condition == c & intra$time == t, m]
    pro = prot[prot$condition == c & prot$time == t, m]
    ipp = IPP[IPP$condition == c & IPP$time == t, m]
    
    rep1 = RNA[RNA$condition == c & RNA$time ==t & !is.null(RNA[[m]]), "replicate"]
    rep2 = 4+intra[intra$condition == c & intra$time == t & !is.null(intra[[m]]), "replicate"]
    rep3 = 7+prot[prot$condition == c & prot$time == t & !is.null(prot[[m]]), "replicate"]
    rep4 = 10+IPP[IPP$condition == c & IPP$time == t & !is.null(IPP[[m]]), "replicate"]
    
    all = c(rna, int, pro, ipp)
    rep = c(rep1,rep2,rep3,rep4)
    
    n = length(all)
    
    if(n > 0){
      means[i,m] =  mean(all)
      stderr[i,m] = sd(all)/sqrt(n)
      n.s[i,m] = n
      #cat(all, "n = ", n, "time = ", t, "cond = ", c, "Met = ", m, "\n")
      
      sets[[m]] = list(time = c(sets[[m]]$time, rep(t, n)), condition = c(sets[[m]]$condition, rep(c, n)), replicate = c(sets[[m]]$replicate, rep) , signal = c(sets[[m]]$signal, all), thirty = c(sets[[m]]$thirty, rep(thirty, n)))
      #str(sets)
    }
    
    
  }
}

#str(sets)

means = t(means)
stderr = t(stderr)
n.s = t(n.s)


output = matrix(c(0,1,5,10,15,30,45,60,120,rep(NA,36)), nrow = 9, ncol = 5, byrow = F)
An = c(1,6,8,9)

for(m in metabolites){

  output[,2] = means[rownames(means)==m][5:13]
  output[,3] = stderr[rownames(stderr)==m][5:13]
  output[An,4] = means[rownames(means)==m][1:4]
  output[An,5] = stderr[rownames(stderr)==m][1:4]
  write.table(output, paste(m,"_prism.txt",sep = ""),sep = "\t")
}

write.table(means, "metabolite-means.txt", sep = "\t")
write.table(stderr, "metabolite-sderr.txt", sep = "\t")
write.table(n.s, "metabolite-n.txt", sep = "\t")
####################  STATISTICAL ANALYSIS  ######################
#str(sets)

require(lme4); require(lmerTest); require(qqman)


n = length(metabolites)

pvals <- rep(NA,n)



for(i in 1:n){

  now.metabolite = as.data.frame(sets[[metabolites[i]]])
  treatment = now.metabolite[now.metabolite$condition == "O2",]
  treatment$time = as.factor(treatment$time)
  model = lmer(signal ~ time + (1|replicate), data = treatment)
  anova = anova(model)
  pvals[i]  = anova$`Pr(>F)`
  if(pvals[i] == 0) {pvals[i] = 2.2e-16}

}



FDR = p.adjust(pvals, method = "BH")
print(pvals)
results = cbind(metabolites, pvals, FDR)
results = results[order(FDR),]

write.csv(as.data.frame(results), "significant_metabolites.csv")
write.csv(means, "significant_heat_map.csv")



