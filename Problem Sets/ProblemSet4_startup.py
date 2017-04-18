
# coding: utf-8

# In[337]:

import math
import pandas
import numpy
import scipy.stats
import matplotlib  #for part B plot
import seaborn # for partB plot as well


# In[338]:

## -----------------------------------------------------------------
## load expression data
## -----------------------------------------------------------------
logRPMCountsFile = "GSE67196_Petrucelli2015_ALS_genes.logRpm.tsv"
logRPM= pandas.read_csv(
    logRPMCountsFile,
    sep = "\t",
    header = 0,       ## first row has all the header information
    index_col = 0,    ## first column has gene names as row index info
)


# In[339]:

logRPM # print log RPM
logRPM.head


# In[340]:

## -----------------------------------------------------------------
## load sample annotations
## -----------------------------------------------------------------
annotFile = "gse67196_annot.tsv"
annot = pandas.read_csv(annotFile, sep="\t", header=0, index_col=0)


# In[341]:

## pull tissue type info from annot DataFrame
tissue = annot["!Sample_characteristics_ch1"]
## simplify tissue names
tissue = tissue.str.replace("^tissue: ", "")


# In[342]:

## pull genotype info from annot DataFrame
##genotype=condition
genotype = annot["!Sample_characteristics_ch1.1"]
genotype = genotype.str.replace("^genotype: ", "")


# In[343]:

genotype.head()
c9ALSSampales = genotype[genotype == 'c9ALS'] ## only get c9 & healthy samples from the genotype meta data
healthySamples = genotype[genotype == 'Healthy']
healthySet = set(healthySamples.index) ## created set for each
setC9 = set(c9ALSSampales.index)
setC9


# In[344]:

#Select from the data frame, the samples pertaining to the conditions 
#and brain regions assigned to you 
c9want = tissue[tissue == 'Cerebellum'] ## chose cerbeulm from meta data
setCereb = set(c9want.index) # created set 
setCereb


# In[345]:

samplesC9 = setCereb.intersection(setC9) ## got the intersections of c9 and cereb
samplesHealth = setCereb.intersection(healthySet) ## got intersection of healthy and cereb 


# In[346]:

print(samplesC9)
print(samplesHealth)


# In[347]:

sC9 = list(samplesC9) # create them to lists
sH = list(samplesHealth)
sC9


# In[348]:

samples = sC9 + sH # add list together
print(samples)


# In[349]:

sortLogRPM = logRPM[samples] ## stores samples and perpective gene data
sortLogRPM


# In[350]:

coeff = numpy.array([-1,-1,-1,-1,-1,-1,-1,-1,1,1,1,1,1,1,1,1])
coeff ## creates data to compute corraltion coeff


# In[351]:

numpy.corrcoef(coeff,test)[0,1] ## simple test


# In[352]:

cCoeff = [] # creates list to corraltation coeff for each gene in the list =
for i in sortLogRPM.index: # for loop goes through each index of the log data
    cCoeff.append(numpy.corrcoef(numpy.array(sortLogRPM.loc[i]),coeff)[0,1]) # adds corrCoef to list


# In[353]:

def getT(pStat): # define function to determine the t-value
    return (math.sqrt(14)/math.sqrt(1 - math.pow(pStat,2)))*pStat


# In[354]:

def getP(tStat): # define function to detemine the p-value from the t-value
    return 2*(1-scipy.stats.t.cdf(numpy.abs(tStat), 14))


# In[355]:

tVals = [] # creates array to store values for tVals and pVals
pVals = []
for p in cCoeff: # goes through each value of cCoeff of the genes
    tVals.append(getT(p)) ## calulates and adds it to the gene
    
for t in tVals: # goes through each value of tVals of the genes
    pVals.append(getP(t)) ## calulates and adds it to the gene


# In[356]:

#Cedit given to group : Nick Dawes, Megan Chan, Shrey Desai, Elias Sanchez, Chirs Apgar, Jose Augilar
FC = [] #list of fold change values for each gene
for i in range(0, len(sortLogRPM)):
    currentC9 = sortLogRPM.ix[i, sC9] #c9 cerebellum expression values from logRPM for each gene
    currentH = sortLogRPM.ix[i, sH] #healthy cerebellum expression values from logRPM for each gene
    meanC9 = (numpy.mean(currentC9)) 
    meanHealthy = (numpy.mean(currentH))
    FC.append(meanC9 / meanHealthy) #fold change calculation for each genes c9 and healthy expression mean
logFC = [] #list for log fold change
#Checks for unworkable values and fills the log fold change list with a 0 or the correct log fold change value
for i in range(0, len(sortLogRPM)):
    if FC[i] == 0 or numpy.isnan(FC[i]) or numpy.isinf(FC[i]):
        logFC.append(0) 
    else:
        logFC.append(math.log(FC[i],2))


# In[357]:

# creates dataFrame with listed value
resultTable = pandas.DataFrame({'Correlation Coefficient': cCoeff, 't-statistic': tVals, 'p-value': pVals, 'log-FC': logFC}, index =sortLogRPM.index)


# In[358]:

resultTable


# In[359]:

#Genes that holds and does not hold the null hypothesis
#Null hypthesis: there is no correlation between c9 and gene expression w/ p >= 0.05

#Holds: PHLPP2 (CC)0.389297 (p-val)0.136117 (t-value)1.581364
#since the corrlation coefficient is very low towards 0, indicating very little correlation, thus creating a t-stat that
# which states that there is little evidence to how relation, reulting in a p-value greater then 0.05 confidence value
# ultimate stating that the gene is not abnormal in indicating gene expression

#Does not hold: (CC)0.600022 (p-value)0.014003 (t-value)2.806403
#since the corellation coefficient is very high and close to one, there is alot of correlation, thus getting a t-stat of
# a greater value since there is a lot of evidence to show the correlation, finally the p-value is way less than 0.05
# indicating that the genes abnormal expression is a valuble to be diffrently expressed.


# In[360]:

# Credit given to group: Nick Dawes, Megan Chan, Shrey Desai, Elias Sanchez, Chirs Apgar, Jose Augilar
sig_gene_list = []##init list holding genes with criteria: p-val <= 0.05 & |log2(fc) >=2|
for i in range(0, len(resultTable.index)):##iterate through all genes
    if resultTable.ix[i, 'p-value'] <= .05 and math.fabs(resultTable.ix[i, 'log-FC'] >= 2): ##test conditions
        sig_gene_list.append(resultTable.index[i])##append if pass condition
seaborn.clustermap(sortLogRPM.ix[sig_gene_list]).savefig("clustermap.png")##go back to main df(logRPM filtered by relevant samples) and pick out "significant genes")


# In[ ]:



