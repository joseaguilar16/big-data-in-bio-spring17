
# coding: utf-8

# In[134]:

# Jose Aguilar
# jra3528
#Formated Data matrix by hand, no parsing used


# In[135]:

import math
import pandas
import numpy
import scipy.stats
import matplotlib  #for part B plot
import seaborn # for partB plot as well


# In[136]:

## -----------------------------------------------------------------
## load expression data
## -----------------------------------------------------------------
logRPMCountsFile = "GSE60424_GEOSubmit_FC1to11_normalized_counts.txt"
logRPM= pandas.read_csv(
    logRPMCountsFile,
    sep = "\t",
    header = 0,       ## first row has all the header information
    index_col = 0,    ## first column has gene names as row index info
)


# In[137]:

logRPM # print log RPM
logRPM.head()


# In[138]:

## -----------------------------------------------------------------
## load sample annotations
## -----------------------------------------------------------------
annotFile = "GSE60424_series_matrix.txt"
annot = pandas.read_csv(annotFile, sep="\t", header=0, index_col=0)


# In[139]:

annot


# In[140]:

## pull tissue type info from annot DataFrame
cellType = annot.ix['Celltype']
## simplify tissue names
cellType = cellType.str.replace("^celltype: ", "")


# In[141]:

## pull genotype info from annot DataFrame
##genotype=condition
status = annot.ix['Status']
status = status.str.replace("^diseasestatus: ", "")


# In[142]:

status


# In[143]:


cellSamples = cellType[cellType == 'Whole Blood'] ## only get c9 & healthy samples from the genotype meta data
healthySamples = status[status == 'Healthy Control']
msSample = status[status == 'MS pretreatment']
ALSSample = status[status == 'ALS']
healthySet = set(healthySamples.index) ## created set for each
cellSet = set(cellSamples.index)
MSSet = set(msSample.index)
ALSSet = set(ALSSample.index)
cellSet


# In[144]:

msSamp = MSSet.intersection(cellSet) ## got the intersections of c9 and cereb
ALSSamp = ALSSet.intersection(cellSet) ## got intersection of healthy and cereb 
HSamp = healthySet.intersection(cellSet)


# In[145]:

sALS = list(ALSSamp) # create them to lists
sH = list(HSamp)
sMS = list(msSamp)
print(sMS)


# In[146]:

samples = sH + sMS + sALS# add list together
print(samples)


# In[147]:

sortLogRPM = logRPM[samples] ## stores samples and perpective gene data
sortLogRPM


# In[148]:

def numZeros(arry):
    count = 0
    for i in range(len(arry)):
        if arry[i] == 0:
            count = count + 1
    return count


# In[163]:

sortLogRPM = sortLogRPM[sortLogRPM.apply(numZeros,1) < 2]
sortLogRPM


# In[164]:

sortLogRPM.to_csv('TAMLSFilteredCount.csv', sep = '\t')


# In[165]:

nova = []
for gene in sortLogRPM.index:
    arry = numpy.array(sortLogRPM.loc[gene])
    nova.append(scipy.stats.f_oneway(arry[0:4],arry[4:7],arry[7:]))
nova


# In[173]:

def getNova(arr):
    return scipy.stats.f_oneway(arr[0:4],arr[4:7],arr[7:])[1]


# In[174]:

getNova(numpy.array(sortLogRPM.loc['ENSG00000000419']))


# In[193]:

pValues = sortLogRPM.apply(getNova, 1)
pValues = pValues[pValues <=0.05]
pValsBool = sortLogRPM.apply(getNova, 1) <=0.05
log = sortLogRPM[pValsBool]
pValues


# In[187]:

#Cedit given to group : Nick Dawes, Megan Chan, Shrey Desai, Elias Sanchez, Chirs Apgar, Jose Augilar
FC = [] #list of fold change values for each gene
for i in range(0, len(log)):
    currentC9 = log.ix[i, sALS] #c9 cerebellum expression values from logRPM for each gene
    currentH = log.ix[i, sH] #healthy cerebellum expression values from logRPM for each gene
    meanC9 = (numpy.mean(currentC9)) 
    meanHealthy = (numpy.mean(currentH))
    FC.append(meanC9 / meanHealthy) #fold change calculation for each genes c9 and healthy expression mean
logFC = [] #list for log fold change
#Checks for unworkable values and fills the log fold change list with a 0 or the correct log fold change value
for i in range(0, len(log)):
    if FC[i] == 0 or numpy.isnan(FC[i]) or numpy.isinf(FC[i]):
        logFC.append(0) 
    else:
        logFC.append(math.log(FC[i],2))


# In[195]:

# creates dataFrame with listed value
resultTable = pandas.DataFrame({'p-value': pValues, 'log-FC': logFC}, index = log.index)


# In[203]:

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


# In[206]:

# Credit given to group: Nick Dawes, Megan Chan, Shrey Desai, Elias Sanchez, Chirs Apgar, Jose Augilar
sig_gene_list = []##init list holding genes with criteria: p-val <= 0.05 & |log2(fc) >=2|
for i in range(0, len(resultTable.index)):##iterate through all genes
    if resultTable.ix[i, 'p-value'] <= .05 and math.fabs(resultTable.ix[i, 'log-FC'] >= 1): ##test conditions
        sig_gene_list.append(resultTable.index[i])##append if pass condition
sig_gene_list
seaborn.clustermap(log.ix[sig_gene_list]).savefig("clustermapTAMALS.png")##go back to main df(logRPM filtered by relevant samples) and pick out "significant genes")


# In[199]:

# Writen by group : Nick Dawes, Megan Chan, Shrey Desai, Elias Sanchez, Chirs Apgar, Jose Augilar
#The cluster map found by our group for c9ALS cerebellum samples as compared to healthy cerebellum samples seems 
#to mostly agree with the equivalent cluster map as shown in Figure 1g in the Petrucelli ALS paper. A small set 
#of genes that are significantly differentially expressed in the c9ALS samples show clustering,whereas the healthy 
#samples do not show anywhere close to the same levels of differential expression in the clustered genes. There is a 
#difference between the group cluster map and the Petrucelli cluster map, however, in the rotation of the map, where
#it appears that the genes were represented on the y-axis, and the samples on the x-axis in our clustermap. 
#The Petrucelli cluster map appears to have the genes on the x-axis while the various samples are represented on the 
#y-axis.


# In[ ]:



