
# coding: utf-8

# In[27]:

#Jose Aguilar, jra3528
#Read in DNA sequence file
dnafile="dna.mod.fasta" #change to input file dna.mod.fasta
dnadata=open(dnafile,'r') #open file in read mode
dnadata.readline() #ignore fasta header
dnaseq=dnadata.readline() #collect the sequence in a variable
#dnaseq


# In[28]:

#List of lists with start (0 indexed) and end coordinates for coding sequences
#First entry indicates that the first coding sequence starts at 162 and ends at base 892
example_cds=[[6,11],[15,23]]
cds=[[162, 892],
 [38756, 38792],
 [45459, 45518],
 [46699, 46770],
 [47246, 47396],
 [74360, 74484],
 [78703, 78794],
 [79600, 79797],
 [81249, 81427],
 [83313, 83522],
 [86137, 86319],
 [89094, 89291],
 [89678, 89785],
 [90057, 90127],
 [110896, 111200],
 [112852, 112915],
 [118812, 118964],
 [113811, 114405],
 [128934, 129118],
 [129436, 129568],
 [134961, 135014],
 [142317, 142462],
 [143420, 143647],
 [145831, 145999],
 [146836, 146864]]

cds_new=[[162, 892],
 [43757, 43793],
 [45459, 45518],
 [46699, 46770],
 [47246, 47396],
 [74360, 74484],
 [78703, 78794],
 [79600, 79797],
 [81249, 81427],
 [83313, 83522],
 [86137, 86319],
 [89094, 89291],
 [89678, 89785],
 [90057, 90127],
 [110896, 111200],
 [112852, 112915],
 [113811, 113963],
 [114345, 114405],
 [128934, 129118],
 [129436, 129568],
 [134961, 135014],
 [142317, 142462],
 [143420, 143647],
 [145831, 145999],
 [146836, 146864]]


# In[29]:

#Dictionary containing the codon table for translation
#Key is the codon and value is the corresponding amino acid
codontable={'AAA': 'K',
 'AAC': 'N',
 'AAG': 'K',
 'AAU': 'N',
 'ACA': 'T',
 'ACC': 'T',
 'ACG': 'T',
 'ACU': 'T',
 'AGA': 'R',
 'AGC': 'S',
 'AGG': 'R',
 'AGU': 'S',
 'AUA': 'I',
 'AUC': 'I',
 'AUG': 'M',
 'AUU': 'I',
 'CAA': 'Q',
 'CAC': 'H',
 'CAG': 'Q',
 'CAU': 'H',
 'CCA': 'P',
 'CCC': 'P',
 'CCG': 'P',
 'CCU': 'P',
 'CGA': 'R',
 'CGC': 'R',
 'CGG': 'R',
 'CGU': 'R',
 'CUA': 'L',
 'CUC': 'L',
 'CUG': 'L',
 'CUU': 'L',
 'GAA': 'E',
 'GAC': 'D',
 'GAG': 'E',
 'GAU': 'D',
 'GCA': 'A',
 'GCC': 'A',
 'GCG': 'A',
 'GCU': 'A',
 'GGA': 'G',
 'GGC': 'G',
 'GGG': 'G',
 'GGU': 'G',
 'GUA': 'V',
 'GUC': 'V',
 'GUG': 'V',
 'GUU': 'V',
 'UAA': 'STOP',
 'UAC': 'Y',
 'UAG': 'STOP',
 'UAU': 'Y',
 'UCA': 'S',
 'UCC': 'S',
 'UCG': 'S',
 'UCU': 'S',
 'UGA': 'STOP',
 'UGC': 'C',
 'UGG': 'W',
 'UGU': 'C',
 'UUA': 'L',
 'UUC': 'F',
 'UUG': 'L',
 'UUU': 'F'}


# In[30]:

# Function will count the deseried characters in the string seq using the string count method
def countBases(seq):
    a = seq.count('A')
    t = seq.count('T')
    g = seq.count('G')
    c = seq.count('C')
    print("A:", a, "T:", t, "G:", g, "C:", c) #prints the repective counts of each base
    GCper = ((g+c)/(a+t+g+c)*100) #calculates the GC percrent and prints it out
    print("%GC:",GCper)
        


# In[31]:

# function will convert any 'T' char. into 'U' and keep the rest of the chracters the same 
def mRNAconvert(seq):
    output = ""           #creates blank string
    for i in seq:         #for loop goes through each char of the string
        if(i== 'T'):      #if char is 'T' then adds 'U' to output string otherwise it adds the correspondent string
            output += 'U'
        else:
            output += i
    return output         #returns the string output 


# In[32]:

# function converts the mRNA string to the corresponding protien seq based on a codon dictionary
def convertProtein(seq):
    output = "";
    i = 0                  #starts index of string at 0
    while i < len(seq):    #while loop that stops if i is greater than the length of seq
        output += codontable[seq[i:i+3]] #adds the corresponding codon given the three base sequence to the output string
        i += 3             #increments i to to skip already processed codons
    return output


# In[33]:

# uses the 2D array to slpice the unneeded seq of the mRNA 
def splicer(RNAseq): 
    output = ""
    for n in range(len(cds_new)): #for loop goes through the list of items 
          start= cds_new[n][0]  #start index of coding region is the first index of the array element 
          end= cds_new[n][1]+1  #end index of the coding region is the last index of element
          output += RNAseq[start:end] #add it to the string
    return output


# In[34]:

#calls methods needed to run expierment, using the dnaseq for the file
countBases(dnaseq)
preRNA = mRNAconvert(dnaseq)            
print("i.  pre-mRNA sequence:", preRNA)
RNA = splicer(preRNA)
print("ii. spliced mRNA sequence:",RNA)
print("")
print("Protein Seq:", convertProtein(RNA))


# In[35]:

# The protein sequence is coded by the ATXN2 Gene (Ataxin 2), this gene is part of a group of genes related to 
# neurlogical and neuromuscular diseases causes by a repetitive mutation sequence of bases. Due to this mutation
# of the GAC repeat of 23 and more effecting the behavior of motor neurons. This neurodegenerative dosorder effects
# motor neurons in the brain and spinal cord leading the ALS (Amyotrophic lateral sclerosis) or ALS 13. 


# In[ ]:



