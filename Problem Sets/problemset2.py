
# coding: utf-8

# In[118]:

# Jose Aguilar jra3528
import pandas as pd


# In[119]:

exampleDNA = "GTAGAGC"
exampleSeq = ["ATAGACC","ATAGAGC","GTAGACC","ATAGACC","ATAGAGC","GTAGCCC"]


# In[120]:

# fillMatrix counts the base in the postion and fills it in an array
def fillMatrix(ex):
    for sequence in ex: #gets each read from the file
        for j in range(len(sequence)): #loop through each char of the sequence
            if(sequence[j] == 'A'): #checks the character of each postion
                DNAmatrix.loc['A',j] += 1 #adds one to the base row and j column 
            if(sequence[j] == 'C'):
                DNAmatrix.loc['C',j] += 1
            if(sequence[j] == 'G'):
                DNAmatrix.loc['G',j] += 1
            if(sequence[j] == 'T'):
                DNAmatrix.loc['T',j] += 1


# In[121]:

# returns the positions of calcualted SNPs from the reads
def getSNPLoc(percent,DNA):
    snpLocs = [] #creates temp var to store lists
    for i in range(len(DNA)): # loops through each char in the DNA seq
        temp = getLetter(i,DNA) #calls temp method to get the highest percent of a base and lowest base represeted
        if(not(temp[0] == DNA[i]) and temp[1] <= percent): #adds it to the list if the lowest char does not match the DNA seq
            snpLocs.append([temp[0],i])                    # and if the percent is less then the (100%-error rate)
    return snpLocs


# In[122]:

#get the highest percent of a base and lowest base represeted in the mathMatrix
def getLetter(col,DNA):
    sumCol = DNAmatrix[col].sum() # gets the sum of the column 
    maxPer = 0                    # creates vars
    minCount = 1000
    minChar = ''
    for j in ['A','C','G','T']: #goes through each row of the mathMatrix
        temp = (DNAmatrix.loc[j,col]/sumCol)  #gets percent of each row and count for each row
        tempC = DNAmatrix.loc[j,col]
        if(temp > maxPer):                    #finds out the higest percent and lowest represented char
            maxPer = DNAmatrix.loc[j,col]/sumCol
        if(tempC < minCount):
            minCount = tempC
            minChar = j
    return [minChar,maxPer] 


# In[123]:

# prints the sequence changes and changes the sequence to match the SNPs
def printNewSeq(changes, refrence):
    newDNA = []                    #temp var to store the new changed DNA seq
    print("SNPs:")
    for i in range(len(changes)): #goes through each value of the change list of positions of SNPs
        tempSeq = ""   #prints pos and base changes
        print("SNP at position",changes[i][1],":", refrence[changes[i][1]], "is replaced by", changes[i][0])
        for j in range(len(refrence)): # for loop goes through each char anc changes the sequence based on SNPs
            if j == changes[i][1]:
                tempSeq = tempSeq + changes[i][0]
            else:
                tempSeq = tempSeq + refrence[j]
        newDNA.append(tempSeq) #adds it to the var list 
    print()
    print("Sequence variants:")  #prints the new sequences
    for k in newDNA:
        print(k)
    return newDNA


# In[124]:

#calls methods on example DNA sequence
exampleCols = range(0,7)
DNAmatrix = pd.DataFrame(index=['A','C','G','T'],columns=exampleCols) #creates dataFrame
DNAmatrix = DNAmatrix.fillna(0) #fills them to zero
fillMatrix(exampleSeq) #counts the bases in the sequence list
print("Profile Matrix:") #calls other methods
print(DNAmatrix)
temp = getSNPLoc((2/3), exampleDNA)
storeDNA = printNewSeq(temp,exampleDNA)


# In[125]:

#import your problem set I to reuse the transcribe and translate functions
#when using a python script (.py file) without an IDE like jupyter hub, 
#you can do this by using from problemset1 import function1, function2
#but this seems more complicated when using a jupyter file (.ipynb file)

#You may be able to use the bash load (%%load) command to load the script.


# In[126]:

#creates dataFram with readchromosome information 
chr12_gtf = pd.read_table('genes_ucsc.chr12.mod.gtf', names=['sequene','source','feature','start','end','score','strand','frame','attributes'])


# In[127]:

chr12_gtf.head()


# In[128]:

#creates a series based on the if the attributes column contains the ATXN2_partial gene
atxn2_cds = chr12_gtf[chr12_gtf['attributes'].str.contains('gene_id "ATXN2_partial') & (chr12_gtf['feature'] == 'CDS')]


# In[129]:

atxn2_cds


# In[130]:

#sets starts and stop index to splice
start = atxn2_cds['start'][1]
stop = atxn2_cds['end'][1]


# In[131]:

# function will convert any 'T' char. into 'U' and keep the rest of the chracters the same 
def mRNAconvert(seq):
    output = ""           #creates blank string
    for i in seq:         #for loop goes through each char of the string
        if(i== 'T'):      #if char is 'T' then adds 'U' to output string otherwise it adds the correspondent string
            output += 'U'
        else:
            output += i
    return output         #returns the string output 


# In[132]:

# function converts the mRNA string to the corresponding protien seq based on a codon dictionary
def convertProtein(seq):
    output = "";
    i = 0                  #starts index of string at 0
    while i < len(seq):    #while loop that stops if i is greater than the length of seq
        output += codontable[seq[i:i+3]] #adds the corresponding codon given the three base sequence to the output string
        i += 3             #increments i to to skip already processed codons
    return output


# In[133]:

referenceDNA='ACCCCCGAGAAAGCAACCCAGCGCGCCGCCCGCTCCTCACGTGTCCCTCCCGGCCCCGGGGCCACCTCACGTTCTGCTTCCGTCTGACCCCTCCGACTTCCGGTAAAGAGTCCCTATCCGCACCTCCGCTCCCACCCGGCGCCTCGGCGCGCCCGCCCTCCGATGCGCTCAGCGGCCGCAGCTCCTCGGAGTCCCGCGGTGGCCACCGAGTCTCGCCGCTTCGCCGCAGCCAGGTGGCCCGGGTGGCGCTCGCTCCAGCGGCCGGCGCGGCGGAGCGGGCGGGGCGGCGGTGGCGCGGCCCCGGGACCGTATCCCTCCGCCGCCCCTCCCCCGCCCGGCCCCGGCCCCCCTCCCTCCCGGCAGAGCTCGCCTCCCTCCGCCTCAGACTGTTTTGGTAGCAACGGCAACGGCGGCGGCGCGTTTCGGCCCGGCTCCCGGCGGCTCCTTGGTCTCGGCGGGCCTCCCCGCCCCTTCGTCGTCCTCCTTCTCCCCCTCGCCAGCCCGGGCGCCCCTCCGGCCGCGCCAACCCGCGCCTCCCCGCTCGGCGCCCGCGCGTCCCCGCCGCGTTCCGGCGTCTCCTTGGCGCGCCCGGCTCCCGGCTGTCCCCGCCCGGCGTGCGAGCCGGTGTATGGGCCCCTCACCATGTCGCTGAAGCCCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCAGCAGCAGCAGCAGCAGCAGCAGCCGCCGCCCGCGGCTGCCAATGTCCGCAAGCCCGGCGGCAGCGGCCTTCTAGCGTCGCCCGCCGCCGCGCCTTCGCCGTCCTCGTCCTCGGTCTCCTCGTCCTCGGCCACGGCTCCCTCCTCGGTGGTCGCGGCGACCTCCGGCGGCGGGAGGCCCGGCCTGGGCAGGTGGGTGTCGGCACCCCAGCCCCCTCCGCTCCGGGCCCGGCGTCCCCTCCCCCGCGGCCCGCGCCGCCGTCCCCGCCCCGTGACCCGCCGGGCTACCCGGGGTGGGC'
refProtSeq='MRSAAAAPRSPAVATESRRFAAARWPGWRSLQRPARRSGRGGGGAAPGPYPSAAPPPPGPGPPPSRQSSPPSASDCFGSNGNGGGAFRPGSRRLLGLGGPPRPFVVLLLPLASPGAPPAAPTRASPLGARASPPRSGVSLARPAPGCPRPACEPVYGPLTMSLKPQQQQQQQQQQQQQQQQQQQQQQQPPPAAANVRKPGGSGLLASPAAAPSPSSSSVSSSSATAPSSVVAATSGGGRPGLG'


# In[134]:

#creates an array with all reads from the data file
def stripLiner(seq):
    perSeqs = []
    a=open(seq,"r")
    for line in a:
        pSeq=line.strip()
        perSeqs.append(pSeq)
    return perSeqs
    #strip off newline character
    #rest of the logic to build a profile matrix
    #to loop through characters of every sequence, consider using enumerate


# In[135]:

#person 1 main run through
pSeqA = stripLiner("person1seqs")
print("Person 1:")
cols = range(len(pSeqA[0]))
DNAmatrix = pd.DataFrame(index=['A','C','G','T'],columns=cols)
DNAmatrix = DNAmatrix.fillna(0)
fillMatrix(pSeqA)
print("Profile Matrix:")
print(DNAmatrix)
tempListA = getSNPLoc(0.98,referenceDNA)
perA = printNewSeq(tempListA,referenceDNA)


# In[136]:

#person 2 main run through
pSeqB = stripLiner("person2seqs")
print("Person 2:")
cols = range(len(pSeqB[0]))
DNAmatrix = pd.DataFrame(index=['A','C','G','T'],columns=cols)
DNAmatrix = DNAmatrix.fillna(0)
fillMatrix(pSeqB)
print("Profile Matrix:")
print(DNAmatrix)
tempListB = getSNPLoc(0.98,referenceDNA)
perB = printNewSeq(tempList,referenceDNA)


# In[137]:

#person 3 main run through
pSeqC = stripLiner("person3seqs")
print("Person 3:")
cols = range(len(pSeqC[0]))
DNAmatrix = pd.DataFrame(index=['A','C','G','T'],columns=cols)
DNAmatrix = DNAmatrix.fillna(0)
fillMatrix(pSeqC)
print("Profile Matrix:")
print(DNAmatrix)
tempListC = getSNPLoc(0.98,referenceDNA)
perC = printNewSeq(tempListC,referenceDNA)


# In[139]:

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


# In[152]:

def convertProtein(seq):
    output = "";
    i = 0                  #starts index of string at 0
    while i < len(seq):    #while loop that stops if i is greater than the length of seq
        output += codontable[seq[i:i+3]] #adds the corresponding codon given the three base sequence to the output string
        i += 3             #increments i to to skip already processed codons
    return output


# In[153]:

#splices the protien seq to [start,stop)
def toProtien(sequence):
    output = []
    #goes through each converted read
    for seq in sequence:
        tempRNA = mRNAconvert(seq[start:stop-1])
        #print(tempRNA)
        protienSeq = convertProtein(tempRNA) #calls convert protein seq
        output.append(protienSeq)
    return output


# In[154]:

#calls on each splice method and prints them
protien1 = toProtien(perA)
protien2 = toProtien(perB)
protien3 = toProtien(perC)
print("Person 1 spliced protein sequence:")
print(protien1)
print("Person 2 spliced protein sequence:")
print(protien2)
print("Person 3 spliced protein sequence:")
print(protien3)


# In[155]:

#compares the refrence protein seq to the actual seq
def compareTo(listSeq):
    output = [] #creates var to store unequal seqs
    for seq in listSeq: #goes through seq in the list
        for i in range(len(seq)): # goes through each char
            if(not(seq[i] == refProtSeq[i])): 
                output.append([seq,i])
    return output


# In[157]:

print(compareTo(protien3))


# In[158]:

person3changes = compareTo(protien3)


# In[159]:

#returns the protein changes and the indicies of the changes in the DNA
def returnDNAchanges(list):
    for i in list:
        tempString = i[0]
        tempInt = i[1]
        print("Protien Sequence chnage from:", refProtSeq[tempInt],"to", tempString[tempInt],"the location is", tempInt)
        print("The DNA change goes from the codon:", referenceDNA[start+tempInt*3:start+tempInt*3+3],"and the SNP index in the persons DNA is between",  start+tempInt*3, "to",start+tempInt*3+3)


# In[160]:

returnDNAchanges(person3changes)


# In[138]:

#Some of the simplifications in the porblem set is that the reads have been aligned for us, and spliced from the whole genome,
#another simplication includes the gtf file already creates and us knowing the length of the read sequence. One huge simplification 
#was the only having to do one splice on the sequence. Finally some assumptions include the consitant number of reads and length of 
# the data file as well the same base error rate for each person. This also did not inlude the cost of each SNP allignemnt
#and the effect of the protein sequence.


# In[ ]:



