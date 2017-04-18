
# coding: utf-8

# In[164]:

#Jose Aguilar
#eid jra3528
v = "TCAGGGAACACCATACATACGATT"
w = "TTGGGAACACCAAACTACATTAGT"


# In[165]:

#output the table in a fromated table, does not work.
def printTable(table):
    output = "X w"
    for i in range(len(w)):
        output = output + " " + str(w[i])
    output = output + '\n'
    for row in range(len(table)):
        if(row == 0):
            output = output + "v"
        else:
            output = output + str(v[row-1])
        for col in range(len(table[0])):
            output = output + " " +str(table[row][col]);
        output = output + '\n'
    print(output)


# In[166]:

#Creates a blank matrix of specified size, and adds the indexed values from 0 to len(sting)+1 on the first row and col
#Blank Matrix for the edit distance of each base
matrixDist = [[0 for x in range(len(v)+1)] for y in range(len(w)+1)]
for i in range(len(v)+1):
    matrixDist[0][i] = i
for i in range(len(w)+1):
    matrixDist[i][0] = i


# In[167]:

#Creates a blank matrix with the same condtions as above, however it stores the edit opperation between each character of the string
#uses for loop to add the value of the basic indexes of a blank and non blank string
matrixEdit = [[0 for x in range(len(v)+1)] for y in range(len(w)+1)]
for i in range(len(v)+1):
    matrixEdit[0][i] = i
for i in range(len(w)+1):
    matrixEdit[i][0] = i


# In[168]:

#Method determines the lowest edit distance of the current row and col
def getMinOfEdit(row,col,cost):
    #DIS
    sub = matrixDist[row-1][col-1]+cost #detmines prevous edit distance of the cells before the row and col
    ins = matrixDist[row-1][col]+1  # by calling the sepecficed row and col change, above, left or diagonal
    delete = matrixDist[row][col-1]+1 # adds one to isertion/deletion score, add cost to the substition socrre based on match or mismatch
    #print("(",row,",",col,") ","D: ",delete,"I: ",ins,"S: ",sub)
    if(delete <= ins and delete <= sub): #gets the min value an returns value and character
        #print("DEL")
        return [delete, 'D']
    elif(ins <= sub and ins <= delete):
        #print("INS")
        return [ins, 'I']
    else:
        #print("SUB")
        return [sub, 'M']
    


# In[169]:

# Method caluclates the edit distance at each sepcified row and col, be comparing the values of each string at the mapped index
def editDist(row, col):
    if(v[col-1] == w[row-1]): #compares the string, subtracts one to match string index, if chars are equal
        costOfSub = 0         #set costOfSub to zero meaning subsitution should not be used MATCH
    else:
        costOfSub = 1        #MISMATCH of costOdSub
    temp = getMinOfEdit(row,col,costOfSub) #sets temp, and calls method of current row, col and costOfSub, returns array
    matrixDist[row][col] = temp[0] # sets 0th element of array to the edit distance
    matrixEdit[row][col] = temp[1] # sets 1st element of array to the edit operation
    return matrixDist[row][col] 


# In[170]:

# for loop goes through each row and col and calls edit dsit method
for row in range(1,len(w)+1):
    for col in range(1,len(v)+1):
        editDist(row,col)


# In[171]:

matrixDist


# In[172]:

matrixEdit


# In[173]:

#method starts from the lowest right cell and determines the most proabable prevous call of min
def backTrack(row,col):
    output = [] #creates empty array
    while(row != 0 and col !=0): #calls loop while row and col do not equal zero 
        temp = matrixEdit[row][col] # stores the current value of row and call (char)
        output.append(temp) # adds character to the array
        if(temp == 'M'): #determines which cell to go to next based on current char
            row = row-1
            col = col-1
        elif(temp == 'D'):
            col = col-1
        else:
            row = row-1
    return output #return array


# In[174]:

list = backTrack(24,24) #calls back track method and stores it 
list.reverse() # revereses the array to map to sequence


# In[175]:

print(list)


# In[176]:

#gets edited string from the backTrack list
def getEditOperations(string):
    output = [] #creates empty array
    temp = '' #creates empty char var
    index = 0 #sets index wIndex, vIndex and backtrack index to zerp
    vI = 0
    i = 0
    while(i < len(string) and index < len(w)): #loop that continues if the index of List, index of w string are below there sizes
        if(string[i] == 'M'):      # checks current char in list, and increments needed values and adds string to output 
            temp = 'M:'+ w[index]
            index += 1
            vI +=1
        elif(string[i] == 'I'):
            temp = 'I:'+ w[index]
            index += 1
            vI += 1
        elif(string[i] == 'D'):
            temp = 'D:'+ v[vI]
            vI += 1
            #index = index - 1
        i += 1
        #print(temp)
        output.append(temp)
    return output # retuns data


# In[177]:

def createCigar(data):
    data.append('X') #adds charcter to the end of the list to include last cigar score
    counter = 0 # sets counter to 0
    current = data[0] #sets the current char value of the string
    index = 0  #sets index of the list to zero
    output = '' #creates retuned string value
    while(index < len(data)): # loop if index is less then the length of the list of backTrack
        if(data[index] == current): #if the current char is equal to the char at the index increase counter
            counter += 1
        else:
            temp = '' + str(counter) + current #otherwise add the character and couter to the string
            #print(temp)
            output = output + temp #adds it to the string
            current = data[index] #sets current to char at index
            counter = 1 #resets counter
        index += 1 # increases index
    return output


# In[178]:

#prints everthing
print("Final Edit Distance:", matrixDist[len(w)][len(v)])
print("Edit operations:",getEditOperations(list))
print("Cigar string:", createCigar(list))


# In[ ]:




# In[ ]:



