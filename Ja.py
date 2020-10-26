def JaccardSimilarity(a,b):
	intersection=len(set(a).intersection(set(b)))
	union=len(a)+len(b)-intersection
	output= float(intersection)/float(union)
	return output

import random
import mmh3
import time

############################################################
##Problem 3.1
##just to test the function JaccardSimilarity
#listof10=[[1,3,5,7,9],[1,2,3,4,5],[3,2,3,4,11],[2,2,3,4,23],[1,3,3,4,2],[1,2,23,4,5],[1,2,6,4,5],[1,2,3,7,5],[1,2,413,4,8],[6,2,8,5,7]]
# temp=[]
# for i in range(0,len(listof10)-1):
# 	for j in range(i+1,len(listof10)):
# 		print(JaccardSimilarity(listof10[i],listof10[j]))
# print("--- %s seconds ---" % (time.time() - start_time))
#(!Need to be checked!)took 0.1705s for above calculation of C(10,2)=45 pairs. We expect 0.1705/45=0.00379 for one pair, so 0.00379*C(N,2) seconds=O(N^2).
############################################################


##########################################################################
#import the data.txt file
file = open('data.txt', 'r')
lines=file.readlines()
lines = [x.strip().split() for x in lines]#convert the input from string to 2-nested list of substrings
##############################################################################

#######################################################################
##These are pairs whose JS is acctually >0.85 after inspection
# print(lines[2242])
# print(lines[49824])
# print(JaccardSimilarity(lines[2242],lines[49824]))
# print(lines[77683])
# print(lines[124849])
# print(JaccardSimilarity(lines[77683],lines[124849]))
# print(lines[215149])
# print(lines[274946])
# print(JaccardSimilarity(lines[215149],lines[274946]))
# print(lines[167121])
# print(lines[479963])
# print(JaccardSimilarity(lines[167121],lines[479963]))
# print(lines[318540])
# print(lines[752337])
# print(JaccardSimilarity(lines[318540],lines[752337]))
######################################################################


##########################################################################
##Problem 3.2
##This is to take a random sample from the dataset and get a rough estimate about the distribution
##We observed that the distribution of JS is around 0.05, and this information is useful determining the parameters later on
#sample=random.sample(lines, 5000)
#file1 = open("sample_distribution_test.txt","w")#write the distribution of jaccardsimiarity of our random sample
# for i in range(len(sample)-1):
# 	for j in range(i+1,len(sample)):
# 		file1.write(str(JaccardSimilarity(sample[i],sample[j]))+"\n")
#########################################################################


def minHash(mylist,i):
	temp=[]
	for element in mylist:
		temp.append(mmh3.hash(element,i))#using ith seed in the family of minhash
	index=temp.index(min(temp))#get the index of the element with the smallest hashing value
	return mylist[index]

def getSignature(mylist,k):
	signature_list=[]
	for i in range(k):
		signature_list.append(minHash(mylist,i))
	return signature_list

def creatBands(mylist,r,b):#to divide the signatures into bands based on the value of r and b
	temp=[]
	for i in range(b):
		row=[]
		for j in range(r):
			row.append(mylist[i*r+j])
		temp.append(row)
	return temp

start_time = time.time()

#############################################################################
#We choose k=20, r=4, =5, because of the following reasons. Firstly, k=20 is reasonably
#small. It took about 2-3 minutes to get the signatures via minhash and 10-15 minutes for
#comparing the bands and selecting potential pairs. Although the candidates pool is of size
#20,000,000+, it only took 5 minutes to inspect the pool and extract actual required pairs.
#Secondly, we analyzed the funciton 1-(1-s^r)^b. It reached the steepest point at around 
#(1/b)^(1/r)=0.66874, which is good enough for us to filter out many pairs whose JS are
#distributed around 0.05 based on our observation from the sample taken previously. Also,
#the funciton reached about 0.965 when s is about 0.84 so this is ideal for us to include
#our actual pairs.
k=20
r=4
b=5
#############################################################################



signature_list=[getSignature(element,k) for element in lines]#create signatures for each element in lines
t0=time.time() - start_time
print("until getSignature using minHash: "+"--- %s seconds ---" % (t0))
bands_list=[creatBands(element,r,b) for element in signature_list]#create bands for each signature list


totalDict = {} # key is the hash value of a string converted from the numbers
# in a band, value is a list containing the index of set where the hash value is
# found
pairList = [] #candidate list
file2 = open("potential_pairs.txt","w")#to save our potential list

for i in range(len(bands_list)):
	if i%10000==0:
		print("progress:"+" "+str(i))#just to monitor the progress of hashing bands
	for j in range(len(bands_list[i])):#iterate the bands list
		bandString = str(bands_list[i][j])
		hashedString = mmh3.hash(bandString, j)#hash a band to a dictionary with the hash value as the key
		if hashedString in totalDict:
			for index in totalDict[hashedString]:
				file2.write(str(index)+','+str(i)+"\n")#those pairs fall into the same bucket are our potential pairs because one of their bands has the same hashing value
				pairList.append([index,i])#the value in the dictionary is a list of indexes whose bands have the hash value as the key
			totalDict[hashedString].append(i)
		else:
			totalDict[hashedString] = [i]
print("lenth of potential pair list is "+str(len(pairList)))
print("until selecting potential_pairs using: "+"--- %s seconds ---" % (time.time() - start_time))

###############################################################
#Write up actual pairs satisfying our threshhold after inspected all candidates
file3 = open("actual_pairs.txt","w")
for pair in pairList:#here the pairList might have duplicates because two elements can have more than one band matched and extracted by us
	if JaccardSimilarity(lines[pair[0]],lines[pair[1]])>0.85:
		file3.write(str(pair)+"\n")
		print(pair)
###############################################################

print("total time: ""--- %s seconds ---" % (time.time() - start_time))
