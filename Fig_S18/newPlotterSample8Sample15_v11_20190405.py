import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
from scipy.stats import rankdata

#ver16 denovoClustering code onwards definition of clusterObject
class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same

def returnTestDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)
    seqDict={}
    for eachEntry in testDict:
        seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]
    sortedDictionary=dictSorter(seqDict)
    return [testDict,sortedDictionary]

inputSample='8'
outputSample='15'
[inputTest,inputSorted]=returnTestDict(inputSample)
[outputTest,outputSorted]=returnTestDict(outputSample)

def returnClusterDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_clusterDict_subrefCollapsed_v1.pckl','rb') as f:
        clusterDict=cPickle.load(f)
    return clusterDict

inputClusterDict=returnClusterDict(inputSample)
outputClusterDict=returnClusterDict(outputSample)

def returnFinalRefCollapsed(sampleName):
    with open(sampleName+'_finalRefCollapseOutput_v1_20190313.pckl','rb') as f:
        loadedList=cPickle.load(f)
    return loadedList

def returnFinalRefCollapsed2(sampleName):
    with open(sampleName+'_finalRefCollapseOutput_v2_20190318.pckl','rb') as f:
        loadedList=cPickle.load(f)
    return loadedList

[inputLinkingRefDict,inputCollapsedRefList]=returnFinalRefCollapsed(inputSample)
[outputLinkingRefDict,outputCollapsedRefList]=returnFinalRefCollapsed2(outputSample)

def ATHomopolymer(sequence):
    if 'C' in sequence or 'G' in sequence:
        return False
    else:
        return True

def isNotATDict(dict1):
    totalCounts=sum(dict1.values())
    totalATCount=0
    for eachSeq in dict1:
        if ATHomopolymer(eachSeq):
            totalATCount+=dict1[eachSeq]
    if (float(totalATCount)/float(totalCounts))>=0.8: #threshold for considering something an AT homopolymeric cluster
        return False
    else:
        return True
        
    
def returnCountsSeqVectors(clusterDict,collapsedRefDict,seqSorted,testDict):
    keysConsider=collapsedRefDict.keys()
    countsAssociated=[]
    seqAssociated=[]
    newKeys=[]
    for eachEl in keysConsider:
        totalCounts=0
        seqDictComb={}
        for eachTup in collapsedRefDict[eachEl]:
            totalCounts+=clusterDict[eachTup[0]][eachTup[1]].numberOfReads
            for eachIndex in clusterDict[eachTup[0]][eachTup[1]].listSeqIndices:
                seqDictComb[seqSorted[eachIndex]]=testDict[seqSorted[eachIndex]][sortWithPCRDuplicates]
        if isNotATDict(seqDictComb): 
            newKeys.append(eachEl)
            countsAssociated.append(totalCounts)       
            seqAssociated.append(seqDictComb) 
    return [newKeys,countsAssociated,seqAssociated]

[inputKeys,inputAssCounts,inputAssSeq]=returnCountsSeqVectors(inputClusterDict,inputLinkingRefDict,inputSorted,inputTest)
[outputKeys,outputAssCounts,outputAssSeq]=returnCountsSeqVectors(outputClusterDict,outputLinkingRefDict,outputSorted,outputTest)

with open('output_newPlotSample8Sample15_v10_20190330.pckl','rb') as f:
    [superClusterList,inputCrossContVal,outputCrossContVal]=cPickle.load(f)

inputCountVector=[]
outputCountVector=[]

def returnSumCounts(countList,refList,keysList):
    toReturn=0
    for eachEl in refList:
        toReturn+=countList[keysList.index(eachEl)]
    return toReturn

for eachEl in superClusterList:
    inputCount=returnSumCounts(inputAssCounts,eachEl[0],inputKeys)
    outputCount=returnSumCounts(outputAssCounts,eachEl[1],outputKeys)

    inputCountVector.append(inputCount)
    outputCountVector.append(outputCount)

#inputRank=rankdata(inputCountVector,method='ordinal')
#outputRank=rankdata(outputCountVector,method='ordinal')

inputTotalReadCount=631236 # totalReadsUsedUMINonduplicity
outputTotalReadCount=1001197 # totalReadsUsedUMINonduplicity

sortedIncreasingIndex=np.argsort(np.array(outputCountVector))
sortedDecreasingIndex=sortedIncreasingIndex[::-1]

inputColorVector=[]
outputColorVector=[]
colorsToUse=[np.array([0.0,88.0,0.0])/255.0,np.array([246.0,132.0,9.0])/255.0,np.array([0.0,70.0,149.0])/255.0,np.array([158.0,0.0,0.0])/255.0,np.array([185.0,79.0,211.0])/255.0]
colorsUsed=0
outputCountPlot=[]
inputCountPlot=[]
for counter, i in enumerate(sortedDecreasingIndex):
    if counter<5:
        outputColorVector.append(np.array(colorsToUse[colorsUsed]))
        inputColorVector.append(np.array([0.0,0.0,255.0])/255.0)
        colorsUsed+=1
    else:
        outputColorVector.append(np.array([125.0,125.0,0.0])/255.0)
        inputColorVector.append(np.array([0.0,0.0,255.0])/255.0)        
    outputCountPlot.append(outputCountVector[i])
    inputCountPlot.append(inputCountVector[i])
   
print 'Number of species plotted: '+str(len(inputCountPlot))
sortedIncreasingIndex1=np.argsort(np.array(inputCountPlot))
sortedDecreasingIndex1=sortedIncreasingIndex1[::-1]
finalInput=[]
finalOutput=[]
finalInputColor=[]
finalOutputColor=[]
for counter, i in enumerate(sortedDecreasingIndex1):
    finalInput.append(inputCountPlot[i])
    finalOutput.append(outputCountPlot[i])
    finalInputColor.append(inputColorVector[i])
    finalOutputColor.append(outputColorVector[i])

f,ax=plt.subplots(nrows=1,ncols=1)
f.set_size_inches((7,3.5))
for i in range(1,len(inputCountPlot)+1):
    ax.plot(i,finalInput[i-1],markerfacecolor=tuple(finalInputColor[i-1]),marker='o',markersize=4,markeredgecolor='none',linestyle='')
    ax.plot(i,finalOutput[i-1],markerfacecolor=tuple(finalOutputColor[i-1]),marker='o',markersize=4,markeredgecolor='none',linestyle='')

ax.plot(range(1,len(inputCountPlot)+1),finalInput,color=tuple(np.array([0.0,0.0,255.0])/255.0),linestyle='-',linewidth=0.1)
ax.plot(range(1,len(inputCountPlot)+1),finalOutput,color=tuple(np.array([125.0,125.0,0.0])/255.0),linestyle='-',linewidth=0.1)

ax.set_xlim(left=0,right=len(inputCountPlot)+1)
ax.set_ylim(bottom=-0.5, top=1000000)
ax.set_yscale('symlog',linthreshy=1.0)
ax.set_xlabel('Species')
ax.set_ylabel('Counts')
ax.set_xscale('symlog',linthreshx=50.0)

'''
plt.figure(2)
[n,bins,patches]=plt.hist(inputCrossContVal,bins=np.arange(0.0,1.05,0.05))
print n 
print bins

plt.figure(3)
[n,bins,patches]=plt.hist(outputCrossContVal,bins=np.arange(0.0,1.05,0.05))
print n 
print bins
'''
plt.savefig('outputPlot_newPlotterSample8Sample15_v11_20190405.pdf',dpi=300)
