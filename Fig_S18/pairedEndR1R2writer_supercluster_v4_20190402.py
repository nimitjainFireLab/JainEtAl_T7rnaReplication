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

inputRefToProbe='TTTCCTATTTAAATCAATTTTTGGG'

allSeqDict={}
for eachEl in superClusterList:
    if inputRefToProbe in eachEl[0]: #work with this eachEl
        for eachEl1 in eachEl[0]:
            for eachSeq in inputAssSeq[inputKeys.index(eachEl1)]:
                if eachSeq in allSeqDict: #should never be launched
                    #print 'What??'
                    allSeqDict[eachSeq]+=inputAssSeq[inputKeys.index(eachEl1)][eachSeq]
                else:
                    allSeqDict[eachSeq]=inputAssSeq[inputKeys.index(eachEl1)][eachSeq]
        for eachEl1 in eachEl[1]:
            for eachSeq in outputAssSeq[outputKeys.index(eachEl1)]:
                if eachSeq in allSeqDict:
                    #print 'Expected'
                    allSeqDict[eachSeq]+=outputAssSeq[outputKeys.index(eachEl1)][eachSeq]
                else:
                    allSeqDict[eachSeq]=outputAssSeq[outputKeys.index(eachEl1)][eachSeq]
'''
for eachSeq in dictSorter(allSeqDict):
    print eachSeq+'\t'+str(allSeqDict[eachSeq])
'''

with open('dummyR1_supercluster_v4_20190402.fastq','w') as f1, open('dummyR2_supercluster_v4_20190402.fastq','w') as f2:
    numberSeqPrinted=1
    for eachSeq in dictSorter(allSeqDict):
        for i in range(allSeqDict[eachSeq]):
            f1.write('@'+str(numberSeqPrinted)+'\n')
            f2.write('@'+str(numberSeqPrinted)+'\n')

            f1.write(eachSeq+'\n')
            f2.write(revComp(eachSeq)+'\n')
            
            f1.write('+\n')
            f2.write('+\n')

            f1.write('K'*len(eachSeq)+'\n')
            f2.write('K'*len(eachSeq)+'\n')
            numberSeqPrinted+=1

