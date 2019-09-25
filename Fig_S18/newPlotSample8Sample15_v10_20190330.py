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

def findAssociatedRefs(keysList1,refList,seqDictList1,keysList2,seqDictList2):
    newRefs={}
    for counter, eachEl in enumerate(refList):
        dictOfInterest=seqDictList1[keysList1.index(eachEl)]
        for counter1, eachEl1 in enumerate(keysList2):
            for eachSeq in seqDictList2[counter1]:
                if eachSeq in dictOfInterest:
                    newRefs[eachEl1]=0
                    break
    return newRefs.keys()

inputRefsAlreadySeen={}
outputRefsAlreadySeen={}
numberSuperClusters=0
superClusterList=[]
for counter, eachEl in enumerate(inputKeys):
    print counter
    if eachEl not in inputRefsAlreadySeen: #is a collapsed ref
        toggle=0
        superClusterList.append([[],[]]) #first list is for input reference/clusters, second list is for output reference/clusters
        associatedInputRefs=[eachEl]
        superClusterList[numberSuperClusters][0].append(eachEl)
        inputRefsAlreadySeen[eachEl]=0
        while toggle!=2:
            if toggle==0:
                associatedOutputRefs=findAssociatedRefs(inputKeys,associatedInputRefs,inputAssSeq,outputKeys,outputAssSeq)
                newRefs={}
                for counter1, eachEl1 in enumerate(associatedOutputRefs):
                    if eachEl1 not in outputRefsAlreadySeen:
                        newRefs[eachEl1]=0
                if len(newRefs)==0:
                    toggle=2
                else:
                    for eachEl1 in newRefs:
                        outputRefsAlreadySeen[eachEl1]=0
                        superClusterList[numberSuperClusters][1].append(eachEl1)
                    toggle=1
            elif toggle==1:
                associatedInputRefs=findAssociatedRefs(outputKeys,associatedOutputRefs,outputAssSeq,inputKeys,inputAssSeq)
                newRefs={}
                for counter1, eachEl1 in enumerate(associatedInputRefs):
                    if eachEl1 not in inputRefsAlreadySeen:
                        newRefs[eachEl1]=0
                if len(newRefs)==0:
                    toggle=2
                else:
                    for eachEl1 in newRefs:
                        inputRefsAlreadySeen[eachEl1]=0
                        superClusterList[numberSuperClusters][0].append(eachEl1)
                    toggle=0
        numberSuperClusters+=1
print numberSuperClusters
#print superClusterList
crossContaminationThreshold=0.8
indexInputCrossMatrix=7
indexOutputCrossMatrix=13

def returnCombinedListDictSeq(listSuperClusters,indexReturn,keysList,seqDictList):
    listToReturn=[]
    for counter, eachEl in enumerate(listSuperClusters):
        tempDict={}
        for counter1, eachEl1 in enumerate(eachEl[indexReturn]):
            dictToFill=seqDictList[keysList.index(eachEl1)]
            for eachSeq in dictToFill:
                tempDict[eachSeq]=0
        listToReturn.append(tempDict)
    return listToReturn

inputCombinedListDictSeq=returnCombinedListDictSeq(superClusterList,0,inputKeys,inputAssSeq)
outputCombinedListDictSeq=returnCombinedListDictSeq(superClusterList,1,outputKeys,outputAssSeq)

def returnCrossContMatrix(listLibraries,combinedListDict):
    resultsMatrix=[]
    for counter, i in enumerate(listLibraries):
        print i
        with open(str(i)+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
            toCheckDict=cPickle.load(f)
        resultsMatrix.append([])
        for counter1, eachEl in enumerate(combinedListDict):
            totalCount=0.0
            for eachSeq in eachEl:
                if eachSeq in toCheckDict:
                    totalCount+=toCheckDict[eachSeq][sortWithPCRDuplicates]
            resultsMatrix[counter].append(totalCount)
    resultsMatrix=np.array(resultsMatrix)
    return resultsMatrix

inputCrossContMatrix=returnCrossContMatrix([1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18],inputCombinedListDictSeq)
outputCrossContMatrix=returnCrossContMatrix([1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18],outputCombinedListDictSeq)


def checkWhetherNotCont(resultsMatrix,idName,refID):
    if resultsMatrix[idName,refID]==0:
        crossContVal=1 #this would basically correspond to cases where the output has no sequences corresponding to the input
    else:
        crossContVal=float(resultsMatrix[idName,refID])/float(np.sum(resultsMatrix[:,refID]))
    if crossContVal>crossContaminationThreshold:
        return [True,crossContVal]
    else:
        return [False,crossContVal]

inputCrossContVector=[]
outputCrossContVector=[]
inputCrossContVectorVal=[]
outputCrossContVectorVal=[]

for counter, eachEl in enumerate(superClusterList):
    inputCrossContVector.append(checkWhetherNotCont(inputCrossContMatrix,indexInputCrossMatrix,counter)[0])
    outputCrossContVector.append(checkWhetherNotCont(outputCrossContMatrix,indexOutputCrossMatrix,counter)[0])
    inputCrossContVectorVal.append(checkWhetherNotCont(inputCrossContMatrix,indexInputCrossMatrix,counter)[1])
    outputCrossContVectorVal.append(checkWhetherNotCont(outputCrossContMatrix,indexOutputCrossMatrix,counter)[1])

newSuperClusterList=[]
for counter, eachEl in enumerate(superClusterList):
    if inputCrossContVector[counter]==True and outputCrossContVector[counter]==True:        
        newSuperClusterList.append(eachEl)

with open('output_newPlotSample8Sample15_v10_20190330.pckl','wb') as f:
    cPickle.dump([newSuperClusterList,inputCrossContVectorVal,outputCrossContVectorVal],f,protocol=cPickle.HIGHEST_PROTOCOL)
