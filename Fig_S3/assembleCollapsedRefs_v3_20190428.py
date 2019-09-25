from moduleUsefulFunctions_20180215 import *
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.drawing.nx_agraph import to_agraph
from scipy.interpolate import interp1d
import cPickle
import sys

#ver16 denovoClustering code onwards definition of clusterObject
class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

sortWithPCRDuplicates=0

def returnTestDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)
    seqDict={}
    for eachEntry in testDict:
        seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]
    sortedDictionary=dictSorter(seqDict)
    return [testDict,sortedDictionary,seqDict]

inputSample=sys.argv[1]
[inputTest,inputSorted,mergedTestDict]=returnTestDict(inputSample)

def returnClusterDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_clusterDict_subrefCollapsed_v1.pckl','rb') as f:
        clusterDict=cPickle.load(f)
    return clusterDict

inputClusterDict=returnClusterDict(inputSample)

def returnFinalRefCollapsed(sampleName):
    with open(sampleName+'_finalRefCollapseOutput_v1_20190313.pckl','rb') as f:
        loadedList=cPickle.load(f)
    return loadedList

[inputLinkingRefDict,inputCollapsedRefList]=returnFinalRefCollapsed(inputSample)

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
        newKeys.append(eachEl)
        countsAssociated.append(totalCounts)       
        seqAssociated.append(seqDictComb) 
    return [newKeys,countsAssociated,seqAssociated]

[inputKeys,inputAssCounts,inputAssSeq]=returnCountsSeqVectors(inputClusterDict,inputLinkingRefDict,inputSorted,inputTest)

def returnRep(seq):
    seqRep=0
    revCompRep=0
    maxSeqRep=0
    maxRevCompRep=0
    for eachEl in mergedTestDict: 
        if seq in eachEl:
            seqRep+=mergedTestDict[eachEl]
            if mergedTestDict[eachEl]>maxSeqRep:
                maxSeqRep=mergedTestDict[eachEl]            
        if revComp(seq) in eachEl:
            revCompRep+=mergedTestDict[eachEl]
            if mergedTestDict[eachEl]>maxRevCompRep:
                maxRevCompRep=mergedTestDict[eachEl]
    return [seqRep,revCompRep,maxSeqRep,maxRevCompRep] 

def returnExtSequence(graphG,startingNode,startingSeq,direction,growthThresh):    
    currentNode=startingNode
    currentSeq=startingSeq
    addedSeq=''
    while 1:
        weightList=[]
        repList=[]
        repSeqList=[]
        repRevCompList=[]
        nodeList=[]
        seqList=[]
        print 'Grow to the '+direction+' now...'
        if direction=='left':
            iteratorList=G.predecessors(currentNode)
        elif direction=='right':
            iteratorList=G.successors(currentNode)
        for eachNode in iteratorList:
            if direction=='left':
                tempSeq=eachNode[0]+currentSeq
                weightList.append(G[eachNode][currentNode]['weight'])
            elif direction=='right':
                tempSeq=currentSeq+eachNode[-1]
                weightList.append(G[currentNode][eachNode]['weight'])
            resultRep=returnRep(tempSeq)
            repList.append(min(resultRep[2:]))
            repSeqList.append(resultRep[0])
            repRevCompList.append(resultRep[1])
            nodeList.append(eachNode)
            seqList.append(tempSeq)
        if len(weightList)==0:
            bestVal=-1
        else:
            tempWeights=[]
            tempRepSeqList=[]
            tempRepRevCompList=[]
            tempNodeList=[]
            tempSeqList=[]
            tempRepList=[]
            indexSorted=np.argsort(np.array(repList))
            bestRepVal=repList[indexSorted[-1]]
            for counter, eachEl in enumerate(repList):
                if eachEl==bestRepVal: 
                    tempWeights.append(weightList[counter])
                    tempRepSeqList.append(repSeqList[counter])
                    tempRepRevCompList.append(repRevCompList[counter])
                    tempNodeList.append(nodeList[counter])
                    tempSeqList.append(seqList[counter])
                    tempRepList.append(min([repSeqList[counter],repRevCompList[counter]]))
            indexSorted1=np.argsort(np.array(tempRepList))
            bestVal=tempRepList[indexSorted1[-1]]                               
        if bestVal<growthThresh:
            return addedSeq
        else:
            tempWeights2=[]
            tempRepSeqList2=[]
            tempRepRevCompList2=[]
            tempNodeList2=[]
            tempSeqList2=[]
            for counter, eachEl in enumerate(tempRepList):
                if eachEl==bestVal:
                    tempWeights2.append(tempWeights[counter])
                    tempRepSeqList2.append(tempRepSeqList[counter])
                    tempRepRevCompList2.append(tempRepRevCompList[counter])
                    tempNodeList2.append(tempNodeList[counter])
                    tempSeqList2.append(tempSeqList[counter])
            indexSorted2=np.argsort(np.array(tempWeights2))
            bestIndex=indexSorted2[-1]
            print tempSeqList2[bestIndex]+'\t'+str(tempRepSeqList2[bestIndex])+'\t'+str(tempRepRevCompList2[bestIndex])
            currentSeq=tempSeqList2[bestIndex]
            currentNode=tempNodeList2[bestIndex]
            if direction=='left':
                addedSeq=currentNode[0]+addedSeq
            elif direction=='right':
                addedSeq+=currentNode[-1]

def returnGraph(dictSequences,kmer):
    k=kmer
    seqDict=dictSequences
    G=nx.DiGraph() #make an empty directed graph, may not be ACYCLIC!
    for counter, eachSeq in enumerate(seqDict.keys()):
        if len(eachSeq)>=k:
            for i in range(0,len(eachSeq)+1-k,1):
                oneNode=eachSeq[i:i+k-1]
                secondNode=eachSeq[i+1:i+k]
                if G.has_edge(oneNode,secondNode):
                    G[oneNode][secondNode]['weight']+=seqDict[eachSeq]
                else:
                    G.add_edge(oneNode,secondNode,weight=seqDict[eachSeq])
    return G

def returnNumberMutations(alignInfoList):
    numberMutations=0
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            startPosition=eachEl.position
        elif eachEl.typeOfEvent=='End':
            endPosition=eachEl.position
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion':
            if eachEl.position==startPosition or eachEl.position==endPosition: #only count insertions in the middle of the sequence as mutations
                pass
            else:
                numberMutations+=1 #as insertions are not being collapsed together (collapseInsertions=False in the smithWaterman calls below), each base inserted is counted separately
        elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
            numberMutations+=1
    return numberMutations

def returnStartPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            return eachEl.position

def returnEndPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='End':
            return eachEl.position

def FivePrimeClipSequence(alignInfoList):
    startPosition=returnStartPosition(alignInfoList)
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' and eachEl.position==startPosition:
            return eachEl.notes
    return ''

def ThreePrimeClipSequence(alignInfoList):
    endPosition=returnEndPosition(alignInfoList)
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Insertion' and eachEl.position==endPosition:
            return eachEl.notes
    return ''  

def collapseRefs(listSeqToAlign,listSeqToWhichAlign,par1,par2,par3):
    toleranceNumberMutations=par1
    ntLengthUpdateFullLength=par2
    softClippingThresh=par3
    indexListToReturn=[]
    for eachRef in listSeqToAlign:
        alignmentInfo=[]
        alignmentScore=[]
        for eachRef1 in listSeqToWhichAlign:
            score=-1
            listRefs=[eachRef1,revComp(eachRef1)]
            for eachSeq in listRefs:
                alignOut=smithWaterman(eachSeq,eachRef,collapseInsertions=False)
                if alignOut[0]>score:
                    score=alignOut[0]
                    info=alignedObject(alignOut,eachRef1,-1) #note eachRef is specified as refUsed, not eachSeq; artificial count of -1 passed
            alignmentInfo.append(info)
            alignmentScore.append(float(score)/min([len(eachRef1),len(eachRef)]))
        alignmentScore=np.array(alignmentScore)
        indicesSort=np.argsort(alignmentScore)        
        subclusterIndex=indicesSort[-1] #note there may be clashes between two references. Here, whatever sorts to the -1 position is chosen
        numberMutations=returnNumberMutations(alignmentInfo[subclusterIndex].alignmentInfo[1])        
        endPositionBestAlignment=returnEndPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
        startPositionBestAlignment=returnStartPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
        maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(alignmentInfo[subclusterIndex].refUsed),len(eachRef)]))
        totalNumberMutationsCalc=len(FivePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+len(ThreePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+numberMutations
        if numberMutations<=np.ceil(toleranceNumberMutations*(endPositionBestAlignment-startPositionBestAlignment)):
            if (totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength)) or (startPositionBestAlignment<softClippingThresh and endPositionBestAlignment>(len(alignmentInfo[subclusterIndex].refUsed)-softClippingThresh)):
                indexListToReturn.append(subclusterIndex)
            else:
                indexListToReturn.append('Left')
        else:
            indexListToReturn.append('Left')
    return indexListToReturn

def returnSumCountsAcrossDict(listDict):
    toReturnList=[]
    for eachEl in listDict:
        toReturnList.append(sum(eachEl.values()))
    return toReturnList

finalCountSortToUse=np.argsort(inputAssCounts)
finalCountSortToUse=finalCountSortToUse[::-1] #descending order
listRefsToSave=[]
toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
softClippingThresh=2 #for length upgrading
countsToSave=[]
dictSequencesToSave=[]
errorCounter=0
for eachI in finalCountSortToUse:
    k=20
    print inputKeys[eachI]+'\t'+str(inputAssCounts[eachI])
    seqDict=inputAssSeq[eachI]    
    G=returnGraph(seqDict,k)
    print k
    maxWeight=0
    for eachEl in G.edges():
        if G[eachEl[0]][eachEl[1]]['weight']>maxWeight:
            maxWeight=G[eachEl[0]][eachEl[1]]['weight']
            maxEdge=(eachEl[0],eachEl[1])
            maxSeq=eachEl[0]+eachEl[1][-1]

    print maxEdge
    print maxWeight
    print maxSeq

    growthThreshold=float(maxWeight)/10.0 #counts

    leftAdded=returnExtSequence(G,maxEdge[0],maxSeq,'left',growthThreshold)
    rightAdded=returnExtSequence(G,maxEdge[1],leftAdded+maxSeq,'right',growthThreshold)

    currentContig=leftAdded+maxSeq+rightAdded

    newSeqDict={}
    moreToLeft=[]
    for eachEl in mergedTestDict: 
        if currentContig in eachEl or revComp(currentContig) in eachEl:
            newSeqDict[eachEl]=mergedTestDict[eachEl]
        if currentContig in eachEl:
            if eachEl.find(currentContig)>=(len(eachEl)-len(currentContig)-eachEl.find(currentContig)):
                moreToLeft.append(True)
            else:
                moreToLeft.append(False)

    k=len(currentContig)
    G=returnGraph(newSeqDict,k)

    maxEdge=(currentContig[:-1],currentContig[1:])
    maxSeq=currentContig
    newGrowthThreshold=5

    #Decide which direction to grow in first based on location of currentContig in reads
    if 2*sum(moreToLeft)>len(moreToLeft): #more bases are on left, extend left first
        leftAdded=returnExtSequence(G,maxEdge[0],maxSeq,'left',newGrowthThreshold)
        rightAdded=returnExtSequence(G,maxEdge[1],leftAdded+maxSeq,'right',newGrowthThreshold)    
    else:
        rightAdded=returnExtSequence(G,maxEdge[1],maxSeq,'right',newGrowthThreshold)    
        leftAdded=returnExtSequence(G,maxEdge[0],maxSeq+rightAdded,'left',newGrowthThreshold)

    finalRefToTest=leftAdded+maxSeq+rightAdded
    if len(listRefsToSave)==0:
        listRefsToSave.append(finalRefToTest)
        countsToSave.append(inputAssCounts[eachI])
        dictSequencesToSave.append(inputAssSeq[eachI])
    else:
        alignmentMappingList=collapseRefs([finalRefToTest],listRefsToSave,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
        alignmentMappingIndex=alignmentMappingList[0]
        if alignmentMappingIndex=='Left':
            listRefsToSave.append(finalRefToTest)
            countsToSave.append(inputAssCounts[eachI])
            dictSequencesToSave.append(inputAssSeq[eachI])
        else:
            countsToSave[alignmentMappingIndex]+=inputAssCounts[eachI]
            for eachNovelSeq in inputAssSeq[eachI]:
                if eachNovelSeq in dictSequencesToSave[alignmentMappingIndex]:
                    dictSequencesToSave[alignmentMappingIndex][eachNovelSeq]+=inputAssSeq[eachI][eachNovelSeq]
                    errorCounter+=1
                else:
                    dictSequencesToSave[alignmentMappingIndex][eachNovelSeq]=inputAssSeq[eachI][eachNovelSeq]
                    
    print ''
    print ''
    print ''
    print ''
    print listRefsToSave
    print countsToSave
    print returnSumCountsAcrossDict(dictSequencesToSave)
    print 'New cluster beginning...'

print 'Error counter: '+str(errorCounter)
with open(inputSample+'_assembleCollapsedRefs_v3_20190428.pckl','wb') as f:
    cPickle.dump([listRefsToSave,countsToSave,dictSequencesToSave],f,protocol=cPickle.HIGHEST_PROTOCOL)

