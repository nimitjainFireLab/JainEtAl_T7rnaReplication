import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess

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

#ver16 denovoClustering code onwards definition of clusterObject
class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

def diNucAtEnds(seq,dotBrack,diNucSeq,distEnds):
    fiveIndex=dotBrack.find('(')
    seqFive=seq[:fiveIndex]
    seqFive=seqFive[:distEnds]
    
    threeIndex=dotBrack.rfind(')')+1
    seqThree=seq[threeIndex:]
    seqThree=seqThree[::-1]
    seqThree=seqThree[:distEnds]
    
    if diNucSeq in seqFive and diNucSeq in seqThree:
        return True
    else:
        return False


def getBestRef(lengthV,indexV,countV,diGV,diCV):
    dictRefs={}
    indexUse=indexV[lengthV]
    countUse=countV[lengthV]
    diGUse=diGV[lengthV]
    diCUse=diCV[lengthV]

    diGIndex=indexUse[diGUse]
    diGCount=countUse[diGUse]

    diCIndex=indexUse[diCUse]
    diCCount=countUse[diCUse]

    for counter, eachEl in enumerate(diGIndex):
        listSeq=getAllCombs(sortedDictionary[eachEl],dotBracketVector[eachEl],'GG',distanceFromEnds)
        for eachSeq in listSeq:
            if eachSeq in dictRefs:
                dictRefs[eachSeq][0]+=diGCount[counter]
            else:
                dictRefs[eachSeq]=[diGCount[counter],0]
                
    for counter, eachEl in enumerate(diCIndex):
        listSeq=getAllCombs(sortedDictionary[eachEl],dotBracketVector[eachEl],'CC',distanceFromEnds)
        for eachSeq in listSeq:
            if revComp(eachSeq) in dictRefs:
                dictRefs[revComp(eachSeq)][1]+=diCCount[counter]

    newDict={}
    for eachSeq in dictRefs:
        newDict[eachSeq]=min(dictRefs[eachSeq])

    sortedNewDict=dictSorter(newDict)
    if newDict[sortedNewDict[0]]==0:
        return ''
    else:
        return sortedNewDict[0]
    
def getAllCombs(seq,dotBrack,diNucSeq,distEnds):
    fiveIndex=dotBrack.find('(')
    seqFive=seq[:fiveIndex]
    seqFive=seqFive[:distEnds]
    listIndices5=[]
    for i in range(len(seqFive)+1-len(diNucSeq)):
        if seqFive[i:i+len(diNucSeq)]==diNucSeq:
            listIndices5.append(i)
            
    threeIndex=dotBrack.rfind(')')+1
    seqThree=seq[threeIndex:]
    seqThree=seqThree[::-1]
    seqThree=seqThree[:distEnds]
    diNucInv=diNucSeq[::-1]
    listIndices3=[]
    for i in range(len(seqThree)+1-len(diNucSeq)):
        if seqThree[i:i+len(diNucSeq)]==diNucInv:
            listIndices3.append(len(seq)-i)

    listAllCombs=[]
    for i in listIndices5:
        for j in listIndices3:
            listAllCombs.append(seq[i:j])
    return listAllCombs

sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same

testDictName=sys.argv[1]
with open(testDictName,'rb') as f:
    testDict=cPickle.load(f)

sampleName=testDictName
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_testDict.pckl')]

seqDict={}
for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)

with open(sampleName+'_clusterDict_subrefCollapsed_v1.pckl','rb') as f:
    clusterDict=cPickle.load(f)

with open(sampleName+'_combinedList_structureMetricsPreCalc_v2.pckl','rb') as f:
    combinedList=cPickle.load(f)

[seqCount,oneHairpin,seqLength,diGReq,diCReq,dotBracketVector,isATstretch,hairpinScoreVector]=combinedList

#sortedDictionary list should be in same order as all these lists
seqCount=np.array(seqCount)
oneHairpin=np.array(oneHairpin)
seqLength=np.array(seqLength)
diGReq=np.array(diGReq)
diCReq=np.array(diCReq)
isATstretch=np.array(isATstretch)
hairpinScoreVector=np.array(hairpinScoreVector)
onlyDiGOrDiCReq=np.logical_xor(diGReq,diCReq)
onlyDiGReq=np.logical_and(diGReq,onlyDiGOrDiCReq)
onlyDiCReq=np.logical_and(diCReq,onlyDiGOrDiCReq)

distanceFromEnds=5 #in nucleotides
newOnlyDiGReq=[]
for counter, eachEl in enumerate(onlyDiGReq):
    if eachEl==False:
        newOnlyDiGReq.append(False)
    else:
        newOnlyDiGReq.append(diNucAtEnds(sortedDictionary[counter],dotBracketVector[counter],'GG',distanceFromEnds))
newOnlyDiGReq=np.array(newOnlyDiGReq)
onlyDiGReq=newOnlyDiGReq

newOnlyDiCReq=[]
for counter, eachEl in enumerate(onlyDiCReq):
    if eachEl==False:
        newOnlyDiCReq.append(False)
    else:
        newOnlyDiCReq.append(diNucAtEnds(sortedDictionary[counter],dotBracketVector[counter],'CC',distanceFromEnds))
newOnlyDiCReq=np.array(newOnlyDiCReq)
onlyDiCReq=newOnlyDiCReq

slidingWindowLength=10 #in nucleotides
xToPlot=[]
yToPlot=[]
listAllRefs=[]
listAllClusters=[]
listdiGRef=[]
paramDistinguishWithinRef=4
maxNumberMutationsAllowed=2

for counter, eachCluster in enumerate(clusterDict):
    for eachRef in clusterDict[eachCluster]:
        if clusterDict[eachCluster][eachRef].forAlignment==1: #subref collapsed dicts should have all subreferences set to forAlignment=1 anyway
            arrayIndices=np.array(clusterDict[eachCluster][eachRef].listSeqIndices)
            countsToUse=seqCount[arrayIndices]
            lengthsToUse=seqLength[arrayIndices]
            diGToUse=onlyDiGReq[arrayIndices]
            diCToUse=onlyDiCReq[arrayIndices]
            #note that don't need to use oneHairpin vector because currently by definition, if oneHairpin is False, diG and diC vectors would be false
            
            readNumberSubcluster=np.sum(countsToUse)
            xToPlot.append(readNumberSubcluster)
            listAllRefs.append(eachRef)
            listAllClusters.append(eachCluster)
 
            minLength=np.min(lengthsToUse)
            maxLength=np.max(lengthsToUse)
            currentWindowLow=minLength
            currentWindowHigh=minLength+slidingWindowLength
            toBreakLoop=0
            bestIndexVal=-1
            while toBreakLoop==0:
                if currentWindowHigh>maxLength:
                    toBreakLoop=1
                logicalIndexingArray=np.logical_and(lengthsToUse>=currentWindowLow,lengthsToUse<currentWindowHigh)
                countsWindow=countsToUse[logicalIndexingArray]
                diGWindow=diGToUse[logicalIndexingArray]
                diCWindow=diCToUse[logicalIndexingArray]
                tempIndex=min([np.sum(countsWindow[diGWindow]),np.sum(countsWindow[diCWindow])])
                if tempIndex>=bestIndexVal: #note that by using >= in this, in case of ties for bestIndexVal--the longest set of sequences will be used
                    bestIndexVal=tempIndex
                    bestLengthArray=logicalIndexingArray
                currentWindowLow+=1
                currentWindowHigh+=1
            if bestIndexVal<0.5: #really looking for bestIndexVal being 0, comparing to 0.5 in case float comparisons cause problems
                yToPlot.append(0)
                listdiGRef.append('')
            else:
                #try and establish a best reference for the species
                bestRef=getBestRef(bestLengthArray,arrayIndices,countsToUse,diGToUse,diCToUse)
                if bestRef=='':
                    yToPlot.append(0)
                    listdiGRef.append('')
                else:
                    #have a best ref now, compute bestIndexVal through alignment
                    listdiGRef.append(bestRef)
                    print bestRef
                    refList=[bestRef,revComp(bestRef)]
                    numberStoreDict={}
                    for eachEl in refList:
                        numberStoreDict[eachEl]=0.0 
                    for eachIndexVal in arrayIndices[bestLengthArray]:
                        eachEl=sortedDictionary[eachIndexVal]
                        #some of code below is taken from aligner_257Dimer.py
                        alignmentInfo=[]
                        alignmentScore=[]
                        for eachSeq in refList:
                            alignOut=smithWaterman(eachSeq,eachEl,collapseInsertions=False)
                            alignmentScore.append(alignOut[0]) #don't need length normalization as all references have same length!
                            alignmentInfo.append(alignedObject(alignOut,eachSeq,seqCount[eachIndexVal]))
                        alignmentScore=np.array(alignmentScore)
                        indicesSort=np.argsort(alignmentScore)
                        bestAlignment=indicesSort[-1]
                        numberMutations=returnNumberMutations(alignmentInfo[bestAlignment].alignmentInfo[1])
                        endPositionBestAlignment=returnEndPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
                        startPositionBestAlignment=returnStartPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
                        totalNumberMutationsCalc=startPositionBestAlignment+(len(alignmentInfo[bestAlignment].refUsed)-endPositionBestAlignment)+numberMutations
                        if totalNumberMutationsCalc<=maxNumberMutationsAllowed:
                            if alignmentScore[bestAlignment]>=alignmentScore[indicesSort[-2]]+paramDistinguishWithinRef:                    
                                numberStoreDict[alignmentInfo[bestAlignment].refUsed]+=seqCount[eachIndexVal]
                    print numberStoreDict
                    yToPlot.append(min(numberStoreDict.values()))

with open(sampleName+'_bestIndexValCalculatorV2Code_20190218_allLists.pckl','wb') as f:
    cPickle.dump([xToPlot,yToPlot,listAllRefs,listAllClusters,listdiGRef],f,protocol=cPickle.HIGHEST_PROTOCOL)


