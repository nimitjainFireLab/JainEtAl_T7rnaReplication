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
toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
softClippingThresh=2 #for length upgrading
countThreshold=0

inputSample=sys.argv[1]

def returnClusterDict(sampleName):
    listAllRefs1=[]
    listAllClusters1=[]
    listTotalCounts=[]
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_clusterDict_subrefCollapsed_v1.pckl','rb') as f:
        clusterDict=cPickle.load(f)
    for counter, eachCluster in enumerate(clusterDict):
        for eachRef in clusterDict[eachCluster]:
            if clusterDict[eachCluster][eachRef].forAlignment==1: #subref collapsed dicts should have all subreferences set to forAlignment=1 anyway
                if clusterDict[eachCluster][eachRef].numberOfReads>=countThreshold:
                    listAllRefs1.append(eachRef)
                    listAllClusters1.append(eachCluster)
                    listTotalCounts.append(clusterDict[eachCluster][eachRef].numberOfReads)
    return [clusterDict,listAllRefs1,listAllClusters1,listTotalCounts]

[inputClusterDict,newClusterDictRef,newClusterDictCluster,newTotalCount]=returnClusterDict(inputSample)

print len(newClusterDictRef)

sortedIndices=np.argsort(np.array(newTotalCount))
sortedIndices=sortedIndices[::-1] #decreasing order

newClusterDictRef2=[]
newClusterDictRef2Dict={}

for counter, eachEl in enumerate(sortedIndices):
    print counter
    if counter==0:
        newClusterDictRef2.append(newClusterDictRef[eachEl])
        newClusterDictRef2Dict[newClusterDictRef[eachEl]]={(newClusterDictCluster[eachEl],newClusterDictRef[eachEl]):0}
    else:
        alignmentMappingList=collapseRefs([newClusterDictRef[eachEl]],newClusterDictRef2,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
        alignmentMappingIndex=alignmentMappingList[0]
        if alignmentMappingIndex=='Left':
            newClusterDictRef2.append(newClusterDictRef[eachEl])
            newClusterDictRef2Dict[newClusterDictRef[eachEl]]={(newClusterDictCluster[eachEl],newClusterDictRef[eachEl]):0}
        else:
            newClusterDictRef2Dict[newClusterDictRef2[alignmentMappingIndex]][(newClusterDictCluster[eachEl],newClusterDictRef[eachEl])]=0
            
with open(inputSample+'_finalRefCollapseOutput_v2_20190318.pckl','wb') as f:
    cPickle.dump([newClusterDictRef2Dict,newClusterDictRef2],f,protocol=cPickle.HIGHEST_PROTOCOL)
