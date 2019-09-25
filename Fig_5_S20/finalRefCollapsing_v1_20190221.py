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

sampleName=sys.argv[1]
with open(sampleName+'_plotHairpinFunctionBestIndexVal_v4_20190220_allLists.pckl','rb') as f:
    [xToPlot1,yToPlot1,listAllRefs1,listAllClusters1,listdiGRef1]=cPickle.load(f)

with open(sampleName+'_crossContaminationCheckResultsMatrix_v1_20190220.pckl','rb') as f:
    resultsMatrix=cPickle.load(f)

crossContaminationThreshold=0.8
toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
softClippingThresh=2 #for length upgrading

sampleID=int(sampleName[:sampleName.find('_')])-1

newRefList=[]
newBestIndexValList=[]
newTotalCount=[]
newClusterDictRef=[]
newClusterDictCluster=[]

for counter, eachEl in enumerate(listdiGRef1):
    crossContVal=float(resultsMatrix[sampleID,counter])/float(np.sum(resultsMatrix[:,counter]))
    if crossContVal>crossContaminationThreshold:
        newRefList.append(eachEl)
        newBestIndexValList.append(yToPlot1[counter])
        newTotalCount.append(xToPlot1[counter])
        newClusterDictRef.append(listAllRefs1[counter])
        newClusterDictCluster.append(listAllClusters1[counter])

sortedIndices=np.argsort(np.array(newBestIndexValList))
sortedIndices=sortedIndices[::-1] #decreasing order

newRefList2=[]
newBestIndexValList2=[]
newTotalCount2=[]
newClusterDictRef2=[]
newClusterDictCluster2=[]

for counter, eachEl in enumerate(sortedIndices):
    if counter==0:
        newRefList2.append(newRefList[eachEl])
        newBestIndexValList2.append(newBestIndexValList[eachEl])
        newTotalCount2.append(newTotalCount[eachEl])
        newClusterDictRef2.append(newClusterDictRef[eachEl])
        newClusterDictCluster2.append(newClusterDictCluster[eachEl])
    else:
        alignmentMappingList=collapseRefs([newRefList[eachEl]],newRefList2,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
        alignmentMappingIndex=alignmentMappingList[0]
        if alignmentMappingIndex=='Left':
            newRefList2.append(newRefList[eachEl])
            newBestIndexValList2.append(newBestIndexValList[eachEl])
            newTotalCount2.append(newTotalCount[eachEl])
            newClusterDictRef2.append(newClusterDictRef[eachEl])
            newClusterDictCluster2.append(newClusterDictCluster[eachEl])    

print sampleName+', Number of references!='+str(len(newRefList2))
with open(sampleName+'_finalRefCollapsing_v1_20190221_allLists.pckl','wb') as f:
    cPickle.dump([newTotalCount2,newBestIndexValList2,newClusterDictRef2,newClusterDictCluster2,newRefList2],f,protocol=cPickle.HIGHEST_PROTOCOL)
