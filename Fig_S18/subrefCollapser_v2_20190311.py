from random import randint
import sys
import resource
from moduleUsefulFunctions_20180215 import *

#ver16 denovoClustering code onwards definition of clusterObject
class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

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

with open(sampleName+'_clusterDict_v16.pckl','rb') as f:
    clusterDict=cPickle.load(f)

newClusterDict={} #would not contain the Clashed and Unaligned cluster, also will remove Unaligned subcluster for all clusters

toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
softClippingThresh=2 #for length upgrading

for counter, eachCluster in enumerate(clusterDict):
    if eachCluster!='Clashed' and eachCluster!='Unaligned':
        dictForAlignment1={}
        dictForAlignment0={}
        for eachRef in clusterDict[eachCluster]:
            if eachRef!='Unaligned':
                if clusterDict[eachCluster][eachRef].forAlignment==0:
                    if len(clusterDict[eachCluster][eachRef].listSeqIndices)>1:
                        print 'What??!'
                    dictForAlignment0[eachRef]=clusterDict[eachCluster][eachRef].numberOfReads
                elif clusterDict[eachCluster][eachRef].forAlignment==1:
                    dictForAlignment1[eachRef]=clusterDict[eachCluster][eachRef].numberOfReads
        sortedList1=dictSorter(dictForAlignment1)
        sortedList0=dictSorter(dictForAlignment0)
        if counter%1000==0:
            print 'Working on cluster '+str(counter+1)+' of '+str(len(clusterDict))+' clusters'
        #print 'Step 1: Align subrefs with forAlignment=0 to subrefs with forAlignment=1'
        subrefs0Left=[]
        if len(sortedList0)>0:
            alignmentMappingList=collapseRefs(sortedList0,sortedList1,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
            for counter1, eachEl in enumerate(alignmentMappingList):
                if eachEl=='Left':
                    subrefs0Left.append(sortedList0[counter1])
                else:
                    for eachKmer in clusterDict[eachCluster][sortedList0[counter1]].dictKmers:
                        clusterDict[eachCluster][sortedList1[eachEl]].dictKmers[eachKmer]=0
                    clusterDict[eachCluster][sortedList1[eachEl]].numberOfReads+=clusterDict[eachCluster][sortedList0[counter1]].numberOfReads
                    clusterDict[eachCluster][sortedList1[eachEl]].listSeqIndices=list(clusterDict[eachCluster][sortedList1[eachEl]].listSeqIndices+clusterDict[eachCluster][sortedList0[counter1]].listSeqIndices)
            for eachEl in sortedList1:
                clusterDict[eachCluster][eachEl].listSeqIndices.sort()
        #print 'Step 1 done'
        
        #print 'Step 2: Collapsing together subrefs with forAlignment=1'
        #new sorted order for current forAlignment=1 subreferences based on numberOfReads may be different than the order of sortedList1, but still use the old sorted order
        list1=[]
        for counter1, eachEl in enumerate(sortedList1):
            if counter1==0:
                list1.append(eachEl)
            else:
                alignmentMappingList=collapseRefs([eachEl],list1,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
                alignmentMappingIndex=alignmentMappingList[0]
                if alignmentMappingIndex=='Left':
                    list1.append(eachEl)
                else:
                    for eachKmer in clusterDict[eachCluster][eachEl].dictKmers:
                        clusterDict[eachCluster][list1[alignmentMappingIndex]].dictKmers[eachKmer]=0
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].numberOfReads+=clusterDict[eachCluster][eachEl].numberOfReads
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices=list(clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices+clusterDict[eachCluster][eachEl].listSeqIndices)
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices.sort()
                    del clusterDict[eachCluster][eachEl]
        #print 'Step 2 done'

        #print 'Step 3: Collapsing together subrefs with forAlignment=0'
        list1=[]
        for counter1, eachEl in enumerate(subrefs0Left): #note that subrefs0Left should already be in sorted order of counts based on how it is made
            if counter1==0:
                list1.append(eachEl)
                clusterDict[eachCluster][eachEl].forAlignment=1
            else:
                alignmentMappingList=collapseRefs([eachEl],list1,toleranceNumberMutations,ntLengthUpdateFullLength,softClippingThresh)
                alignmentMappingIndex=alignmentMappingList[0]
                if alignmentMappingIndex=='Left':
                    list1.append(eachEl)
                    clusterDict[eachCluster][eachEl].forAlignment=1
                else:
                    for eachKmer in clusterDict[eachCluster][eachEl].dictKmers:
                        clusterDict[eachCluster][list1[alignmentMappingIndex]].dictKmers[eachKmer]=0
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].numberOfReads+=clusterDict[eachCluster][eachEl].numberOfReads
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices=list(clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices+clusterDict[eachCluster][eachEl].listSeqIndices)
                    clusterDict[eachCluster][list1[alignmentMappingIndex]].listSeqIndices.sort()
        #print 'Step 3 done'

        #print 'Step 4: Writing new cluster'
        newClusterDict[eachCluster]={}
        for eachRef in clusterDict[eachCluster]:
            if clusterDict[eachCluster][eachRef].forAlignment==1:
                newClusterDict[eachCluster][eachRef]=clusterObject(clusterDict[eachCluster][eachRef].reference,dict(clusterDict[eachCluster][eachRef].dictKmers),clusterDict[eachCluster][eachRef].numberOfReads,clusterDict[eachCluster][eachRef].forAlignment,list(clusterDict[eachCluster][eachRef].listSeqIndices))
        #print 'Step 4 done'
            
with open(sampleName+'_clusterDict_subrefCollapsed_v1.pckl','wb') as f:
    cPickle.dump(newClusterDict,f,protocol=cPickle.HIGHEST_PROTOCOL)
