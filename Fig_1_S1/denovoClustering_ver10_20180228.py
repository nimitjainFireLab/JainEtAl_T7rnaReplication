from random import randint
import sys
import itertools

from moduleUsefulFunctions_20180215 import *

###########################

#FUNCTION SECTION

###########################

class clusterObject:
    def __init__(self,reference,dictKmers,numberOfReads):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads

def printClusterDict(dictionaryToPrint):
    for counter, eachCluster in enumerate(dictionaryToPrint):
        print 'Cluster '+str(counter+1)+': '+eachCluster
        for eachRef in dictionaryToPrint[eachCluster]:
            print '\t'+eachRef+'\t'+str(dictionaryToPrint[eachCluster][eachRef].numberOfReads)

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

def upgradeMappingDict(kmerDictToAdd,mappingDict,clusterRefSeq):
    for eachKmer in kmerDictToAdd:
        mappingDict[eachKmer]=clusterRefSeq

def returnAllTuplesClashes(dictionaryToUse):
    listToReturn=[]
    allPerms=itertools.permutations(dictionaryToUse.keys())
    for eachEl in allPerms:
        listToReturn.append(eachEl)
    return listToReturn


###########################

#END OF FUNCTION SECTION

###########################

###########################

#PARAMETER SECTION

###########################

sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same
k=20 #used to define identity of a reference sequence
toleranceNumberMutations=1.0/float(k) #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
fractionCountsUpdate=0.05 #if counts for a particular sequence are lower than this fraction compared to a particular cluster or subcluster, then don't change things based on this sequence
readThresholdClustering=2 #if counts less than or equal to this value, then don't start new subclusters or clusters--instead assign things to an unaligned bin
softClippingThresh=2 #for length upgrading
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
readCutoffUpgrades=9 #for length updates or cluster merging
###########################

#END OF PARAMETER SECTION

###########################


###########################

#DATA SECTION

###########################

#Sort testDict results for clustering
with open(sys.argv[1],'rb') as f:
    testDict=cPickle.load(f)

sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_testDict.pckl')]

clusterDictName=sampleName+'_clusterDict.pckl'
clusterMetricsName=sampleName+'_clusterMetrics.pckl'
unalignedClusterDictName=sampleName+'_unalignedCluster.pckl'
readsShorterThanKDictName=sampleName+'_clusterReadsShorterThanK.pckl'
clashedClusterDictName=sampleName+'_clashedCluster.pckl'
seqDict={}
readsShorterThanK={}
clusterDict={} #keys are "reference" sequences, these go into subclusters whose keys are subreference sequences and values include dictionary of kmers, number of reads that map to the reference (and potentially other things in the future)
clusterMappingDict={}  #keys are kmers and values are cluster reference sequences
unalignedSequencesDict={}
clashedSequencesDict={}

for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)
print len(seqDict)
numberSeqProcessed=0
for eachEl in sortedDictionary:
    numberSeqProcessed+=1
    insertKmers={}
    if len(eachEl)>=k:
        for i in range(0,len(eachEl)+1-k):
            insertKmers[eachEl[i:i+k]]=0
            insertKmers[revComp(eachEl[i:i+k])]=0
        clustersFound={} #keys are cluster "reference" sequences
        for eachKmer in insertKmers:
            if eachKmer in clusterMappingDict:
                clustersFound[clusterMappingDict[eachKmer]]=0
        if len(clustersFound)!=0:
            if len(clustersFound)==1: #only one matching reference found, this is the candidate reference, now align with all subclusters and if alignment is good update with kmers from eachEl, otherwise start new subcluster--subclusters are to be able to distinguish sequences with long insertions and deletions with respect to each other
                keysList=clustersFound.keys()
                desiredKey=keysList[0]
                alignmentInfo=[]
                alignmentScore=[]
                for counter, eachRef in enumerate(clusterDict[desiredKey]):
                    if eachRef!='Unaligned':
                        score=-1
                        listRefs=[eachRef,revComp(eachRef)]
                        for eachSeq in listRefs:
                            alignOut=smithWaterman(eachSeq,eachEl,collapseInsertions=False)
                            if alignOut[0]>score:
                                score=alignOut[0]
                                info=alignedObject(alignOut,eachRef,-1) #note eachRef is specified as refUsed, not eachSeq; this is useful in downstream code below; artificial count of -1 passed
                        alignmentInfo.append(info)
                        alignmentScore.append(float(score)/min([len(eachEl),len(eachRef)])) #normalization by length is useful in downstream code below
                    else:
                        alignmentInfo.append('')
                        alignmentScore.append(0)
                alignmentScore=np.array(alignmentScore)
                indicesSort=np.argsort(alignmentScore)
                subclusterIndex=indicesSort[-1] #note there may be clashes between two subcluster references to receive the counts of the current eachEl. Here, by default, the best subcluster reference gets the counts even if the next best subcluster reference has the same numberMutations in reality (whatever sorts to the -1 position)
                numberMutations=returnNumberMutations(alignmentInfo[subclusterIndex].alignmentInfo[1])
                startNewSubcluster=0
                justAssign=0
                lengthUpdate=0
                endPositionBestAlignment=returnEndPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
                startPositionBestAlignment=returnStartPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
                maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(alignmentInfo[subclusterIndex].refUsed),len(eachEl)]))
                totalNumberMutationsCalc=len(FivePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+len(ThreePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+numberMutations
                if numberMutations<=np.ceil(toleranceNumberMutations*(endPositionBestAlignment-startPositionBestAlignment)):
                    if startPositionBestAlignment<softClippingThresh and endPositionBestAlignment>(len(alignmentInfo[subclusterIndex].refUsed)-softClippingThresh):
                        if (len(eachEl)>(len(alignmentInfo[subclusterIndex].refUsed)+ntLengthUpdateFullLength)) and (seqDict[eachEl]>=fractionCountsUpdate*clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].numberOfReads) and seqDict[eachEl]>readCutoffUpgrades: #update reference length only if length criterion satisfied and counts are above threshold                            
                            lengthUpdate=1
                        justAssign=1                        
                    else:
                        if totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength):
                            justAssign=1
                        else:
                            startNewSubcluster=1
                    if justAssign==1:
                        #update clusterObject
                        for eachKmer in insertKmers:
                            clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].dictKmers[eachKmer]=0
                        clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].numberOfReads+=seqDict[eachEl]
                        upgradeMappingDict(insertKmers,clusterMappingDict,desiredKey)                        
                    if lengthUpdate==1:                    
                        clusterDict[desiredKey][eachEl]=clusterObject(eachEl,dict(clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].dictKmers),clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].numberOfReads)
                        del clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed]                    
                    ###
                    '''
                    There is a drawback of incremeneting length on the fly like this. Sequences that came previously may have aligned differently to the previous reference as compared to the new reference.
                    In addition, say a sequence which served as a reference was actually common to 2 longer references. Then, only one of the longer references will "inherit" counts from the shorter reference, when in reality it is unknown whether the counts for the shorter reference correspond to longer reference 1 or to longer reference 2.
                    However, I still decided to implement increase in length on the fly because it does increase accuracy a little bit (refer to the case of the 2 longer references above; if length was not increased, then both references will be merged into the shorter reference). 
                    '''
                    ###
                else:
                    startNewSubcluster=1
                if startNewSubcluster==1:
                    # new subcluster
                    if seqDict[eachEl]>readThresholdClustering:
                        clusterDict[desiredKey][eachEl]=clusterObject(eachEl,insertKmers,seqDict[eachEl])
                        upgradeMappingDict(insertKmers,clusterMappingDict,desiredKey)                        
                    else:
                        if 'Unaligned' in clusterDict[desiredKey]:
                            clusterDict[desiredKey]['Unaligned'].numberOfReads+=seqDict[eachEl]
                        else:
                            clusterDict[desiredKey]['Unaligned']=clusterObject('Unaligned',{},seqDict[eachEl])                        
            else: #multiple matching references found, merge these together                    
                maxNumberClusterReads=0
                for eachClusterSeq in clustersFound:                        
                    totalNumberReadsAssociatedWithCluster=0
                    for eachRef in clusterDict[eachClusterSeq]:
                        totalNumberReadsAssociatedWithCluster+=clusterDict[eachClusterSeq][eachRef].numberOfReads
                    if totalNumberReadsAssociatedWithCluster>maxNumberClusterReads:
                        maxNumberClusterReads=totalNumberReadsAssociatedWithCluster
                if seqDict[eachEl]>=fractionCountsUpdate*maxNumberClusterReads and seqDict[eachEl]>readCutoffUpgrades:
                #new cluster for the bridging sequence eachEl
                    clusterDict[eachEl]={eachEl:clusterObject(eachEl,insertKmers,seqDict[eachEl])}
                    upgradeMappingDict(insertKmers,clusterMappingDict,eachEl)                        
                    for eachClusterSeq in clustersFound:                            
                        for counter, eachRef in enumerate(clusterDict[eachClusterSeq]):
                            if eachRef!='Unaligned':
                                score=-1
                                listRefs=[eachRef,revComp(eachRef)]
                                for eachSeq in listRefs:
                                    alignOut=smithWaterman(eachEl,eachSeq,collapseInsertions=False) #note how order of eachEl and eachSeq is flipped as compared to previous invocation of smithWaterman above
                                    if alignOut[0]>score:
                                        score=alignOut[0]
                                        info=alignOut[1]
                                numberMutations=returnNumberMutations(info)
                                maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(eachRef),len(eachEl)]))
                                totalNumberMutationsCalc=len(FivePrimeClipSequence(info))+len(ThreePrimeClipSequence(info))+numberMutations
                                startNewSubcluster=0
                                if numberMutations<=np.ceil(toleranceNumberMutations*(returnEndPosition(info)-returnStartPosition(info))):
                                    #update the eachEl subcluster in eachEl cluster
                                    if totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength):
                                        for eachKmer in clusterDict[eachClusterSeq][eachRef].dictKmers:
                                            clusterDict[eachEl][eachEl].dictKmers[eachKmer]=0
                                        clusterDict[eachEl][eachEl].numberOfReads+=clusterDict[eachClusterSeq][eachRef].numberOfReads
                                    else:
                                        startNewSubcluster=1
                                else: #add a new subcluster to the eachEl cluster
                                    startNewSubcluster=1
                                if startNewSubcluster==1:
                                    clusterDict[eachEl][eachRef]=clusterObject(eachRef,dict(clusterDict[eachClusterSeq][eachRef].dictKmers),clusterDict[eachClusterSeq][eachRef].numberOfReads)
                                upgradeMappingDict(clusterDict[eachClusterSeq][eachRef].dictKmers,clusterMappingDict,eachEl)                        
                            else:
                                if 'Unaligned' in clusterDict[eachEl]:
                                    clusterDict[eachEl]['Unaligned'].numberOfReads+=clusterDict[eachClusterSeq]['Unaligned'].numberOfReads
                                else:
                                    clusterDict[eachEl]['Unaligned']=clusterObject('Unaligned',{},clusterDict[eachClusterSeq]['Unaligned'].numberOfReads)
                        del clusterDict[eachClusterSeq]
                else: #Put in a clash cluster
                    clustersThatMatchList=returnAllTuplesClashes(clustersFound) #list of tuples
                    matchingCombClusterFound=0
                    if 'Clashed' in clusterDict:
                        for eachPossibleComb in clustersThatMatchList:
                            if str(eachPossibleComb) in clusterDict['Clashed']:
                                matchingCombClusterFound=1
                                clusterDict['Clashed'][str(eachPossibleComb)].numberOfReads+=seqDict[eachEl]
                                break
                        if matchingCombClusterFound==0:
                            clusterDict['Clashed'][str(clustersThatMatchList[0])]=clusterObject(str(clustersThatMatchList[0]),{},seqDict[eachEl])                  
                    else:
                        clusterDict['Clashed']={str(clustersThatMatchList[0]):clusterObject(str(clustersThatMatchList[0]),{},seqDict[eachEl])}
                    clashedSequencesDict[eachEl]=seqDict[eachEl]
        else:
            if seqDict[eachEl]>readThresholdClustering: #Start a new cluster
                clusterDict[eachEl]={eachEl:clusterObject(eachEl,insertKmers,seqDict[eachEl])}
                upgradeMappingDict(insertKmers,clusterMappingDict,eachEl)                    
            else:
                if 'Unaligned' in clusterDict:
                    clusterDict['Unaligned']['Unaligned'].numberOfReads+=seqDict[eachEl]
                else:
                    clusterDict['Unaligned']={'Unaligned':clusterObject('Unaligned',{},seqDict[eachEl])}
                unalignedSequencesDict[eachEl]=seqDict[eachEl]
    else:
        readsShorterThanK[eachEl]=seqDict[eachEl]
    
##    print eachEl+'\t'+str(numberSeqProcessed)
##    printClusterDict(clusterDict)
##    print ''
##    print ''
##    print ''

clusterMetricsDict={}
clusterMetricsDict['codeUsed']=sys.argv[0]
clusterMetricsDict['testDictUsed']=sys.argv[1]
clusterMetricsDict['sortWithPCRDuplicates']=sortWithPCRDuplicates
clusterMetricsDict['k']=k
clusterMetricsDict['toleranceNumberMutations']=toleranceNumberMutations
clusterMetricsDict['softClippingThresh']=softClippingThresh
clusterMetricsDict['fractionCountsUpdate']=fractionCountsUpdate
clusterMetricsDict['readThresholdClustering']=readThresholdClustering
clusterMetricsDict['numberReadsShorterThanK']=len(readsShorterThanK)
clusterMetricsDict['len(clusterDict)']=len(clusterDict)
clusterMetricsDict['numberSeqProcessed']=numberSeqProcessed
clusterMetricsDict['len(unalignedSequencesDict)']=len(unalignedSequencesDict)
clusterMetricsDict['ntLengthUpdateFullLength']=ntLengthUpdateFullLength
clusterMetricsDict['readCutoffUpgrades']=readCutoffUpgrades

with open(clusterMetricsName,'wb') as f:
    cPickle.dump(clusterMetricsDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(clusterDictName,'wb') as f:
    cPickle.dump(clusterDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(unalignedClusterDictName,'wb') as f:
    cPickle.dump(unalignedSequencesDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(clashedClusterDictName,'wb') as f:
    cPickle.dump(clashedSequencesDict,f,protocol=cPickle.HIGHEST_PROTOCOL)
    
with open(readsShorterThanKDictName,'wb') as f:
    cPickle.dump(readsShorterThanK,f,protocol=cPickle.HIGHEST_PROTOCOL)
    
###########################

#END OF DATA SECTION

###########################
