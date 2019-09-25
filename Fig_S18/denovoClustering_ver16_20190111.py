from random import randint
import sys
import resource

from moduleUsefulFunctions_20180215 import *

###########################

#FUNCTION SECTION

###########################

class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

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
readThresholdClustering=0 #if counts less than or equal to this value, then don't start new subclusters or clusters--instead assign things to an unaligned bin
softClippingThresh=2 #for length upgrading
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
readCutoffUpgrades=9 #for length updates or cluster merging
numberSubclusterAlignMax=5 #maximum number of "best" subclusters to use for alignment of a sequence within a cluster--parameter must be >=1
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

clusterDictName=sampleName+'_clusterDict_v16.pckl'
clusterMetricsName=sampleName+'_clusterMetrics_v16.pckl'
unalignedClusterDictName=sampleName+'_unalignedCluster_v16.pckl'
readsShorterThanKDictName=sampleName+'_clusterReadsShorterThanK_v16.pckl'
clashedClusterDictName=sampleName+'_clashedCluster_v16.pckl'
numberAlignRefPerClusterName=sampleName+'_numberAlignRefPerCluster_v16.pckl'
seqDict={}
readsShorterThanK={}
clusterDict={} #keys are "reference" sequences, these go into subclusters whose keys are subreference sequences and values include dictionary of kmers, number of reads that map to the reference (and potentially other things in the future)
clusterMappingDict={}  #keys are kmers and values are cluster reference sequences
unalignedSequencesDict={}
clashedSequencesDict={}
numberAlignRefPerCluster={}

for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)
print len(seqDict)
numberSeqProcessed=0
for eachEl in sortedDictionary:
    numberSeqProcessed+=1
    if numberSeqProcessed%10000==0:
        print 'Processed '+str(numberSeqProcessed)+' of '+str(len(seqDict))
        print 'Memory used: '+str(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        with open(clusterDictName,'wb') as f:
            cPickle.dump(clusterDict,f,protocol=cPickle.HIGHEST_PROTOCOL)
        print ''
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
            if len(clustersFound)==1: #only one matching reference found, this is the candidate reference, now align with up to numberSubclusterAlignMax subclusters and if alignment is good update with kmers from eachEl, otherwise start new subcluster--subclusters are to be able to distinguish sequences with long insertions and deletions with respect to each other
                keysList=clustersFound.keys()
                desiredKey=keysList[0]
                alignmentInfo=[]
                alignmentScore=[]
                for counter, eachRef in enumerate(clusterDict[desiredKey]):
                    if eachRef!='Unaligned' and clusterDict[desiredKey][eachRef].forAlignment==1:
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
                        clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].listSeqIndices.append(numberSeqProcessed-1)
                        upgradeMappingDict(insertKmers,clusterMappingDict,desiredKey)                        
                    if lengthUpdate==1:                    
                        clusterDict[desiredKey][eachEl]=clusterObject(eachEl,dict(clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].dictKmers),clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].numberOfReads,1,list(clusterDict[desiredKey][alignmentInfo[subclusterIndex].refUsed].listSeqIndices)) #as updating length, means alignment must have happened--hence forAlignment=1
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
                        if desiredKey in numberAlignRefPerCluster:
                            if numberAlignRefPerCluster[desiredKey]>=numberSubclusterAlignMax:
                                clusterDict[desiredKey][eachEl]=clusterObject(eachEl,insertKmers,seqDict[eachEl],0,[numberSeqProcessed-1])
                            else:
                                clusterDict[desiredKey][eachEl]=clusterObject(eachEl,insertKmers,seqDict[eachEl],1,[numberSeqProcessed-1])
                                numberAlignRefPerCluster[desiredKey]+=1
                        else:
                            clusterDict[desiredKey][eachEl]=clusterObject(eachEl,insertKmers,seqDict[eachEl],1,[numberSeqProcessed-1])
                            numberAlignRefPerCluster[desiredKey]=1
                        upgradeMappingDict(insertKmers,clusterMappingDict,desiredKey)                        
                    else:
                        if 'Unaligned' in clusterDict[desiredKey]:
                            clusterDict[desiredKey]['Unaligned'].numberOfReads+=seqDict[eachEl]
                        else:
                            clusterDict[desiredKey]['Unaligned']=clusterObject('Unaligned',{},seqDict[eachEl],0,[])                        
            else: #multiple matching references found, merge these together                    
                maxNumberClusterReads=0
                for eachClusterSeq in clustersFound:                        
                    totalNumberReadsAssociatedWithCluster=0
                    for eachRef in clusterDict[eachClusterSeq]:
                        totalNumberReadsAssociatedWithCluster+=clusterDict[eachClusterSeq][eachRef].numberOfReads #includes 'Unaligned' subcluster
                    if totalNumberReadsAssociatedWithCluster>maxNumberClusterReads:
                        maxNumberClusterReads=totalNumberReadsAssociatedWithCluster
                if seqDict[eachEl]>=fractionCountsUpdate*maxNumberClusterReads and seqDict[eachEl]>readCutoffUpgrades:
                #new cluster for the bridging sequence eachEl
                    clusterDict[eachEl]={eachEl:clusterObject(eachEl,insertKmers,seqDict[eachEl],1,[numberSeqProcessed-1])} #its possible that after merging, number of subclusters used for alignment is greater than numberSubclusterAlignMax
                    numberAlignRefPerCluster[eachEl]=1
                    upgradeMappingDict(insertKmers,clusterMappingDict,eachEl)                        
                    for eachClusterSeq in clustersFound:                            
                        for counter, eachRef in enumerate(clusterDict[eachClusterSeq]):
                            if eachRef!='Unaligned':
                                if clusterDict[eachClusterSeq][eachRef].forAlignment==1:
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
                                        if (totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength)) or (returnStartPosition(info)<softClippingThresh and returnEndPosition(info)>(len(eachEl)-softClippingThresh)):
                                            for eachKmer in clusterDict[eachClusterSeq][eachRef].dictKmers:
                                                clusterDict[eachEl][eachEl].dictKmers[eachKmer]=0
                                            clusterDict[eachEl][eachEl].numberOfReads+=clusterDict[eachClusterSeq][eachRef].numberOfReads
                                            clusterDict[eachEl][eachEl].listSeqIndices=list(clusterDict[eachEl][eachEl].listSeqIndices+clusterDict[eachClusterSeq][eachRef].listSeqIndices)
                                            clusterDict[eachEl][eachEl].listSeqIndices.sort()
                                        else:
                                            startNewSubcluster=1
                                    else: #add a new subcluster to the eachEl cluster
                                        startNewSubcluster=1
                                    if startNewSubcluster==1:
                                        clusterDict[eachEl][eachRef]=clusterObject(eachRef,dict(clusterDict[eachClusterSeq][eachRef].dictKmers),clusterDict[eachClusterSeq][eachRef].numberOfReads,1,list(clusterDict[eachClusterSeq][eachRef].listSeqIndices)) #forAlignment was already 1 here--keep it as is
                                        numberAlignRefPerCluster[eachEl]+=1
                                else:
                                    clusterDict[eachEl][eachRef]=clusterObject(eachRef,dict(clusterDict[eachClusterSeq][eachRef].dictKmers),clusterDict[eachClusterSeq][eachRef].numberOfReads,0,list(clusterDict[eachClusterSeq][eachRef].listSeqIndices))
                                upgradeMappingDict(clusterDict[eachClusterSeq][eachRef].dictKmers,clusterMappingDict,eachEl)                        
                            else:
                                if 'Unaligned' in clusterDict[eachEl]:
                                    clusterDict[eachEl]['Unaligned'].numberOfReads+=clusterDict[eachClusterSeq]['Unaligned'].numberOfReads
                                else:
                                    clusterDict[eachEl]['Unaligned']=clusterObject('Unaligned',{},clusterDict[eachClusterSeq]['Unaligned'].numberOfReads,0,[])
                        del clusterDict[eachClusterSeq]
                        del numberAlignRefPerCluster[eachClusterSeq]
                    #check if number of subclusters used for alignment is less than numberSubclusterAlignMax but still there are distinct references with forAlignment=0
                    interimDictRefAlignment0={}
                    interimListRefAlignment1=[]
                    for eachRef in clusterDict[eachEl]:
                        if eachRef!='Unaligned':
                            if clusterDict[eachEl][eachRef].forAlignment==0:
                                interimDictRefAlignment0[eachRef]=clusterDict[eachEl][eachRef].numberOfReads
                            elif clusterDict[eachEl][eachRef].forAlignment==1:
                                interimListRefAlignment1.append(eachRef)
                    interimListRefAlignment0=dictSorter(interimDictRefAlignment0)
                    for eachRef in interimListRefAlignment0:
                        if len(interimListRefAlignment1)>=numberSubclusterAlignMax:
                            break
                        alignmentInfo=[]
                        alignmentScore=[]
                        for eachRef1 in interimListRefAlignment1:
                            score=-1
                            listRefs=[eachRef1,revComp(eachRef1)]
                            for eachSeq in listRefs:
                                alignOut=smithWaterman(eachSeq,eachRef,collapseInsertions=False)
                                if alignOut[0]>score:
                                    score=alignOut[0]
                                    info=alignedObject(alignOut,eachRef1,-1) #note eachRef1 is specified as refUsed
                            alignmentInfo.append(info)
                            alignmentScore.append(float(score)/min([len(eachRef),len(eachRef1)])) #normalization by length is useful in downstream code below
                        alignmentScore=np.array(alignmentScore)
                        indicesSort=np.argsort(alignmentScore)
                        subclusterIndex=indicesSort[-1] 
                        numberMutations=returnNumberMutations(alignmentInfo[subclusterIndex].alignmentInfo[1])
                        endPositionBestAlignment=returnEndPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
                        startPositionBestAlignment=returnStartPosition(alignmentInfo[subclusterIndex].alignmentInfo[1])
                        maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(alignmentInfo[subclusterIndex].refUsed),len(eachRef)]))
                        totalNumberMutationsCalc=len(FivePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+len(ThreePrimeClipSequence(alignmentInfo[subclusterIndex].alignmentInfo[1]))+numberMutations
                        if numberMutations<=np.ceil(toleranceNumberMutations*(endPositionBestAlignment-startPositionBestAlignment)):                        
                            if (totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength)) or (startPositionBestAlignment<softClippingThresh and endPositionBestAlignment>(len(alignmentInfo[subclusterIndex].refUsed)-softClippingThresh)):
                                for eachKmer in clusterDict[eachEl][eachRef].dictKmers:
                                    clusterDict[eachEl][alignmentInfo[subclusterIndex].refUsed].dictKmers[eachKmer]=0
                                clusterDict[eachEl][alignmentInfo[subclusterIndex].refUsed].numberOfReads+=clusterDict[eachEl][eachRef].numberOfReads
                                clusterDict[eachEl][alignmentInfo[subclusterIndex].refUsed].listSeqIndices=list(clusterDict[eachEl][alignmentInfo[subclusterIndex].refUsed].listSeqIndices+clusterDict[eachEl][eachRef].listSeqIndices)
                                clusterDict[eachEl][alignmentInfo[subclusterIndex].refUsed].listSeqIndices.sort()
                                del clusterDict[eachEl][eachRef]
                            else:
                                clusterDict[eachEl][eachRef].forAlignment=1
                                interimListRefAlignment1.append(eachRef)
                                numberAlignRefPerCluster[eachEl]+=1
                        else:
                            clusterDict[eachEl][eachRef].forAlignment=1
                            interimListRefAlignment1.append(eachRef)
                            numberAlignRefPerCluster[eachEl]+=1
                else: #Put in a clash cluster
                    clustersThatMatchList=clustersFound.keys()
                    clustersThatMatchList.sort()
                    clustersThatMatchList=str(tuple(clustersThatMatchList))
                    if 'Clashed' in clusterDict:
                        if clustersThatMatchList in clusterDict['Clashed']:
                            clusterDict['Clashed'][clustersThatMatchList].numberOfReads+=seqDict[eachEl]
                        else:
                            clusterDict['Clashed'][clustersThatMatchList]=clusterObject(clustersThatMatchList,{},seqDict[eachEl],0,[])
                    else:
                        clusterDict['Clashed']={clustersThatMatchList:clusterObject(clustersThatMatchList,{},seqDict[eachEl],0,[])}                    
                    clashedSequencesDict[eachEl]=seqDict[eachEl]
        else:
            if seqDict[eachEl]>readThresholdClustering: #Start a new cluster
                clusterDict[eachEl]={eachEl:clusterObject(eachEl,insertKmers,seqDict[eachEl],1,[numberSeqProcessed-1])} # as first subcluster, forAlignment=1
                upgradeMappingDict(insertKmers,clusterMappingDict,eachEl)                    
                numberAlignRefPerCluster[eachEl]=1
            else:
                if 'Unaligned' in clusterDict:
                    clusterDict['Unaligned']['Unaligned'].numberOfReads+=seqDict[eachEl]
                else:
                    clusterDict['Unaligned']={'Unaligned':clusterObject('Unaligned',{},seqDict[eachEl],0,[])}
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
clusterMetricsDict['numberSubclusterAlignMax']=numberSubclusterAlignMax

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
    
with open(numberAlignRefPerClusterName,'wb') as f:
    cPickle.dump(numberAlignRefPerCluster,f,protocol=cPickle.HIGHEST_PROTOCOL)


###########################

#END OF DATA SECTION

###########################
