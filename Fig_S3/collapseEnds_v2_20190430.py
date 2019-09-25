from moduleUsefulFunctions_20180215 import *

def internalMutationsNotThere(alignInfoList):
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
                numberMutations+=1 #for the purposes here, it doesn't matter how many bases inserted as action item is based on whether there was any insertion at all or not
        elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
            numberMutations+=1
    if numberMutations>0:
        return False
    else:
        return True

def returnStartPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            return eachEl.position

def returnEndPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='End':
            return eachEl.position


sortWithPCRDuplicates=0
k=20 #used to define identity of a reference sequence

def returnTestDict(sampleName):
    with open(sampleName+'_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_testDict.pckl','rb') as f:
        testDict=cPickle.load(f)
    seqDict={}
    for eachEntry in testDict:
        seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]
    sortedDictionary=dictSorter(seqDict)
    return [testDict,sortedDictionary,seqDict]

inputSample=sys.argv[1]
[inputTest,inputSorted,inputSeqDict]=returnTestDict(inputSample)

totalCounts=0
for eachEl in inputSeqDict:
    if len(eachEl)>=k: #corresponds to length condition check in denovoClustering_ver16_20190111.py
        totalCounts+=inputSeqDict[eachEl]

with open(inputSample+'_assembleCollapsedRefs_v3_20190428.pckl','rb') as f:
    [listRefs,countRefs,dictSequencesRefs]=cPickle.load(f)

dictRefs={}
for counter, eachEl in enumerate(listRefs):
    dictRefs[eachEl]=(float(countRefs[counter])/float(totalCounts))*100.0

possibleClip=5 #up to how many bases can be clipped at either end

newClusterRefs={}
for counter, eachRef in enumerate(dictSorter(dictRefs)):
    print eachRef+'\t'+str(dictRefs[eachRef])
    if counter==0:
        startPositions=np.zeros(possibleClip)
        endPositions=np.zeros(possibleClip)
        for eachSeq in inputSorted:
            alignOut=smithWaterman(eachRef,eachSeq)
            startPosition=returnStartPosition(alignOut[1])
            endPosition=returnEndPosition(alignOut[1])
            if internalMutationsNotThere(alignOut[1]) and startPosition<possibleClip and endPosition>len(eachRef)-possibleClip:
                startPositions[startPosition]+=inputSeqDict[eachSeq]
                endPositions[len(eachRef)-endPosition]+=inputSeqDict[eachSeq]
        bestStart=np.argsort(startPositions)[-1]
        bestEnd=len(eachRef)-np.argsort(endPositions)[-1]
        newClusterRefs[eachRef]=(eachRef[bestStart:bestEnd],dictRefs[eachRef])

print newClusterRefs

with open(inputSample+'_collapseEnds_v2_20190430.pckl','wb') as f:
    cPickle.dump(newClusterRefs,f,protocol=cPickle.HIGHEST_PROTOCOL)

    
        



