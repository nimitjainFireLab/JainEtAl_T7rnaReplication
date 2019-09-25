import sys

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

with open(sys.argv[1],'rb') as f:
    swCollapsedDict=cPickle.load(f)

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.find('_smithWatermanCollapsedClusterDict.pckl')]
with open('../'+sampleName+'_testDict.pckl','rb') as f:
    testDict=cPickle.load(f)

pcrDuplicates=0
seqDict={}
possibleClip=5 #up to how many bases can be clipped at either end

for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][pcrDuplicates]

sortedList=dictSorter(seqDict)
newClusterRefs={}
for eachRef in swCollapsedDict:
    #print eachRef
    if len(eachRef)>2*possibleClip:
        startPositions=np.zeros(possibleClip)
        endPositions=np.zeros(possibleClip)
        for eachSeq in sortedList:
            alignOut=smithWaterman(eachRef,eachSeq)
            startPosition=returnStartPosition(alignOut[1])
            endPosition=returnEndPosition(alignOut[1])
            #print eachSeq
            if internalMutationsNotThere(alignOut[1]) and startPosition<possibleClip and endPosition>len(eachRef)-possibleClip:
                #print (eachSeq,startPosition,endPosition)
                startPositions[startPosition]+=seqDict[eachSeq]
                endPositions[len(eachRef)-endPosition]+=seqDict[eachSeq]
        #print startPositions
        #print endPositions
        bestStart=np.argsort(startPositions)[-1]
        bestEnd=len(eachRef)-np.argsort(endPositions)[-1]
        #print eachRef[bestStart:bestEnd]
        newClusterRefs[eachRef]=(eachRef[bestStart:bestEnd],swCollapsedDict[eachRef])
    else:
        newClusterRefs[eachRef]=(eachRef,swCollapsedDict[eachRef])

print newClusterRefs

with open(sampleName+'_collapsedEndsClusterDict.pckl','wb') as f:
    cPickle.dump(newClusterRefs,f,protocol=cPickle.HIGHEST_PROTOCOL)

