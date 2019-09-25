from moduleUsefulFunctions_20180215 import *
plt.close()

class dimerSixBaseObject:
    def __init__(self):
        self.readSeqDict={} #keys are tuples of the form ('Mutations in 1st half','Dimer junction','Mutations in 2nd half','Three Prime Base addition')

def returnStartPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            return eachEl.position

def returnEndPosition(alignInfoList):
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='End':
            return eachEl.position

def isFullLength(alignmentObject):
    refSeq=alignmentObject.refUsed
    startPos=returnStartPosition(alignmentObject.alignmentInfo[1])
    endPos=returnEndPosition(alignmentObject.alignmentInfo[1])
    if startPos==0 and endPos==len(refSeq):
        return True
    else:
        return False

def checkWithinBounds(position,boundsList):
    withinBounds=0
    for eachTup in boundsList:
        if position>=eachTup[0] and position<=eachTup[1]:
            withinBounds+=1
    if withinBounds==1:
        return True
    elif withinBounds>1:
        print 'What??!!'
    elif withinBounds==0:
        return False

def returnNumberMutations(alignInfoList,boundsList):
    numberMutations=0
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            startPosition=eachEl.position
        elif eachEl.typeOfEvent=='End':
            endPosition=eachEl.position
    for eachEl in alignInfoList:
        if checkWithinBounds(eachEl.position,boundsList):
            if eachEl.typeOfEvent=='Insertion':
                if eachEl.position==startPosition or eachEl.position==endPosition: #only count insertions in the middle of the sequence as mutations
                    pass
                else:
                    numberMutations+=1 #single base insertions counted as one mutation
            elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
                numberMutations+=1
    return numberMutations

def returnMutationString(alignInfoList,boundsList,k):
    mutString=''
    for eachEl in alignInfoList:
        if eachEl.typeOfEvent=='Start':
            startPosition=eachEl.position
        elif eachEl.typeOfEvent=='End':
            endPosition=eachEl.position
    for eachEl in alignInfoList:
        if checkWithinBounds(eachEl.position,boundsList):
            if eachEl.typeOfEvent=='Insertion':
                if eachEl.position==startPosition or eachEl.position==endPosition: #only count insertions in the middle of the sequence as mutations
                    pass
                else:
                    mutString+='I'+str(eachEl.position-k)+eachEl.notes
            elif eachEl.typeOfEvent=='Mismatch':
                mutString+='M'+str(eachEl.position-k)+eachEl.notes
            elif eachEl.typeOfEvent=='Deletion':
                mutString+='D'+str(eachEl.position-k)
    return mutString

def diversityStringReturn(alignmentInfoList,boundsList,k):
    stringDiversity=''
    for eachEl in alignmentInfoList:
        if eachEl.typeOfEvent=='Diversity' and checkWithinBounds(eachEl.position,[boundsList[k]]):
            stringDiversity+=eachEl.notes
    return stringDiversity

def junctionExtractor(alignmentInfoList):
    startRecording=0
    junctionSeq=''
    for eachEl in alignmentInfoList:
        if startRecording==1:
            junctionSeq+=eachEl.notes
            if eachEl.typeOfEvent=='Match' and eachEl.position==66:
                startRecording=0
                break
        if eachEl.typeOfEvent=='Match' and eachEl.position==62:
            startRecording=1
    return junctionSeq

def threePrimeEndExtractor(alignmentInfoList):
    startRecording=0
    threePrimeEndSeq=''
    for eachEl in alignmentInfoList:
        if startRecording==1:
            threePrimeEndSeq+=eachEl.notes
        if eachEl.typeOfEvent=='Match' and eachEl.position==126:
            startRecording=1
    return threePrimeEndSeq


numberMutationMax=1
boundsList=[(2,61),(66,125)] #made boundsList symmetric in both halves of dimer
boundsFirstHalf=[boundsList[0]]
boundsSecondHalf=[boundsList[1]]
secondAppendPos=64
diversityPosNumber=6
refsToCount=['GG','CC']
pcrDuplicates=0
with open(sys.argv[1],'rb') as f:
    newObjDict=cPickle.load(f)

dimerPairDict={}
#filter for reads that have <=1 mutation in 1st dimer half and <=1 mutation in 2nd dimer half; do not collapse multiple base insertions together
for eachRef in newObjDict:
    letterToAppend=eachRef[0]
    secondLetterToAppend=eachRef[secondAppendPos]
    if letterToAppend+secondLetterToAppend in refsToCount:
        for newObj in newObjDict[eachRef]:
            if isFullLength(newObj):               
                if returnNumberMutations(newObj.alignmentInfo[1],boundsFirstHalf)<=numberMutationMax and returnNumberMutations(newObj.alignmentInfo[1],boundsSecondHalf)<=numberMutationMax:
                    firstDiv=diversityStringReturn(newObj.alignmentInfo[1],boundsList,0)
                    secondDiv=diversityStringReturn(newObj.alignmentInfo[1],boundsList,1)
                    if len(firstDiv)==diversityPosNumber and len(secondDiv)==diversityPosNumber:
                        firstDiv=letterToAppend+firstDiv
                        secondDiv=secondLetterToAppend+secondDiv
                        if firstDiv==secondDiv:
                            if (firstDiv,secondDiv) in dimerPairDict:
                                pass
                            else:
                                dimerPairDict[(firstDiv,secondDiv)]=dimerSixBaseObject()
                            mutStringFirstHalf=returnMutationString(newObj.alignmentInfo[1],boundsFirstHalf,0)
                            mutStringSecondHalf=returnMutationString(newObj.alignmentInfo[1],boundsSecondHalf,64)
                            junctionSeq=junctionExtractor(newObj.alignmentInfo[1])
                            threePrimeSeq=threePrimeEndExtractor(newObj.alignmentInfo[1])
                            tupleForRead=(mutStringFirstHalf,junctionSeq,mutStringSecondHalf,threePrimeSeq)
                            if tupleForRead in dimerPairDict[(firstDiv,secondDiv)].readSeqDict:
                                dimerPairDict[(firstDiv,secondDiv)].readSeqDict[tupleForRead]+=newObj.count[pcrDuplicates]                                
                            else:
                                dimerPairDict[(firstDiv,secondDiv)].readSeqDict[tupleForRead]=newObj.count[pcrDuplicates]                                
sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.find('_trimmomatic')]

with open(sampleName+'_junctions_20180430_ver2.pckl','wb') as f:
    cPickle.dump(dimerPairDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

