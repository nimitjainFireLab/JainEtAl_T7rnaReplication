from moduleUsefulFunctions_20180215 import *
plt.close()

class dimerObject:
    def __init__(self):
        self.firstHalfMutDict={} #keys are various mutations, values are lists [posMutationNotAdjusted,counts]
        self.secondHalfMutDict={} #keys are various mutations, values are lists [posMutationNotAdjusted,counts]
        self.dimerCount=0.0

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

def returnMutationDict(alignInfoList,boundsList,k):
    mutDict={}
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
                    mutDict[('I',eachEl.position-k,eachEl.notes)]=eachEl.position #note this code can be buggy if reads with more than one mutation in one dimer half are allowed as insertion of same 2 bases without collapsing may make same entry in dictionary
            elif eachEl.typeOfEvent=='Mismatch':
                mutDict[('M',eachEl.position-k,eachEl.notes)]=eachEl.position
            elif eachEl.typeOfEvent=='Deletion':
                mutDict[('D',eachEl.position-k)]=eachEl.position           
    return mutDict

def diversityStringReturn(alignmentInfoList,boundsList,k):
    stringDiversity=''
    for eachEl in alignmentInfoList:
        if eachEl.typeOfEvent=='Diversity' and checkWithinBounds(eachEl.position,[boundsList[k]]):
            stringDiversity+=eachEl.notes
    return stringDiversity

def incrementHalfDict(halfDict,mutationDict,counts):
    for eachEl in mutationDict:
        if eachEl in halfDict:
            halfDict[eachEl][1]+=counts
        else:
            halfDict[eachEl]=[mutationDict[eachEl],counts]

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

readsAgreeByPosition=[]
readsDisagreeByPosition=[]
nullHypothesisByPosition=[]
indexDict={}
counterIndex=0
posVector=[]
for eachTup in boundsList:
    for i in range(eachTup[0],eachTup[1]+1):
        readsAgreeByPosition.append(0.0)
        readsDisagreeByPosition.append(0.0)
        nullHypothesisByPosition.append(0.0)
        posVector.append(i)
        indexDict[i]=counterIndex
        counterIndex+=1

dimerPairDict={}
readsAgreeByPosition=np.array(readsAgreeByPosition)
readsDisagreeByPosition=np.array(readsDisagreeByPosition)
nullHypothesisByPosition=np.array(nullHypothesisByPosition)
totalReads=0.0
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
                                dimerPairDict[(firstDiv,secondDiv)]=dimerObject()
                            dimerPairDict[(firstDiv,secondDiv)].dimerCount+=newObj.count[pcrDuplicates]                                
                            totalReads+=newObj.count[pcrDuplicates]
                            mutDictFirstHalf=returnMutationDict(newObj.alignmentInfo[1],boundsFirstHalf,0)
                            mutDictSecondHalf=returnMutationDict(newObj.alignmentInfo[1],boundsSecondHalf,64)                            
                            incrementHalfDict(dimerPairDict[(firstDiv,secondDiv)].firstHalfMutDict,mutDictFirstHalf,newObj.count[pcrDuplicates])
                            incrementHalfDict(dimerPairDict[(firstDiv,secondDiv)].secondHalfMutDict,mutDictSecondHalf,newObj.count[pcrDuplicates])
                            for counter, eachDict in enumerate([mutDictFirstHalf,mutDictSecondHalf]):
                                dictToUse=eachDict
                                if counter==0:                                    
                                    dictToCheck=mutDictSecondHalf
                                elif counter==1:
                                    dictToCheck=mutDictFirstHalf
                                for eachMut in dictToUse:
                                    if eachMut in dictToCheck:
                                        readsAgreeByPosition[indexDict[dictToUse[eachMut]]]+=newObj.count[pcrDuplicates]
                                    else:
                                        readsDisagreeByPosition[indexDict[dictToUse[eachMut]]]+=newObj.count[pcrDuplicates]


for eachComb in dimerPairDict:
    for eachEl in dimerPairDict[eachComb].firstHalfMutDict:
        if eachEl in dimerPairDict[eachComb].secondHalfMutDict:
            firstHalfPos=dimerPairDict[eachComb].firstHalfMutDict[eachEl][0]
            secondHalfPos=dimerPairDict[eachComb].secondHalfMutDict[eachEl][0]
            firstHalfCounts=dimerPairDict[eachComb].firstHalfMutDict[eachEl][1]
            secondHalfCounts=dimerPairDict[eachComb].secondHalfMutDict[eachEl][1]
            expectedCounts=(float(firstHalfCounts)*float(secondHalfCounts))/float(dimerPairDict[eachComb].dimerCount)           
            nullHypothesisByPosition[indexDict[firstHalfPos]]+=expectedCounts
            nullHypothesisByPosition[indexDict[secondHalfPos]]+=expectedCounts        
            
sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.find('_trimmomatic')]
plt.plot(np.arange(1,len(readsAgreeByPosition)+1),readsAgreeByPosition,'r-')
plt.plot(np.arange(1,len(readsDisagreeByPosition)+1),readsDisagreeByPosition,'b-')
plt.plot(np.arange(1,len(nullHypothesisByPosition)+1),nullHypothesisByPosition,'g-')
ax=plt.gca()
ax.set_xlim([0,len(readsAgreeByPosition)+1])
plt.savefig(sampleName+'_wholeDimerMutations_20180417_ver4.pdf',dpi=300)
with open(sampleName+'_wholeDimerMutations_20180417_ver4.pckl','wb') as f:
    cPickle.dump([posVector,readsAgreeByPosition,readsDisagreeByPosition,nullHypothesisByPosition,totalReads],f,protocol=cPickle.HIGHEST_PROTOCOL)

