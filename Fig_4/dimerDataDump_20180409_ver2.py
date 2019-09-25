from moduleUsefulFunctions_20180215 import *

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
                    numberMutations+=1 #as insertions are not being collapsed together (collapseInsertions=False in the smithWaterman calls below), each base inserted is counted separately
            elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
                numberMutations+=1
    return numberMutations

def diversityStringReturn(alignmentInfoList,boundsList,k):
    stringDiversity=''
    for eachEl in alignmentInfoList:
        if eachEl.typeOfEvent=='Diversity' and checkWithinBounds(eachEl.position,[boundsList[k]]):
            stringDiversity+=eachEl.notes
    return stringDiversity


numberMutationLimit=1
boundsList=[(0,61),(66,128)] 
secondAppendPos=64
diversityPosNumber=6

with open(sys.argv[1],'rb') as f:
    newObjDict=cPickle.load(f)

dimerCounts={}

for eachRef in newObjDict:
    letterToAppend=eachRef[0]
    secondLetterToAppend=eachRef[secondAppendPos]
    for newObj in newObjDict[eachRef]:
        if isFullLength(newObj):
            if returnNumberMutations(newObj.alignmentInfo[1],boundsList)<numberMutationLimit:
                firstDiv=diversityStringReturn(newObj.alignmentInfo[1],boundsList,0)
                secondDiv=diversityStringReturn(newObj.alignmentInfo[1],boundsList,1)
                if len(firstDiv)==diversityPosNumber and len(secondDiv)==diversityPosNumber:
                    newFirst=letterToAppend+firstDiv
                    newSecond=secondLetterToAppend+secondDiv
                    if (newFirst,newSecond) in dimerCounts:
                        dimerCounts[(newFirst,newSecond)]+=newObj.count
                    else:
                        dimerCounts[(newFirst,newSecond)]=newObj.count

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.find('_trimmomatic')]

with open(sampleName+'_dimerDump_ver2.pckl','wb') as f:
    cPickle.dump(dimerCounts,f,protocol=cPickle.HIGHEST_PROTOCOL)
    

