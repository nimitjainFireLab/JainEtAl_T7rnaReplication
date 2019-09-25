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
            
def mutationTupleReturn(alignmentInfoList, boundsList):
    listMutations=[]
    startPos=returnStartPosition(alignmentInfoList)
    endPos=returnEndPosition(alignmentInfoList)
    for eachEl in alignmentInfoList:
        if checkWithinBounds(eachEl.position,boundsList):
            if eachEl.typeOfEvent=='Insertion':
                if eachEl.position!=startPos and eachEl.position!=endPos:
                    listMutations.append('I'+str(eachEl.position)+str(eachEl.notes))
            elif eachEl.typeOfEvent=='Deletion':
                listMutations.append('D'+str(eachEl.position))
            elif eachEl.typeOfEvent=='Mismatch':
                listMutations.append('M'+str(eachEl.position)+str(eachEl.notes))
            elif eachEl.typeOfEvent=='Diversity':
                listMutations.append('N'+str(eachEl.position)+str(eachEl.notes))
    return tuple(listMutations)   

def diversityStringReturn(alignmentInfoList,boundsList,k):
    #Return diversity string
    strDiversity=''
    for eachEl in alignmentInfoList:
        if eachEl.typeOfEvent=='Diversity':
            if checkWithinBounds(eachEl.position,[boundsList[k]]):
                strDiversity+='N'+str(eachEl.position-k*64)+str(eachEl.notes)
    return strDiversity
            
with open(sys.argv[1],'rb') as f:
    newObjDict=cPickle.load(f)

boundsList=[(0,61),(66,128)] #count mutations within these boundaries
PCRduplicates=0

countsByMutationDict={}
for refUsed in newObjDict:
    for newObj in newObjDict[refUsed]:
        if isFullLength(newObj):
            mutationTuple=mutationTupleReturn(newObj.alignmentInfo[1],boundsList)
            diversityString=diversityStringReturn(newObj.alignmentInfo[1],boundsList,0)
            diversityStringSecond=diversityStringReturn(newObj.alignmentInfo[1],boundsList,1)
            if diversityString==diversityStringSecond:
                if refUsed in countsByMutationDict:
                    if diversityString in countsByMutationDict[refUsed]:
                        if mutationTuple in countsByMutationDict[refUsed][diversityString]:
                            countsByMutationDict[refUsed][diversityString][mutationTuple]+=newObj.count[PCRduplicates]
                        else:
                            countsByMutationDict[refUsed][diversityString][mutationTuple]=newObj.count[PCRduplicates]            
                    else:
                        countsByMutationDict[refUsed][diversityString]={}
                        countsByMutationDict[refUsed][diversityString][mutationTuple]=newObj.count[PCRduplicates]          
                else:            
                    countsByMutationDict[refUsed]={}
                    countsByMutationDict[refUsed][diversityString]={}
                    countsByMutationDict[refUsed][diversityString][mutationTuple]=newObj.count[PCRduplicates]

sampleName=sys.argv[1]
sampleName=sampleName[:sampleName.rfind('_newObjDict.pckl')]

with open(sampleName+'_extractAlignments_20180413_ver2_out.txt','w') as f:
    for eachRef in countsByMutationDict:
        f.write('REFERENCE\n')
        f.write(eachRef)
        f.write('\n')
        for eachDiv in countsByMutationDict[eachRef]:
            f.write(eachDiv)
            f.write('\n')
            for eachEntry in dictSorter(countsByMutationDict[eachRef][eachDiv]):
                f.write(str(eachEntry)+'\t'+str(countsByMutationDict[eachRef][eachDiv][eachEntry]))
                f.write('\n')
        f.write('\n')
        f.write('\n')
