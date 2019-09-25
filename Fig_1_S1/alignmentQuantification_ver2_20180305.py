import sys
import itertools

from moduleUsefulFunctions_20180215 import *

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

def returnAllTuplesClashes(dictionaryToUse):
    #returns list of tuples
    listToReturn=[]
    allPerms=itertools.permutations(dictionaryToUse.keys())
    for eachEl in allPerms:
        listToReturn.append(eachEl)
    return listToReturn

pcrDuplicates=0
toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
paramDistinguishWithinRef=4
#paramDistinguishAcrossRef=2*toleranceNumberMutations

with open(sys.argv[1],'rb') as f:
    endCollapsedDict=cPickle.load(f)

sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.find('_collapsedEndsClusterDict.pckl')]
with open('../'+sampleName+'_testDict.pckl','rb') as f:
    testDict=cPickle.load(f)

seqDict={}
for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][pcrDuplicates]


totalReads=float(sum(seqDict.values()))

refDict={}
for eachEl in endCollapsedDict:
    refDict[endCollapsedDict[eachEl][0]]=endCollapsedDict[eachEl][1]
numberStoreDict={}
for eachRef in refDict:
    numberStoreDict[eachRef]=np.zeros(3) #0th element is counts for reference, 1st element for rev comp of reference and 2nd element for ambiguous
numberStoreDict['Clashed']={}
numberStoreDict['Unaligned']=0.0

for eachEl in dictSorter(seqDict):
    #print eachEl+'\t'+str(seqDict[eachEl])
    #first figure out if a particular reference is well suited
    alignmentInfo=[]
    alignmentScore=[]
    for counter, eachRef in enumerate(dictSorter(refDict)):
        score=-1
        listRefs=[eachRef,revComp(eachRef)]
        for eachSeq in listRefs:
            alignOut=smithWaterman(eachSeq,eachEl,collapseInsertions=False)
            if alignOut[0]>score:
                score=alignOut[0]
                info=alignedObject(alignOut,eachRef,-1) #note eachRef is specified as refUsed, not eachSeq; artificial count of -1 passed
        alignmentInfo.append(info)
        alignmentScore.append(float(score)/min([len(eachEl),len(eachRef)])) #normalization by length is useful in downstream code below
    alignmentScore=np.array(alignmentScore)
    indicesSort=np.argsort(alignmentScore)
    bestAlignment=indicesSort[-1]
    numberMutations=returnNumberMutations(alignmentInfo[bestAlignment].alignmentInfo[1])
    endPositionBestAlignment=returnEndPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
    startPositionBestAlignment=returnStartPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
    maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(alignmentInfo[bestAlignment].refUsed),len(eachEl)]))
    totalNumberMutationsCalc=len(FivePrimeClipSequence(alignmentInfo[bestAlignment].alignmentInfo[1]))+len(ThreePrimeClipSequence(alignmentInfo[bestAlignment].alignmentInfo[1]))+numberMutations
    addToUnaligned=0
    goodRefFound=0
    if numberMutations<=np.ceil(toleranceNumberMutations*(endPositionBestAlignment-startPositionBestAlignment)):
        if totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength):
            if len(alignmentScore)>1:
                if alignmentScore[bestAlignment]>alignmentScore[indicesSort[-2]]:
                    goodRefFound=1
                else:
                    dictionaryClashingRefs={}
                    for eachAlignmentCounter, eachAlignmentScore in enumerate(alignmentScore):
                        if eachAlignmentScore==alignmentScore[bestAlignment]:
                            dictionaryClashingRefs[alignmentInfo[eachAlignmentCounter].refUsed]=0
                    listAllPossibleClashes=returnAllTuplesClashes(dictionaryClashingRefs)
                    matchingCombClashFound=0
                    for eachPossibleComb in listAllPossibleClashes:
                        if str(eachPossibleComb) in numberStoreDict['Clashed']:
                            matchingCombClashFound=1
                            numberStoreDict['Clashed'][str(eachPossibleComb)]+=float(seqDict[eachEl])
                            break
                    if matchingCombClashFound==0:
                        numberStoreDict['Clashed'][str(listAllPossibleClashes[0])]=float(seqDict[eachEl])             
            elif len(alignmentScore)==1:
                goodRefFound=1         
        else:
            addToUnaligned=1        
    else:
        addToUnaligned=1        
    if addToUnaligned==1:
        numberStoreDict['Unaligned']+=seqDict[eachEl]
    if goodRefFound==1:
        #now align again to get the best alignment to the sequence, rev comp or ambiguous
        goodRef=alignmentInfo[bestAlignment].refUsed
        listRefs=[goodRef,revComp(goodRef)]
        alignmentScore=[]
        for eachSeq in listRefs:
            alignOut=smithWaterman(eachSeq,eachEl,collapseInsertions=False)
            alignmentScore.append(alignOut[0]) #don't need length normalization as seq and revcomp have same length!
        if alignmentScore[0]>=alignmentScore[1]+paramDistinguishWithinRef:
            numberStoreDict[goodRef][0]+=seqDict[eachEl]
        elif alignmentScore[1]>=alignmentScore[0]+paramDistinguishWithinRef:
            numberStoreDict[goodRef][1]+=seqDict[eachEl]
        else:
            numberStoreDict[goodRef][2]+=seqDict[eachEl]                
    #print numberStoreDict
    #print ''
    #print ''

totalClashedPercentage=0.0    
for eachEl in numberStoreDict:
    if eachEl=='Clashed':
        for eachEl1 in numberStoreDict[eachEl]:
            numberStoreDict[eachEl][eachEl1]=(numberStoreDict[eachEl][eachEl1]/totalReads)*100
            totalClashedPercentage+=numberStoreDict[eachEl][eachEl1]
    else:
        numberStoreDict[eachEl]=(numberStoreDict[eachEl]/totalReads)*100
print numberStoreDict
print totalClashedPercentage


with open(sampleName+'_numberStoreDictVer2.pckl','wb') as f:
    cPickle.dump(numberStoreDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

