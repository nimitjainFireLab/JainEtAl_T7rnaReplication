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

pcrDuplicates=0
toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object
paramDistinguishWithinRef=4

sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.find('_testDict.pckl')]
metricsName=sampleName+'_alignmentmetrics.pckl'
alignPickleName=sampleName+'_newObjDict.pckl'
unalignedPickleName=sampleName+'_unalignedDict.pckl'

alignmentMetrics={}
alignmentMetrics['pcrDuplicates']=pcrDuplicates
alignmentMetrics['toleranceNumberMutations']=toleranceNumberMutations
alignmentMetrics['ntLengthUpdateFullLength']=ntLengthUpdateFullLength
alignmentMetrics['paramDistinguishWithinRef']=paramDistinguishWithinRef
alignmentMetrics['testDict']=sys.argv[1]
alignmentMetrics['codeUsed']=sys.argv[0]

newObjDict={} #keys are references, values are lists of alignments 
unalignedDict={}

with open(sys.argv[1],'rb') as f:
    testDict=cPickle.load(f)

totalReads=np.sum(testDict.values(),axis=0)[pcrDuplicates]

gStrandRefSequence='GGAAAATTAGNAGGTNGTAGTNCTAAATTTTCCAAAATTTAGNACTACNACCTNCTAATTTTGG'
refList=[gStrandRefSequence+gStrandRefSequence,gStrandRefSequence+revComp(gStrandRefSequence),revComp(gStrandRefSequence)+gStrandRefSequence,revComp(gStrandRefSequence)+revComp(gStrandRefSequence)]
numberStoreDict={}
for eachRef in refList:
    numberStoreDict[eachRef]=0.0 
numberStoreDict['Clashed']=0.0
numberStoreDict['Unaligned']=0.0

for eachEl in testDict:
    #print eachEl+'\t'+str(testDict[eachEl])
    #first figure out if a particular reference is well suited
    alignmentInfo=[]
    alignmentScore=[]
    for eachSeq in refList:
        alignOut=smithWaterman(eachSeq,eachEl,collapseInsertions=False)
        alignmentScore.append(alignOut[0]) #don't need length normalization as all references have same length!
        alignmentInfo.append(alignedObject(alignOut,eachSeq,testDict[eachEl]))
    alignmentScore=np.array(alignmentScore)
    indicesSort=np.argsort(alignmentScore)
    bestAlignment=indicesSort[-1]
    numberMutations=returnNumberMutations(alignmentInfo[bestAlignment].alignmentInfo[1])
    endPositionBestAlignment=returnEndPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
    startPositionBestAlignment=returnStartPosition(alignmentInfo[bestAlignment].alignmentInfo[1])
    maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(alignmentInfo[bestAlignment].refUsed),len(eachEl)]))
    totalNumberMutationsCalc=len(FivePrimeClipSequence(alignmentInfo[bestAlignment].alignmentInfo[1]))+len(ThreePrimeClipSequence(alignmentInfo[bestAlignment].alignmentInfo[1]))+numberMutations
    addToUnaligned=0
    if numberMutations<=np.ceil(toleranceNumberMutations*(endPositionBestAlignment-startPositionBestAlignment)):
        if totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength):
                if alignmentScore[bestAlignment]>=alignmentScore[indicesSort[-2]]+paramDistinguishWithinRef:                    
                    newObj=alignmentInfo[bestAlignment]
                    numberStoreDict[newObj.refUsed]+=float(testDict[eachEl][pcrDuplicates])
                    if newObj.refUsed in newObjDict:
                        newObjDict[newObj.refUsed].append(newObj)
                    else:
                        newObjDict[newObj.refUsed]=[newObj]                   
                else:
                    numberStoreDict['Clashed']+=float(testDict[eachEl][pcrDuplicates])
        else:
            addToUnaligned=1        
    else:
        addToUnaligned=1        
    if addToUnaligned==1:
        numberStoreDict['Unaligned']+=float(testDict[eachEl][pcrDuplicates])
        unalignedDict[eachEl]=testDict[eachEl]    
            
for eachEl in numberStoreDict:
    numberStoreDict[eachEl]=(numberStoreDict[eachEl]/totalReads)*100
print numberStoreDict

with open(sampleName+'_numberStoreDict.pckl','wb') as f:
    cPickle.dump(numberStoreDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(metricsName,'wb') as f:
    cPickle.dump(alignmentMetrics,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(alignPickleName,'wb') as f:
    cPickle.dump(newObjDict,f,protocol=cPickle.HIGHEST_PROTOCOL)

with open(unalignedPickleName,'wb') as f:
    cPickle.dump(unalignedDict,f,protocol=cPickle.HIGHEST_PROTOCOL)
