import sys

from moduleUsefulFunctions_20180215 import *

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


toleranceNumberMutations=0.1 #for subclusters, this is number of mutations per nucleotide, i.e. for a value of 0.1, 1 mutation is allowed every 10 bases on average, another good acceptable value is 0.06
ntLengthUpdateFullLength=6 #how much should a new sequence be longer by to update the fullLengthSequence attribute for each cluster object

with open(sys.argv[1],'rb') as f:
    clusterDict=cPickle.load(f)

newClusterDict={}
for eachCluster in clusterDict:
    if eachCluster=='Clashed' or eachCluster=='Unaligned':
        pass
    else:
        for eachRef in clusterDict[eachCluster]:
            if eachRef!='Unaligned':
                newClusterDict[eachRef]=clusterDict[eachCluster][eachRef].numberOfReads

finalClusterRefList=[]
associatedCountsList=[]
for counter, eachRef in enumerate(dictSorter(newClusterDict)):    
    if counter==0:
        finalClusterRefList.append(eachRef)
        associatedCountsList.append(newClusterDict[eachRef])
    else:
        #align both ref and rev comp to every previous reference and check whether alignment good or not
        canEliminate=0
        for eachRef1 in finalClusterRefList:
            score=-1
            listRefs=[eachRef1,revComp(eachRef1)]
            for eachSeq in listRefs:
                alignOut=smithWaterman(eachSeq,eachRef,collapseInsertions=False)
                if alignOut[0]>score:
                    score=alignOut[0]
                    info=alignOut[1]
            numberMutations=returnNumberMutations(info)
            maxNumberMutationsAllowed=np.ceil(toleranceNumberMutations*min([len(eachRef1),len(eachRef)]))
            totalNumberMutationsCalc=len(FivePrimeClipSequence(info))+len(ThreePrimeClipSequence(info))+numberMutations
            if numberMutations<=np.ceil(toleranceNumberMutations*(returnEndPosition(info)-returnStartPosition(info))):
                if totalNumberMutationsCalc<=(maxNumberMutationsAllowed+ntLengthUpdateFullLength):
                    canEliminate=1
                    break                
        if canEliminate==0:
            finalClusterRefList.append(eachRef)
            associatedCountsList.append(newClusterDict[eachRef])
    
printClusterDict(clusterDict)

finalClusterRefDict={}
for counter, eachEl in enumerate(finalClusterRefList):
    finalClusterRefDict[eachEl]=associatedCountsList[counter]

print finalClusterRefDict
        
sampleName=sys.argv[1]
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_clusterDict.pckl')]

with open(sampleName+'_smithWatermanCollapsedClusterDict.pckl','wb') as f:
    cPickle.dump(finalClusterRefDict,f,protocol=cPickle.HIGHEST_PROTOCOL)
