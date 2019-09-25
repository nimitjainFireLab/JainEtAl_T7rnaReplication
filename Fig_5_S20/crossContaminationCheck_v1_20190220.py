import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess

#ver16 denovoClustering code onwards definition of clusterObject
class clusterObject:
    __slots__=['reference','dictKmers','numberOfReads','forAlignment','listSeqIndices']    
    def __init__(self,reference,dictKmers,numberOfReads,forAlignment,listSeqIndices):
        self.reference=reference
        self.dictKmers=dictKmers
        self.numberOfReads=numberOfReads
        self.forAlignment=forAlignment
        self.listSeqIndices=listSeqIndices

def diNucAtEnds(seq,dotBrack,diNucSeq,distEnds):
    fiveIndex=dotBrack.find('(')
    seqFive=seq[:fiveIndex]
    seqFive=seqFive[:distEnds]
    
    threeIndex=dotBrack.rfind(')')+1
    seqThree=seq[threeIndex:]
    seqThree=seqThree[::-1]
    seqThree=seqThree[:distEnds]
    
    if diNucSeq in seqFive and diNucSeq in seqThree:
        return True
    else:
        return False

sortWithPCRDuplicates=0 #parameter to dictate whether sorting is happening with or without subtraction of PCR duplicates, 0 means that each barcode is only counted once for each insert sequence, 1 means that each barcode is counted as many times as it is observed even if the insert sequence is the same

testDictName=sys.argv[1]
with open(testDictName,'rb') as f:
    testDict=cPickle.load(f)

sampleName=testDictName
sampleName=sampleName[sampleName.rfind('/')+1:sampleName.rfind('_testDict.pckl')]

seqDict={}
for eachEntry in testDict:
    seqDict[eachEntry]=testDict[eachEntry][sortWithPCRDuplicates]

sortedDictionary=dictSorter(seqDict)

with open(sampleName+'_clusterDict_subrefCollapsed_v1.pckl','rb') as f:
    clusterDict=cPickle.load(f)

with open(sampleName+'_combinedList_structureMetricsPreCalc_v2.pckl','rb') as f:
    combinedList=cPickle.load(f)

[seqCount,oneHairpin,seqLength,diGReq,diCReq,dotBracketVector,isATstretch,hairpinScoreVector]=combinedList

#sortedDictionary list should be in same order as all these lists
seqCount=np.array(seqCount)
oneHairpin=np.array(oneHairpin)
seqLength=np.array(seqLength)
diGReq=np.array(diGReq)
diCReq=np.array(diCReq)
isATstretch=np.array(isATstretch)
hairpinScoreVector=np.array(hairpinScoreVector)
onlyDiGOrDiCReq=np.logical_xor(diGReq,diCReq)
onlyDiGReq=np.logical_and(diGReq,onlyDiGOrDiCReq)
onlyDiCReq=np.logical_and(diCReq,onlyDiGOrDiCReq)

distanceFromEnds=5 #in nucleotides
newOnlyDiGReq=[]
for counter, eachEl in enumerate(onlyDiGReq):
    if eachEl==False:
        newOnlyDiGReq.append(False)
    else:
        newOnlyDiGReq.append(diNucAtEnds(sortedDictionary[counter],dotBracketVector[counter],'GG',distanceFromEnds))
newOnlyDiGReq=np.array(newOnlyDiGReq)
onlyDiGReq=newOnlyDiGReq

newOnlyDiCReq=[]
for counter, eachEl in enumerate(onlyDiCReq):
    if eachEl==False:
        newOnlyDiCReq.append(False)
    else:
        newOnlyDiCReq.append(diNucAtEnds(sortedDictionary[counter],dotBracketVector[counter],'CC',distanceFromEnds))
newOnlyDiCReq=np.array(newOnlyDiCReq)
onlyDiCReq=newOnlyDiCReq

with open(sampleName+'_plotHairpinFunctionBestIndexVal_v4_20190220_allLists.pckl','rb') as f:
    [xToPlot1,yToPlot1,listAllRefs1,listAllClusters1,listdiGRef1]=cPickle.load(f)

resultsMatrix=[]
for counter, i in enumerate(range(1,17)):
    print i
    with open(str(i)+'_afterMS2removal_testDict.pckl','rb') as f:
        toCheckDict=cPickle.load(f)
    resultsMatrix.append([])
    for counter1, eachEl in enumerate(listdiGRef1):
        totalCount=0.0
        arrayIndices=np.array(clusterDict[listAllClusters1[counter1]][listAllRefs1[counter1]].listSeqIndices)
        for eachIndex in arrayIndices:
            if sortedDictionary[eachIndex] in toCheckDict:
                totalCount+=toCheckDict[sortedDictionary[eachIndex]][sortWithPCRDuplicates]
        resultsMatrix[counter].append(totalCount)
resultsMatrix=np.array(resultsMatrix)
print resultsMatrix

with open(sampleName+'_crossContaminationCheckResultsMatrix_v1_20190220.pckl','wb') as f:
    cPickle.dump(resultsMatrix,f,protocol=cPickle.HIGHEST_PROTOCOL)

sampleID=int(sampleName[:sampleName.find('_')])-1
histValues=[]
with open(sampleName+'_crossContaminationCheckByReference_v1_20190220.txt','w') as f:
    f.write('DiGRef\tBestIndexVal\tSubref\tClusterSeq\tTotalCount\tCrossCont')
    f.write('\n')
    for counter, eachEl in enumerate(listdiGRef1):
        f.write(eachEl)
        f.write('\t')
        f.write(str(yToPlot1[counter]))
        f.write('\t')
        f.write(listAllRefs1[counter])
        f.write('\t')
        f.write(listAllClusters1[counter])
        f.write('\t')
        f.write(str(xToPlot1[counter]))
        f.write('\t')
        f.write(str(float(resultsMatrix[sampleID,counter])/float(np.sum(resultsMatrix[:,counter]))))
        f.write('\n')
        histValues.append(float(resultsMatrix[sampleID,counter])/float(np.sum(resultsMatrix[:,counter])))

plt.hist(histValues,bins=np.arange(0.0,1.05,0.05))
plt.xlabel('1 - Cross contamination amount')
plt.ylabel('Number of references')
plt.title(sampleName)
#plt.show()
plt.savefig(sampleName+'_crossContaminationCheckHistogram_v1_20190220.pdf',dpi=300)


