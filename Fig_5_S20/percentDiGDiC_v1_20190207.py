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

readNumberToPlot=[]
totalDiGDiCToPlot=[]
totalReadCount=0
totalDiGCount=0
totalDiCCount=0
for counter, eachCluster in enumerate(clusterDict):
    for eachRef in clusterDict[eachCluster]:
        if clusterDict[eachCluster][eachRef].forAlignment==1: #subref collapsed dicts should have all subreferences set to forAlignment=1 anyway
            arrayIndices=np.array(clusterDict[eachCluster][eachRef].listSeqIndices)
            countsToUse=seqCount[arrayIndices]
            diGToUse=onlyDiGReq[arrayIndices]
            diCToUse=onlyDiCReq[arrayIndices]
            
            readNumberToPlot.append(np.sum(countsToUse))
            totalDiGDiCToPlot.append(np.sum(countsToUse[diGToUse])+np.sum(countsToUse[diCToUse]))

            totalReadCount+=np.sum(countsToUse)
            totalDiGCount+=np.sum(countsToUse[diGToUse])
            totalDiCCount+=np.sum(countsToUse[diCToUse])

totalReadCount=float(totalReadCount)
totalDiGCount=float(totalDiGCount)
totalDiCCount=float(totalDiCCount)

print 'Ratios of diG and diC for reporting in paper'
tReadCount=float(np.sum(seqCount))
tDiGCount=float(np.sum(seqCount[onlyDiGReq]))
tDiCCount=float(np.sum(seqCount[onlyDiCReq]))
print 'diG: '+str(tDiGCount/tReadCount)
print 'diC: '+str(tDiCCount/tReadCount)

print ''
print ''
print ''
print 'Ratios of diG and diC to use for simulation'
print 'diG: '+str(totalDiGCount/totalReadCount)
print 'diC: '+str(totalDiCCount/totalReadCount)


print 'Plotting CDF now...'

groupDict={}
for counter, eachEl in enumerate(readNumberToPlot):
    groupDict[counter]=float(totalDiGDiCToPlot[counter])/(totalDiGCount+totalDiCCount)


cdfX=[]
cdfY1=[]
cdfY2=[]

numberX=0
numberY1=0.0
numberY2=0.0
for eachEl in dictSorter(groupDict):
    numberX+=1
    numberY1+=float(readNumberToPlot[eachEl])/totalReadCount
    numberY2+=groupDict[eachEl]
    cdfX.append(numberX)
    cdfY1.append(numberY1)
    cdfY2.append(numberY2)

f, ax=plt.subplots(nrows=1,ncols=1)

ax.plot(cdfX,cdfY1,'k-')
ax.plot(cdfX,cdfY2,'g-')
ax.legend(['Total counts','GG,CC sequences'],loc=4)
ax.set_xscale('log')
ax.set_xlim(left=0.1)
ax.set_xlabel('Species Index')
ax.set_ylabel('CDF')
ax.set_title(sampleName)
#plt.show()
plt.savefig(sampleName+'_percentDiGDiCV1Code_20190207.pdf',dpi=300)
