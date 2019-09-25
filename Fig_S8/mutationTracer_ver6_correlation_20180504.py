import sys
from scipy.stats import pearsonr
import scikits.bootstrap as boot

from moduleUsefulFunctions_20180215 import *
plt.close()

def deletionReversal(seq):
    posReassignedDict={} #dictionary
    for counter, eachEl in enumerate(seq):
        posReassignedDict[counter+1]=counter+1
        if counter<len(seq)-1:
            for counter1, eachEl1 in enumerate(seq[counter+1:]):
                if eachEl==eachEl1:
                    posReassignedDict[counter+1]+=1
                else:
                    break
    return posReassignedDict

def deletionStraight(seq):
    posDict={}
    for counter, eachEl in enumerate(seq):
        posDict[counter+1]=counter+1
    return posDict

def bootstrapFunction(data):
    #data is an indices list
    a=pearsonr(data1[data],data2[data])
    return a[0]

#Ref 2
newObjDictName='tf2_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl'
reference='GGAAAATATACATATTGAAGGTGTGTATGTATATTTGTATATTCACAAAAATATACATACACACCTTCAATATGTATATTATTGG'

#Y2 RNA 
#newObjDictName='5_trimmomaticP_thresholdTestDict_1_seqErrorCorrect_1_newObjDict.pckl'
#reference='GGAAAATTTCAAGATCAGGGCTTGAAATTTTACAAAATTTCAAGCCCTGATCTTGAAATTTTGG'
refSequences=[reference,revComp(reference)]

with open(newObjDictName,'rb') as f:
    newObjDict=cPickle.load(f)
    
PCRduplicates=0 #0 means PCR duplicates are not there, 1 means PCR duplicates potentially there in counting
numberMutationsTolerate=1
endParameter=2
plotDict=[]
fullLengthReads=np.zeros((len(refSequences),1))
totalReads=np.zeros((len(refSequences),1))
transitionDict={'G':'A','A':'G','C':'T','T':'C'}
for i in range(0,len(refSequences)):
    lengthRefSeq=len(refSequences[i])
    dictToAppend={}
    #dictToAppend['Mismatch']=np.zeros(lengthRefSeq)
    dictToAppend['Transition']=np.zeros(lengthRefSeq)
    dictToAppend['Transversion']=np.zeros(lengthRefSeq)
    dictToAppend['Deletion']=np.zeros(lengthRefSeq)
    plotDict.append(dictToAppend)

for referenceIndex, eachRef in enumerate(refSequences):    
    if referenceIndex==0:
        deletionPosDict=deletionStraight(eachRef)
    elif referenceIndex==1:
        deletionPosDict=deletionReversal(eachRef)
    for newObj in newObjDict[eachRef]:
        totalReads[referenceIndex]+=newObj.count[PCRduplicates]
        numberMutations=0
        for eachEl in newObj.alignmentInfo[1]:
            if eachEl.typeOfEvent=='Start':
                startPosition=eachEl.position
            elif eachEl.typeOfEvent=='End':
                endPosition=eachEl.position
        for eachEl in newObj.alignmentInfo[1]:
            if eachEl.typeOfEvent=='Insertion':
                if eachEl.position<=(startPosition+endParameter) or eachEl.position>=(endPosition-endParameter): #only count internal insertions
                    pass
                else:
                    numberMutations+=len(eachEl.notes) #by doing this, doesn't matter if collapseInsertions=True or False in smithwaterman call
            elif eachEl.typeOfEvent=='Mismatch' or eachEl.typeOfEvent=='Deletion':
                numberMutations+=1
        if startPosition<=endParameter and endPosition>=(len(newObj.refUsed)-endParameter): 
            fullLengthReads[referenceIndex]+=newObj.count[PCRduplicates]
            if numberMutations<=numberMutationsTolerate:
                for eachEl in newObj.alignmentInfo[1]: #seven different event types--start, insertion, end, match, diversity, mismatch, deletion
                    if eachEl.typeOfEvent=='Deletion':
                        plotDict[referenceIndex][eachEl.typeOfEvent][deletionPosDict[eachEl.position]-1]+=newObj.count[PCRduplicates]
                    elif eachEl.typeOfEvent=='Mismatch':
                        if transitionDict[eachRef[eachEl.position-1]]==eachEl.notes:
                            plotDict[referenceIndex]['Transition'][eachEl.position-1]+=newObj.count[PCRduplicates]
                        else:
                            plotDict[referenceIndex]['Transversion'][eachEl.position-1]+=newObj.count[PCRduplicates]

f,ax=plt.subplots(nrows=1,ncols=2)
thresholdForPearson=-1
datasets=[]
for i in range(0,len(refSequences)):
    if i==0:
        a1=plotDict[i]['Transition']/fullLengthReads[i]
        a2=plotDict[i]['Transversion']/fullLengthReads[i]
        a3=plotDict[i]['Deletion']/fullLengthReads[i]
    elif i==1:
        a1=plotDict[i]['Transition'][::-1]/fullLengthReads[i]
        a2=plotDict[i]['Transversion'][::-1]/fullLengthReads[i]
        a3=plotDict[i]['Deletion'][::-1]/fullLengthReads[i]
    datasets.append(np.concatenate((a1,a2,a3)))

finalDatasets=[[],[]]
for counter, eachEl in enumerate(datasets[0]):
    if eachEl<=thresholdForPearson and datasets[1][counter]<=thresholdForPearson:
        print 'Yes'
    else:
        finalDatasets[0].append(eachEl)
        finalDatasets[1].append(datasets[1][counter])

ax[0].plot(finalDatasets[0],finalDatasets[1],'ko')
print pearsonr(finalDatasets[0],finalDatasets[1])
#plt.show()

#Bootstrap for correlation coefficient
numberTrials=100000
correlationValuesBoot=[]
data1=np.array(finalDatasets[0])
data2=np.array(finalDatasets[1])
for i in range(0,numberTrials):
    indexComb=np.random.choice(len(finalDatasets[0]),size=len(finalDatasets[0]),replace=True) #bootstrapping should be done, sampling with replacement (hence replace=True)--this is bootstrapping the index
    a=pearsonr(data1[indexComb],data2[indexComb])
    correlationValuesBoot.append(a[0])

print 'Percentile intervals'
print np.percentile(correlationValuesBoot,2.5)
print np.percentile(correlationValuesBoot,97.5)

print 'Bca intervals'
confidenceIntervals=boot.ci(np.arange(len(data1)), statfunction=bootstrapFunction, alpha=0.05, n_samples=numberTrials, method='bca', output='lowhigh')
#first row of confidenceIntervals is lowPercentile
#second row of confidenceIntervals is highPercentile
print confidenceIntervals


ax[1].hist(correlationValuesBoot)
plt.show()
