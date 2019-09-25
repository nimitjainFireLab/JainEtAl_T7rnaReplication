import sys

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

newObjDictName=sys.argv[1]
reference=sys.argv[2]
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

f, axes = plt.subplots(len(refSequences), 1) #axes is a tuple (ax1,ax2,ax3....)
maxValue=0.0
for i in range(0,len(refSequences)):
    xToPlot=np.array(range(1,len(refSequences[i])+1))
    l1, =axes[i].plot(xToPlot,plotDict[i]['Transition']/fullLengthReads[i],'-o',label='Transition',color='red',linewidth=1.5,markersize=4.25,markeredgewidth=0.0)
    l2, =axes[i].plot(xToPlot,plotDict[i]['Transversion']/fullLengthReads[i],'-o',label='Transversion',color='orange',linewidth=1.5,markersize=4.25,markeredgewidth=0.0)
    l3, =axes[i].plot(xToPlot,plotDict[i]['Deletion']/fullLengthReads[i],'-o',label='Deletion',color='blue',linewidth=1.5,markersize=4.25,markeredgewidth=0.0)
    if max(plotDict[i]['Transition']/fullLengthReads[i])>maxValue:
        maxValue=max(plotDict[i]['Transition']/fullLengthReads[i])
    if max(plotDict[i]['Transversion']/fullLengthReads[i])>maxValue:
        maxValue=max(plotDict[i]['Transversion']/fullLengthReads[i])
    if max(plotDict[i]['Deletion']/fullLengthReads[i])>maxValue:
        maxValue=max(plotDict[i]['Deletion']/fullLengthReads[i])
    axes[i].set_xticks(xToPlot)
    #axes[i].set_yscale('symlog',linthreshy=0.01) #can also set basey parameter, default is 10
    #axes[i].set_yticks(np.concatenate((np.arange(0,0.01,0.001),np.arange(0.01,0.1,0.01),np.arange(0.1,1.1,0.1))))
    #yticks=list(np.arange(0,0.01,0.001))+[0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
    axes[i].set_xticklabels(list(refSequences[i].replace('T','U')))    
    axes[i].set_xlim(left=0.0,right=len(refSequences[i])+1)
    f.legend([l1,l2,l3],['Transition','Transversion','Deletion'],loc=8,ncol=3)

maxValue=(np.ceil(maxValue*1000))/1000
for eachEl in axes:
    yticks=np.arange(0,1.2*maxValue,(maxValue/5.0))    
    eachEl.set_yticks(yticks)
    eachEl.set_yticklabels(yticks)
    eachEl.set_ylim(bottom=0.0,top=1.1*maxValue)
    eachEl.grid(which='major',ls='-',lw=0.3,alpha=0.1)    
    
f.set_size_inches(12.5,12)
sampleName=sys.argv[3]
print 'Ref '+sampleName+':\t'+str(totalReads)
print 'Ref '+sampleName+':\t'+str(fullLengthReads)
print ''
print ''
#plt.show()
plt.savefig('Ref_'+sampleName+'_mutationTracer_ver5.pdf',dpi=300)

