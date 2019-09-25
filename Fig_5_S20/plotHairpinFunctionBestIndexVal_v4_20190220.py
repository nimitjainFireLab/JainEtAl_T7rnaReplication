import cPickle
import sys
from moduleUsefulFunctions_20180215 import *
from random import randint
import subprocess
import forgi
import forgi.graph.bulge_graph as fgb

def returnElNotation(dotBracket):
    bg=fgb.BulgeGraph(dotbracket_str=dotBracket)
    return bg.to_element_string()

def returnLeftRightLengths(elementRep):
    numberHairpin=float(elementRep.count('h'))
    fiveIndex=elementRep.find('h')
    threeIndex=elementRep.rfind('h')+1
    
    leftS=0.0
    leftTotal=0.0
    for eachEl in elementRep[:fiveIndex]:
        leftTotal+=1.0
        if eachEl=='s':
            leftS+=1.0
    leftTotal+=numberHairpin/2.0   

    rightS=0.0
    rightTotal=0.0
    for eachEl in elementRep[threeIndex:]:
        rightTotal+=1.0
        if eachEl=='s':
            rightS+=1.0
    rightTotal+=numberHairpin/2.0   

    return [leftS/leftTotal,rightS/rightTotal]
    

def returnDotBracket(sequence,sampleLabel):
    with open('outputRNAfold_plotHairpinFunctionBestIndexVal_v4_20190220_'+sampleLabel+'.txt','w') as f:
        processCall=subprocess.Popen(['RNAfold','--noPS','--noGU'],stdout=f,stdin=subprocess.PIPE)
        processCall.communicate(input=sequence)

sampleName=sys.argv[1]

with open(sampleName+'_bestIndexValCalculatorV2Code_20190218_allLists.pckl','rb') as f:
    [xToPlot1,yToPlot1,listAllRefs1,listAllClusters1,listdiGRef1]=cPickle.load(f)

newRefList=[]
newBestIndexValList=[]
newTotalCount=[]
newClusterDictRef=[]
newClusterDictCluster=[]
for counter, eachEl in enumerate(yToPlot1):
    if eachEl!=0:
        newRefList.append(listdiGRef1[counter])
        newBestIndexValList.append(eachEl)
        newTotalCount.append(xToPlot1[counter])
        newClusterDictRef.append(listAllRefs1[counter])
        newClusterDictCluster.append(listAllClusters1[counter])
        
fastaString=''
for counter, eachEl in enumerate(newRefList):
    line1='>'+str(counter)+'\n'
    line2=eachEl+'\n'
    fastaString+=line1
    fastaString+=line2         
        
returnDotBracket(fastaString,sampleName)

newDotBrack=[]
with open('outputRNAfold_plotHairpinFunctionBestIndexVal_v4_20190220_'+sampleName+'.txt','rU') as f:
    for counter, eachLine in enumerate(f):
        if counter%3==2:
            dotBracket=eachLine.split()[0]
            newDotBrack.append(dotBracket)
seqLengths=[]
repeatLengths=[]
for counter, eachEl in enumerate(newDotBrack):
    seqLengths.append(float(len(eachEl)))
    repeatLengths.append(float(eachEl.count('(')))

seqLengths=np.array(seqLengths)
repeatLengths=np.array(repeatLengths)
ratioRepeatLength=repeatLengths/seqLengths
xToPlot=[]
yToPlot=[]
colorToPlot=[]
for counter, eachEl in enumerate(newRefList):
    [xcoord,ycoord]=returnLeftRightLengths(returnElNotation(newDotBrack[counter]))
    xToPlot.append(xcoord)
    yToPlot.append(ycoord)
    if ratioRepeatLength[counter]<0.3:
        colorToPlot.append(np.array([0.0,0.0,0.0]))
    elif ratioRepeatLength[counter]<0.375:
        colorToPlot.append(np.array([1.0,0.0,0.0]))
    else:
        colorToPlot.append(np.array([0.0,0.0,1.0]))

xToPlot=np.array(xToPlot)
yToPlot=np.array(yToPlot)
ratioXY=yToPlot/xToPlot
boxList=[]
xTicks=np.arange(0.0,0.9,0.05)
medianList=[]
for i in xTicks:
    logicalArray=np.logical_and(xToPlot>i,yToPlot>i)
    boxList.append(ratioXY[logicalArray])
    medianList.append(np.median(ratioXY[logicalArray]))

hairpinThreshold=0.75
f, ax=plt.subplots(nrows=3,ncols=2)
f.set_figheight(16)
f.set_figwidth(18)
ax[0][0].scatter(xToPlot,yToPlot,c=np.array(colorToPlot),edgecolors='none',s=15)
ax[0][0].plot(np.arange(0.0,1.0,0.1),np.arange(0.0,1.0,0.1),'k--')
ax[0][0].plot(np.arange(0.0,1.0,0.1),[hairpinThreshold]*len(np.arange(0.0,1.0,0.1)),'k--')
ax[0][0].plot([hairpinThreshold]*len(np.arange(0.0,1.0,0.1)),np.arange(0.0,1.0,0.1),'k--')
ax[0][0].set_xlabel('Left hairpin ratio')
ax[0][0].set_ylabel('Right hairpin ratio')
ax[0][0].set_title(sampleName)
ax[0][1].boxplot(boxList)
ax[0][1].set_xticks(np.arange(1,len(xTicks)+1))
ax[0][1].set_xticklabels(np.round(xTicks,2))
ax[0][1].set_xlabel('Hairpin threshold')
ax[0][1].set_ylabel('Right hairpin ratio/Left hairpin ratio')
ax[1][0].plot(xTicks,medianList,'bo')
ax[1][0].plot([hairpinThreshold]*len(np.arange(min(medianList),max(medianList),0.01)),np.arange(min(medianList),max(medianList),0.01),'k--')
ax[1][0].set_xlabel('Hairpin threshold')
ax[1][0].set_ylabel('Median{Right hairpin ratio/Left hairpin ratio}')
h=ax[1][1].hist2d(xToPlot,yToPlot,bins=np.arange(0.0,1.0,0.05))
ax[1][1].plot(np.arange(0.0,1.0,0.1),[hairpinThreshold]*len(np.arange(0.0,1.0,0.1)),'k--')
ax[1][1].plot([hairpinThreshold]*len(np.arange(0.0,1.0,0.1)),np.arange(0.0,1.0,0.1),'k--')
plt.colorbar(h[3],ax=ax[1][1])
ax[1][1].set_xlabel('Left hairpin ratio')
ax[1][1].set_ylabel('Right hairpin ratio')



lengthToPlot=seqLengths[np.logical_and(xToPlot>hairpinThreshold,yToPlot>hairpinThreshold)]
ratiosToPlot=ratioRepeatLength[np.logical_and(xToPlot>hairpinThreshold,yToPlot>hairpinThreshold)]

ax[2][0].plot(lengthToPlot,ratiosToPlot,'bo') #IMPORTANT: THE POINT OF THIS PLOT IS TO ASSESS WHETHER THERE IS ANY LENGTH BIAS OF SETTING A HARD HAIRPINTHRESHOLD--IN GENERAL, I SEE THAT THE RATIO INCREASES WITH LENGTH WHICH IS GREAT, SUGGESTING THAT THERE IS MINIMAL LENGTH BIAS--THE IDEA IS THAT A THRESHOLD OF 0.75 MEANS VERY DIFFERENT THINGS AT 50 BP AS OPPOSED TO 100 BP OF RNA (IN TERMS OF NUMBER OF FREE BASES ETC.) AND THE EXPECTATION IS THAT SETTING A HARD HAIRPINTHRESHOLD IS PERHAPS TOO LOOSE A REQUIREMENT FOR LONGER SEQUENCES--BUT IN MY CASE, THE RATIO INCREASES WITH LENGTH WHICH IS GOOD!
ax[2][0].set_xlabel('Sequence length')
ax[2][0].set_ylabel('Hairpin stem length/Sequence length')
ax[2][1].hist(2*ratioRepeatLength,bins=np.arange(0.0,1.05,0.05)) #multiplication by 2 for direct comparison with the 2d histogram
ax[2][1].set_xlabel('Basepaired bases/Sequence length')
ax[2][1].set_ylabel('Count')

plt.savefig(sampleName+'_plotHairpinFunctionBestIndexVal_v4_20190220_plot1.pdf',dpi=300)

newRefList2=[]
newBestIndexValList2=[]
newTotalCount2=[]
newClusterDictRef2=[]
newClusterDictCluster2=[]


for counter, eachEl in enumerate(xToPlot):
    if eachEl>hairpinThreshold and yToPlot[counter]>hairpinThreshold:
        newRefList2.append(newRefList[counter])
        newBestIndexValList2.append(newBestIndexValList[counter])
        newTotalCount2.append(newTotalCount[counter])
        newClusterDictRef2.append(newClusterDictRef[counter])
        newClusterDictCluster2.append(newClusterDictCluster[counter])



f, ax=plt.subplots(nrows=1,ncols=2)
f.set_figheight(8)
f.set_figwidth(16)
ax[0].plot(newTotalCount,newBestIndexValList,'ko',markersize=2.0) #note several points may be overlapping here, especially for low cluster counts, this is ok for the purposes here
ax[0].plot(10.0*np.arange(1,max(newBestIndexValList)+1),np.arange(1,max(newBestIndexValList)+1),'b--')
ax[0].plot(100.0*np.arange(1,max(newBestIndexValList)+1),np.arange(1,max(newBestIndexValList)+1),'r--')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim(left=0.1)
ax[0].set_ylim(bottom=0.1)
ax[0].set_title(sampleName)
ax[0].set_xlabel('Total counts of reference-TotalCount')
ax[0].set_ylabel('Min. counts for complementary diG and diC seq of reference-bestIndexVal')
        
ax[1].plot(newTotalCount2,newBestIndexValList2,'ko',markersize=2.0) #note several points may be overlapping here, especially for low cluster counts, this is ok for the purposes here
ax[1].plot(10.0*np.arange(1,max(newBestIndexValList2)+1),np.arange(1,max(newBestIndexValList2)+1),'b--')
ax[1].plot(100.0*np.arange(1,max(newBestIndexValList2)+1),np.arange(1,max(newBestIndexValList2)+1),'r--')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim(left=0.1)
ax[1].set_ylim(bottom=0.1)
ax[1].set_title(sampleName+', after hairpin thresholding')
ax[1].set_xlabel('Total counts of reference-TotalCount')
ax[1].set_ylabel('Min. counts for complementary diG and diC seq of reference-bestIndexVal')    

plt.savefig(sampleName+'_plotHairpinFunctionBestIndexVal_v4_20190220_plot2.pdf',dpi=300)

with open(sampleName+'_plotHairpinFunctionBestIndexVal_v4_20190220_allLists.pckl','wb') as f:
    cPickle.dump([newTotalCount2,newBestIndexValList2,newClusterDictRef2,newClusterDictCluster2,newRefList2],f,protocol=cPickle.HIGHEST_PROTOCOL)
    

